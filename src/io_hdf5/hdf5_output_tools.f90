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

MODULE MOD_HDF5_Output_Tools
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

#if USE_QDS_DG
INTERFACE WriteQDSToHDF5
  MODULE PROCEDURE WriteQDSToHDF5
END INTERFACE
PUBLIC :: WriteQDSToHDF5
#endif /*USE_QDS_DG*/

#if !(USE_HDG)
PUBLIC :: WritePMLzetaGlobalToHDF5
#endif /*USE_HDG*/

#ifdef PARTICLES
PUBLIC :: WriteIMDStateToHDF5
#endif /*PARTICLES*/

PUBLIC :: WriteDielectricGlobalToHDF5,WriteBFieldToHDF5,WriteBFieldAnalyticToHDF5
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteDielectricGlobalToHDF5()
!===================================================================================================================================
! write DielectricGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars ,ONLY: DielectricGlobal,DielectricEps
USE MOD_Dielectric_Vars ,ONLY: DielectricMu,isDielectricElem,ElemToDielectric
USE MOD_Mesh_Vars       ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars    ,ONLY: ProjectName
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=2
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
#if USE_MPI
REAL                :: StartT,EndT
#endif
REAL                :: OutputTime!,FutureTime
INTEGER             :: iElem
!===================================================================================================================================
! create global Eps field for parallel output of Eps distribution
ALLOCATE(DielectricGlobal(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='DielectricEpsGlobal'
StrVarNames(2)='DielectricMuGlobal'
DielectricGlobal=0.
DO iElem=1,PP_nElems
  IF(isDielectricElem(iElem))THEN
    DielectricGlobal(1,:,:,:,iElem)=DielectricEps(:,:,:,ElemToDielectric(iElem))
    DielectricGlobal(2,:,:,:,iElem)=DielectricMu( :,:,:,ElemToDielectric(iElem))
  END IF
END DO!iElem
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DielectricGlobal TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#endif
END IF
OutputTime=0.0
!FutureTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DielectricGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('DielectricGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)!,FutureTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
CALL WriteAttributeToHDF5(File_ID,'VarNamesDielectricGlobal',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        PP_N            => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5                                             , &
                          nValGlobal =(/N_variables , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=DielectricGlobal)
END ASSOCIATE
#if USE_MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif
SDEALLOCATE(DielectricGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteDielectricGlobalToHDF5



#if !(USE_HDG)
SUBROUTINE WritePMLzetaGlobalToHDF5()
!===================================================================================================================================
! write PMLzetaGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PML_Vars     ,ONLY: PMLzetaGlobal,PMLzeta0,PMLzeta,isPMLElem,ElemToPML
USE MOD_Mesh_Vars    ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars ,ONLY: ProjectName
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: N_variables
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
#if USE_MPI
REAL                :: StartT,EndT
#endif
REAL                :: OutputTime!,FutureTime
INTEGER             :: iElem
!===================================================================================================================================
N_variables=3
! create global zeta field for parallel output of zeta distribution
ALLOCATE(PMLzetaGlobal(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='PMLzetaGlobalX'
StrVarNames(2)='PMLzetaGlobalY'
StrVarNames(3)='PMLzetaGlobalZ'
PMLzetaGlobal=0.
DO iElem=1,PP_nElems
  IF(isPMLElem(iElem))THEN
    IF(ALMOSTZERO(PMLzeta0))THEN
      PMLzetaGlobal(:,:,:,:,iElem)=0.0
    ELSE
      PMLzetaGlobal(:,:,:,:,iElem)=PMLzeta(:,:,:,:,ElemToPML(iElem))/PMLzeta0
    END IF
  END IF
END DO!iElem
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE PMLZetaGlobal TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#endif
END IF
OutputTime=0.0
!FutureTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PMLZetaGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('PMLZetaGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)!,FutureTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
CALL WriteAttributeToHDF5(File_ID,'VarNamesPMLzetaGlobal',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        PP_N            => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5                                             , &
                          nValGlobal =(/N_variables , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=PMLzetaGlobal)
END ASSOCIATE
#if USE_MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif
SDEALLOCATE(PMLzetaGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WritePMLzetaGlobalToHDF5
#endif /*USE_HDG*/


#if USE_QDS_DG
SUBROUTINE WriteQDSToHDF5(OutputTime,PreviousTime)
!===================================================================================================================================
! write QDS field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars ,ONLY: ProjectName
USE MOD_io_HDF5
USE MOD_QDS_DG_Vars  ,ONLY: nQDSElems,QDSSpeciesMass,QDSMacroValues
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: PreviousTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: N_variables
CHARACTER(LEN=255)  :: StrVarNames(1:6)
CHARACTER(LEN=255)  :: FileName
#if USE_MPI
REAL                :: StartT,EndT
#endif
INTEGER             :: iElem,j,k,l
REAL                :: Utemp(1:6,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems)
!===================================================================================================================================
N_variables=6
! create global Eps field for parallel output of Eps distribution
StrVarNames(1) = 'Density'
StrVarNames(2) = 'VeloX'
StrVarNames(3) = 'VeloY'
StrVarNames(4) = 'VeloZ'
StrVarNames(5) = 'Energy'
StrVarNames(6) = 'Temperature'
Utemp=0.
DO iElem =1, nQDSElems
  DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
    IF (QDSMacroValues(1,l,k,j,iElem).GT.0.0) THEN
!      Utemp(1,l,k,j,iElem) = QDSMacroValues(1,l,k,j,iElem)/(Species(QDS_Species)%MassIC*wGP(l)*wGP(k)*wGP(j))*sJ(l,k,j,iElem)
      Utemp(1,l,k,j,iElem) = QDSMacroValues(1,l,k,j,iElem)/QDSSpeciesMass
      IF (Utemp(1,l,k,j,iElem).LT.0.0) then
        print*, 'Utemp(1,l,k,j,iElem).LT.0.0'
        print*, Utemp(1,l,k,j,iElem),iElem, l,k,j, QDSMacroValues(1,l,k,j,iElem)
        print*,"Press ENTER to continue"
        read*
      END IF
      Utemp(2:4,l,k,j,iElem) = QDSMacroValues(2:4,l,k,j,iElem)/QDSMacroValues(1,l,k,j,iElem)
      Utemp(5:6,l,k,j,iElem) = QDSMacroValues(5:6,l,k,j,iElem)
    ELSE
      Utemp(:,l,k,j,iElem) = 0.0
    END IF
  END DO; END DO; END DO
END DO
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE QDS TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#endif
END IF
! generate nextfile info in previous output file
IF(PRESENT(PreviousTime))THEN
  IF(MPIRoot .AND. PreviousTime.LT.OutputTime) CALL GenerateNextFileInfo('QDS',OutputTime,PreviousTime)
END IF
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_QDS',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('QDS',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)!,FutureTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
CALL WriteAttributeToHDF5(File_ID,'VarNamesQDS',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        nQDSElems       => INT(nQDSElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        PP_N            => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName    = 'DG_Solution' , rank = 5                                           , &
                          nValGlobal     = (/N_variables , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                          nVal           = (/N_variables , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nQDSElems   /) , &
                          offset         = (/       0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem  /) , &
                          collective     = .TRUE.        , RealArray = Utemp)
END ASSOCIATE
#if USE_MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif
END SUBROUTINE WriteQDSToHDF5
#endif /*USE_QDS_DG*/


#ifdef PARTICLES
SUBROUTINE WriteIMDStateToHDF5()
!===================================================================================================================================
! Write the particles data aquired from an IMD *.chkpt file to disk and abort the program
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: IMDInputFile,IMDTimeScale,IMDLengthScale,IMDNumber
USE MOD_Mesh_Vars     ,ONLY: MeshFile
USE MOD_Restart_Vars  ,ONLY: DoRestart
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif /*USE_MPI*/
USE MOD_ReadInTools   ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255) :: tempStr
REAL               :: t,tFuture,IMDtimestep
INTEGER            :: iSTATUS,IMDanalyzeIter
!===================================================================================================================================
IF(.NOT.DoRestart)THEN
  IF(IMDTimeScale.GT.0.0)THEN
    SWRITE(UNIT_StdOut,'(A)')'   IMD: calc physical time in seconds for which the IMD *.chkpt file is defined.'
    ! calc physical time in seconds for which the IMD *.chkpt file is defined
    ! t = IMDanalyzeIter * IMDtimestep * IMDTimeScale * IMDNumber
    IMDtimestep=0.0
    CALL GetParameterFromFile(IMDInputFile,'timestep'   , TempStr ,DelimiterSymbolIN=' ',CommentSymbolIN='#')
    CALL str2real(TempStr,IMDtimestep,iSTATUS)
    IF(iSTATUS.NE.0)THEN
      CALL abort(&
      __STAMP__&
      ,'Could not find "timestep" in '//TRIM(IMDInputFile)//' for IMDtimestep!')
    END IF

    IMDanalyzeIter=0
    CALL GetParameterFromFile(IMDInputFile,'checkpt_int', TempStr ,DelimiterSymbolIN=' ',CommentSymbolIN='#')
    CALL str2int(TempStr,IMDanalyzeIter,iSTATUS)
    IF(iSTATUS.NE.0)THEN
      CALL abort(&
      __STAMP__&
      ,'Could not find "checkpt_int" in '//TRIM(IMDInputFile)//' for IMDanalyzeIter!')
    END IF
    CALL PrintOption('IMDtimestep'    , 'OUTPUT' , RealOpt=IMDtimestep)
    CALL PrintOption('IMDanalyzeIter' , 'OUTPUT' , IntOpt=IMDanalyzeIter)
    CALL PrintOption('IMDTimeScale'   , 'OUTPUT' , RealOpt=IMDTimeScale)
    CALL PrintOption('IMDLengthScale' , 'OUTPUT' , RealOpt=IMDLengthScale)
    CALL PrintOption('IMDNumber'      , 'OUTPUT' , IntOpt=IMDNumber)
    t = REAL(IMDanalyzeIter) * IMDtimestep * IMDTimeScale * REAL(IMDNumber)
    CALL PrintOption('t'              , 'OUTPUT' , RealOpt=t)
    SWRITE(UNIT_StdOut,'(A,ES25.14E3,A,F15.3,A)')     '   Calculated time t :',t,' (',t*1e12,' ps)'

    tFuture=t
    CALL WriteStateToHDF5(TRIM(MeshFile),t,tFuture)
    SWRITE(UNIT_StdOut,'(A)')'   Particles: StateFile (IMD MD data) created. Terminating successfully!'
#if USE_MPI
    CALL FinalizeMPI()
    CALL MPI_FINALIZE(iERROR)
    IF(iERROR.NE.0)THEN
      CALL abort(&
      __STAMP__&
      , ' MPI_FINALIZE(iERROR) returned non-zero integer value',iERROR)
    END IF
#endif /*USE_MPI*/
    STOP 0 ! terminate successfully
  ELSE
    CALL abort(&
    __STAMP__&
    , ' IMDLengthScale.LE.0.0 which is not allowed')
  END IF
END IF
END SUBROUTINE WriteIMDStateToHDF5
#endif /*PARTICLES*/


SUBROUTINE WriteBFieldToHDF5()
!===================================================================================================================================
! Subroutine to write the BField numerical solution to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems, nElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5, copy_userblock
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: BGField, NodeType, PsiMag
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal
REAL,ALLOCATABLE               :: outputArray(:,:,:,:,:)
#if USE_MPI
REAL                           :: StartT,EndT
#endif /*USE_MPI*/
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE BG-FIELD TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#endif /*USE_MPI*/

ALLOCATE(outputArray(1:4,0:PP_N,0:PP_N,0:PP_N,1:nElems))
outputArray(1,:,:,:,:) = BGField(1,:,:,:,:)
outputArray(2,:,:,:,:) = BGField(2,:,:,:,:)
outputArray(3,:,:,:,:) = BGField(3,:,:,:,:)
outputArray(4,:,:,:,:) = PsiMag(:,:,:,:)

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:4))
StrVarNames(1)='BG-MagneticFieldX'
StrVarNames(2)='BG-MagneticFieldY'
StrVarNames(3)='BG-MagneticFieldZ'
StrVarNames(4)='PsiMag'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_BField.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BField',File_ID) ! File_Type='BField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=N)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',4,StrArray=StrVarNames)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  PP_N         => INT(PP_N,IK)    ,&
  PP_nElems    => INT(PP_nElems,IK)            ,&
  offsetElem   => INT(offsetElem,IK)           ,&
  nGlobalElems => INT(nGlobalElems,IK)         )
CALL WriteArrayToHDF5(DataSetName='BField', rank=5,&
                      nValGlobal=(/4_IK,PP_N+1_IK,PP_N+1_IK,PP_N+1_IK,nGlobalElems/),&
                      nVal      =(/4_IK,PP_N+1_IK,PP_N+1_IK,PP_N+1_IK,PP_nElems/),&
                      offset    =(/0_IK,     0_IK,     0_IK,     0_IK,offsetElem/),&
                      collective=.false., RealArray=outputArray)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(outputArray)

#if USE_MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif /*USE_MPI*/
END SUBROUTINE WriteBFieldToHDF5


SUBROUTINE WriteBFieldAnalyticToHDF5()
!===================================================================================================================================
! Subroutine to write the BField analytical solution to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems, nElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5, copy_userblock
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: BGFieldAnalytic, NodeType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal
INTEGER,PARAMETER              :: nOutputDim = 3
REAL,ALLOCATABLE               :: outputArray(:,:,:,:,:)
#if USE_MPI
REAL                           :: StartT,EndT
#endif /*USE_MPI*/
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE BG-FIELD Analytic solution TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#endif /*USE_MPI*/

ALLOCATE(outputArray(1:nOutputDim,0:PP_N,0:PP_N,0:PP_N,1:nElems))
outputArray(1,:,:,:,:) = BGFieldAnalytic(1,:,:,:,:)
outputArray(2,:,:,:,:) = BGFieldAnalytic(2,:,:,:,:)
outputArray(3,:,:,:,:) = BGFieldAnalytic(3,:,:,:,:)

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:nOutputDim))
StrVarNames(1)='BG-MagneticFieldX'
StrVarNames(2)='BG-MagneticFieldY'
StrVarNames(3)='BG-MagneticFieldZ'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_BFieldAnalytic.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BField',File_ID) ! File_Type='BField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=N)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',nOutputDim,StrArray=StrVarNames)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  PP_N         => INT(PP_N,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) ,&
  nOutputDim   => INT(nOutputDim,IK)   )
CALL WriteArrayToHDF5(DataSetName='BField', rank=5,&
                      nValGlobal=(/nOutputDim , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                      nVal      =(/nOutputDim , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
                      collective=.false., RealArray=outputArray)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(outputArray)

#if USE_MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif /*USE_MPI*/
END SUBROUTINE WriteBFieldAnalyticToHDF5


END MODULE MOD_HDF5_Output_Tools
