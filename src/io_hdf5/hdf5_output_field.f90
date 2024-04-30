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

MODULE MOD_HDF5_Output_Fields
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

#if USE_HDG
#if defined(PARTICLES)
PUBLIC :: WriteBRAverageElemToHDF5
#endif /*defined(PARTICLES)*/
#else
PUBLIC :: WritePMLzetaGlobalToHDF5
#endif /*USE_HDG*/

PUBLIC :: WriteDielectricGlobalToHDF5,WriteBGFieldToHDF5,WriteBGFieldAnalyticToHDF5
#if (PP_nVar==8)
PUBLIC :: WritePMLDataToHDF5
#endif
PUBLIC :: WriteErrorNormsToHDF5
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
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
REAL                :: StartT,EndT
REAL                :: OutputTime
INTEGER             :: iElem
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
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
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE DielectricGlobal TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DielectricGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('DielectricGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesDielectricGlobal',N_variables,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF ! MPIRoot

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5 , &
                          nValGlobal =(/N_variables , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=DielectricGlobal)
END ASSOCIATE
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SDEALLOCATE(DielectricGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteDielectricGlobalToHDF5


#if USE_HDG
#if defined(PARTICLES)
SUBROUTINE WriteBRAverageElemToHDF5(isBRAverageElem)
!===================================================================================================================================
! write BRAverageElem field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars     ,ONLY: ProjectName
USE MOD_io_HDF5
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: isBRAverageElem(1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER  :: N_variables=1
REAL               :: BRAverageElem(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
CHARACTER(LEN=255) :: StrVarNames(1:N_variables)
CHARACTER(LEN=255) :: FileName
REAL               :: StartT,EndT
REAL               :: OutputTime
INTEGER            :: iElem
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! create global zeta field for parallel output of zeta distribution
StrVarNames(1)='BRAverageElem'
BRAverageElem=0.
DO iElem=1,PP_nElems
  IF(isBRAverageElem(iElem))THEN
    BRAverageElem(:,:,:,:,iElem) = 1.0
  END IF
END DO!iElem
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE BRAverageElem TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BRAverageElem',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('BRAverageElem',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
CALL WriteAttributeToHDF5(File_ID,'VarNamesBRAverageElem',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5 , &
                          nValGlobal =(/N_variables , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=BRAverageElem)
END ASSOCIATE
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBRAverageElemToHDF5
#endif /*defined(PARTICLES)*/


#else
SUBROUTINE WritePMLzetaGlobalToHDF5()
!===================================================================================================================================
! write PMLzetaGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PML_Vars         ,ONLY: PMLzetaGlobal,PMLzeta0,PMLzeta,isPMLElem,ElemToPML
USE MOD_Mesh_Vars        ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars     ,ONLY: ProjectName
USE MOD_io_HDF5
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=3
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
INTEGER             :: iElem
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
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

SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE PMLZetaGlobal TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PMLZetaGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('PMLZetaGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
CALL WriteAttributeToHDF5(File_ID,'VarNamesPMLzetaGlobal',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5 , &
                          nValGlobal =(/N_variables , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=PMLzetaGlobal)
END ASSOCIATE
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SDEALLOCATE(PMLzetaGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WritePMLzetaGlobalToHDF5
#endif /*USE_HDG*/


SUBROUTINE WriteBGFieldToHDF5(OutputTime)
!===================================================================================================================================
! Subroutine to write the BField numerical solution to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems, nElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: copy_userblock
USE MOD_HDF5_Output        ,ONLY: WriteArrayToHDF5
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: BGField, NodeType, NBG, BGDataSize, BGType
USE MOD_SuperB_Vars        ,ONLY: UseTimeDepCoil,nTimePoints,BGFieldTDep,BGFieldFrequency,BGFieldCurrent
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL         :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal
REAL                           :: StartT,EndT
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
IF(PRESENT(OutputTime))THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BGField',OutputTime))//'.h5'
ELSE
  FileName=TRIM(ProjectName)//'_BGField.h5'
END IF ! PRESENT(OutputTime)

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE BG-FIELD ['//TRIM(FileName)//'] TO HDF5 FILE...'
GETTIME(StartT)

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:BGDataSize))
IF(BGType.EQ.1) THEN
  StrVarNames(1)='BG-ElectricFieldX'
  StrVarNames(2)='BG-ElectricFieldY'
  StrVarNames(3)='BG-ElectricFieldZ'
ELSE IF(BGType.EQ.2) THEN
  StrVarNames(1)='BG-MagneticFieldX'
  StrVarNames(2)='BG-MagneticFieldY'
  StrVarNames(3)='BG-MagneticFieldZ'
ELSE IF(BGType.EQ.3) THEN
  StrVarNames(1)='BG-ElectricFieldX'
  StrVarNames(2)='BG-ElectricFieldY'
  StrVarNames(3)='BG-ElectricFieldZ'
  StrVarNames(4)='BG-MagneticFieldX'
  StrVarNames(5)='BG-MagneticFieldY'
  StrVarNames(6)='BG-MagneticFieldZ'
END IF

IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BGField',File_ID) ! File_Type='BGField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID   , 'N'                    , 1          , IntegerScalar=NBG)
  CALL WriteAttributeToHDF5(File_ID   , 'MeshFile'             , 1          , StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID   , 'NodeType'             , 1          , StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID   , 'VarNames'             , BGDataSize , StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID   , 'Time'                 , 1          , RealScalar=0.)
  IF(UseTimeDepCoil)THEN
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldFrequency'     , 1          , RealScalar=BGFieldFrequency)
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldCurrent'       , 1          , RealScalar=BGFieldCurrent)
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldTimeDependent' , 1          , LogicalScalar=.TRUE.)
  ELSE
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldTimeDependent' , 1          , LogicalScalar=.FALSE.)
  END IF
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  BGDataSize   => INT(BGDataSize,IK)    ,&
  N            => INT(PP_N,IK)          ,&
  PP_nElems    => INT(PP_nElems,IK)     ,&
  offsetElem   => INT(offsetElem,IK)    ,&
  nGlobalElems => INT(nGlobalElems,IK)  ,&
  nTimePoints  => INT(nTimePoints,IK)    )

  IF(UseTimeDepCoil)THEN
    CALL WriteArrayToHDF5(DataSetName='DG_Solution'   , rank=6 , &
                          nValGlobal=(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems, nTimePoints/) , &
                          nVal      =(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems   , nTimePoints/) , &
                          offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem  , 0_IK       /) , &
                          collective=.false., RealArray=BGFieldTDep(1:BGDataSize,0:N,0:N,0:N,1:nElems,1:nTimePoints))
  ELSE
    CALL WriteArrayToHDF5(DataSetName='DG_Solution'   , rank=5 , &
                          nValGlobal=(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal      =(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                          offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                          collective=.false., RealArray=BGField(1:BGDataSize,0:N,0:N,0:N,1:nElems))
  END IF ! UseTimeDepCoil

END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBGFieldToHDF5


SUBROUTINE WriteBGFieldAnalyticToHDF5()
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
USE MOD_Interpolation_Vars ,ONLY: BGFieldAnalytic, NodeType, BGDataSize
USE MOD_Restart_Vars       ,ONLY: RestartTime
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
REAL                           :: StartT,EndT
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE BG-FIELD Analytic solution TO HDF5 FILE...'
GETTIME(StartT)

ALLOCATE(outputArray(1:BGDataSize,0:PP_N,0:PP_N,0:PP_N,1:nElems))
outputArray(1,:,:,:,:) = BGFieldAnalytic(1,:,:,:,:)
outputArray(2,:,:,:,:) = BGFieldAnalytic(2,:,:,:,:)
outputArray(3,:,:,:,:) = BGFieldAnalytic(3,:,:,:,:)

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:BGDataSize))
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
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=PP_N)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',BGDataSize,StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'Time'    ,1,RealScalar=RestartTime)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  BGDataSize   => INT(BGDataSize,IK)   ,&
  N            => INT(PP_N,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName='DG_Solution'    , rank=5 , &
                      nValGlobal=(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective=.false., RealArray=outputArray)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(outputArray)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBGFieldAnalyticToHDF5


#if (PP_nVar==8)
SUBROUTINE WritePMLDataToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:offsetElem,nGlobalElems
USE MOD_PML_Vars      ,ONLY:DoPML,PMLToElem,U2,nPMLElems,PMLnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL,ALLOCATABLE               :: UPML(:,:,:,:,:)
INTEGER                        :: iPML
!===================================================================================================================================

IF(DoPML)THEN
  ALLOCATE(StrVarNames(PMLnVar))
  StrVarNames( 1)='PML-ElectricFieldX-P1'
  StrVarNames( 2)='PML-ElectricFieldX-P2'
  StrVarNames( 3)='PML-ElectricFieldX-P3'
  StrVarNames( 4)='PML-ElectricFieldY-P4'
  StrVarNames( 5)='PML-ElectricFieldY-P5'
  StrVarNames( 6)='PML-ElectricFieldY-P6'
  StrVarNames( 7)='PML-ElectricFieldZ-P7'
  StrVarNames( 8)='PML-ElectricFieldZ-P8'
  StrVarNames( 9)='PML-ElectricFieldZ-P9'
  StrVarNames(10)='PML-MagneticFieldX-P10'
  StrVarNames(11)='PML-MagneticFieldX-P11'
  StrVarNames(12)='PML-MagneticFieldX-P12'
  StrVarNames(13)='PML-MagneticFieldY-P13'
  StrVarNames(14)='PML-MagneticFieldY-P14'
  StrVarNames(15)='PML-MagneticFieldY-P15'
  StrVarNames(16)='PML-MagneticFieldZ-P16'
  StrVarNames(17)='PML-MagneticFieldZ-P17'
  StrVarNames(18)='PML-MagneticFieldZ-P18'
  StrVarNames(19)='PML-PhiB-P19'
  StrVarNames(20)='PML-PhiB-P20'
  StrVarNames(21)='PML-PhiB-P21'
  StrVarNames(22)='PML-PsiE-P22'
  StrVarNames(23)='PML-PsiE-P23'
  StrVarNames(24)='PML-PsiE-P24'

  ALLOCATE(UPML(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  UPML=0.0
  DO iPML=1,nPMLElems
    UPML(:,:,:,:,PMLToElem(iPML)) = U2(:,:,:,:,iPML)
  END DO ! iPML

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesPML',PMLnVar,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      N               => INT(PP_N,IK)            ,&
      PMLnVar         => INT(PMLnVar,IK)         ,&
      PP_nElems       => INT(PP_nElems,IK)       ,&
      offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName = 'PML_Solution', rank = 5,&
                          nValGlobal  = (/PMLnVar , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal        = (/PMLnVar , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                          offset      = (/0_IK    , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                          collective  = .TRUE.,RealArray = UPML)
END ASSOCIATE

!  CALL WriteArrayToHDF5(DataSetName='PML_Solution', rank=5,&
!                      nValGlobal=(/5,N+1,N+1,N+1,nGlobalElems/),&
!                      nVal=      (/5,N+1,N+1,N+1,PP_nElems/),&
!                      offset=    (/0,      0,     0,     0,     offsetElem/),&
!                      collective=.TRUE., existing=.FALSE., RealArray=UPML)
!
!  CALL CloseDataFile()
  DEALLOCATE(UPML)
  DEALLOCATE(StrVarNames)
END IF ! DoPML


END SUBROUTINE WritePMLDataToHDF5
#endif


SUBROUTINE WriteErrorNormsToHDF5(OutputTime)
!===================================================================================================================================
! Output the exact solution, the L2 error and LInf error to (in NodeTypeGL = 'GAUSS-LOBATTO') to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5,copy_userblock
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: Uex,NAnalyze,NodeTypeGL
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!USE MOD_Equation_Vars      ,ONLY: StrVarNames
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal
INTEGER,PARAMETER              :: UexDataSize=1
REAL                           :: StartT,EndT
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE Analytic solution and L2/LInf norms TO HDF5 FILE...'
GETTIME(StartT)

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:UexDataSize))
StrVarNames(1)='Phi'
!StrVarNames(2)='BG-MagneticFieldY'
!StrVarNames(3)='BG-MagneticFieldZ'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
!FileName=TRIM(ProjectName)//'_ErrorNorms.h5'
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_ErrorNorms',OutputTime))//'.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('DG_Solution',File_ID)
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=PP_N)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeGL/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',UexDataSize,StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'Time'    ,1,RealScalar=OutputTime)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  UexDataSize  => INT(UexDataSize,IK)   ,&
  N            => INT(NAnalyze,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName='DG_Solution', rank=5 , &
                      nValGlobal=(/UexDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/UexDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK        , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective=.false., RealArray=Uex)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteErrorNormsToHDF5


END MODULE MOD_HDF5_Output_Fields
