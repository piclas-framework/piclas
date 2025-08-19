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

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
PUBLIC :: WriteDielectricGlobalToHDF5
PUBLIC :: WriteErrorNormsToHDF5
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/
#if !(PP_TimeDiscMethod==700)
PUBLIC :: WriteBGFieldToHDF5,WriteBGFieldAnalyticToHDF5
#endif /*!(PP_TimeDiscMethod==700)*/
!===================================================================================================================================

CONTAINS

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
SUBROUTINE WriteDielectricGlobalToHDF5()
!===================================================================================================================================
! write DielectricGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars      ,ONLY: isDielectricElem,DielectricVol,ElemToDielectric
USE MOD_Mesh_Vars            ,ONLY: MeshFile,offsetElem,nElems
USE MOD_io_HDF5
USE MOD_ChangeBasis          ,ONLY: ChangeBasis3D
USE MOD_DG_vars              ,ONLY: N_DG_Mapping,nDofsMapping
USE MOD_HDF5_Output_ElemData ,ONLY: WriteAdditionalElemData
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
INTEGER,PARAMETER   :: nVarOut=2
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
INTEGER             :: iElem, Nloc
! p-adaption output
REAL,ALLOCATABLE    :: U_N_2D_local(:,:)
INTEGER             :: i,j,k,iDOF,nDOFOutput,offsetDOF
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! ---------------------------------------------------------
! Prepare U_N_2D_local array for output as DG_Solution
! ---------------------------------------------------------
! Get the number of output DOFs per processor as the difference between the first and last offset and adding the number of DOFs of the last element
nDOFOutput = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
! Get the offset based on the element-local polynomial degree
IF(offsetElem.GT.0) THEN
  offsetDOF = N_DG_Mapping(1,1+offsetElem)
ELSE
  offsetDOF = 0
END IF
! Allocate local 2D array: create global Eps field for parallel output of Eps distribution
ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))
ALLOCATE(StrVarNames(1:nVarOut))
StrVarNames(1)='DielectricEpsGlobal'
StrVarNames(2)='DielectricMuGlobal'

! Write into 2D array
iDOF = 0
DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  IF(isDielectricElem(iElem))THEN
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1,iDOF) = DielectricVol(ElemToDielectric(iElem))%DielectricEps(i,j,k)
      U_N_2D_local(2,iDOF) = DielectricVol(ElemToDielectric(iElem))%DielectricMu( i,j,k)
    END DO; END DO; END DO
  ELSE
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:nVarOut,iDOF) = 0.0
    END DO; END DO; END DO
  END IF ! isDielectricElem(iElem)
END DO!iElem
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE DielectricGlobal TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
CALL GenerateFileSkeleton('DielectricGlobal',nVarOut,StrVarNames,TRIM(MeshFile),OutputTime,FileNameOut=FileName)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesDielectricGlobal',nVarOut,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF ! MPIRoot

! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
CALL WriteAdditionalElemData(FileName,ElementOutNloc)

! ---------------------------------------------------------
! Output of DG_Solution
! ---------------------------------------------------------
! Associate construct for integer KIND=8 possibility
ASSOCIATE(nVarOut         => INT(nVarOut,IK)           ,&
          nDofsMapping    => INT(nDofsMapping,IK)      ,&
          nDOFOutput      => INT(nDOFOutput,IK)        ,&
          offsetDOF       => INT(offsetDOF,IK)         )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName = 'DG_Solution' , rank = 2                , &
                        nValGlobal  = (/nVarOut     , nDofsMapping/)          , &
                        nVal        = (/nVarOut     , nDOFOutput/)            , &
                        offset      = (/0_IK        , offsetDOF/)             , &
                        collective  = .TRUE.        , RealArray = U_N_2D_local)
END ASSOCIATE
SDEALLOCATE(U_N_2D_local)
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteDielectricGlobalToHDF5


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
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5, copy_userblock
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
CALL WriteArrayToHDF5(DataSetName='DG_Solution'    , rank=5 , &
                      nValGlobal=(/UexDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/UexDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective=.false., RealArray=Uex)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteErrorNormsToHDF5
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/


#if !(PP_TimeDiscMethod==700)
SUBROUTINE WriteBGFieldToHDF5(OutputTime)
!===================================================================================================================================
! Subroutine to write the BField numerical solution to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: copy_userblock
USE MOD_HDF5_Output        ,ONLY: WriteArrayToHDF5
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: NodeType, N_BG, BGDataSize, BGType
USE MOD_SuperB_Vars        ,ONLY: UseTimeDepCoil,nTimePoints,BGFieldFrequency,BGFieldCurrent
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_DG_vars            ,ONLY: N_DG_Mapping,nDofsMapping
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_HDF5_Output_ElemData ,ONLY: WriteAdditionalElemData
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL         :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nVarOut
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: iElem, Nloc
REAL                           :: StartT,EndT
! p-adaption output
REAL,ALLOCATABLE    :: U_N_2D_local(:,:),U_N_3D_local(:,:,:)
INTEGER             :: i,j,k,iDOF,nDOFOutput,offsetDOF
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

! ---------------------------------------------------------
! Prepare U_N_2D_local or U_N_3D_local (time-dependent) array for output as DG_Solution
! ---------------------------------------------------------
! Get the number of output DOFs per processor as the difference between the first and last offset and adding the number of DOFs of the last element
nDOFOutput = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
! Get the offset based on the element-local polynomial degree
IF(offsetElem.GT.0) THEN
  offsetDOF = N_DG_Mapping(1,1+offsetElem)
ELSE
  offsetDOF = 0
END IF
nVarOut = BGDataSize

! Allocate local 2D array: create global Eps field for parallel output of Eps distribution
IF(UseTimeDepCoil)THEN
  ALLOCATE(U_N_3D_local(1:nVarOut,1:nDOFOutput,1:nTimePoints))
ELSE
  ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))
END IF

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

  IF(MPIRoot)THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BGField',File_ID) ! File_Type='BGField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID   , 'N'                    , 1          , IntegerScalar=PP_N)
  CALL WriteAttributeToHDF5(File_ID   , 'MeshFile'             , 1          , StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID   , 'NodeType'             , 1          , StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID   , 'VarNames'             , BGDataSize , StrArray=StrVarNames)
  DEALLOCATE(StrVarNames)
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

! Write into 2D array
iDOF = 0
IF (UseTimeDepCoil) THEN
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_3D_local(1:BGDataSize,iDOF,1:nTimePoints) = N_BG(iElem)%BGFieldTDep(1:BGDataSize,i,j,k,1:nTimePoints)
    END DO; END DO; END DO
  END DO ! iElem = 1, nElems
ELSE
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:BGDataSize,iDOF) = N_BG(iElem)%BGField(1:BGDataSize,i,j,k)
    END DO; END DO; END DO
  END DO ! iElem = 1, nElems
END IF ! UseTimeDepCoil

! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
CALL WriteAdditionalElemData(FileName,ElementOutNloc)

! ---------------------------------------------------------
! Output of DG_Solution
! ---------------------------------------------------------
! Associate construct for integer KIND=8 possibility
ASSOCIATE(nVarOut         => INT(nVarOut,IK)           ,&
          nDofsMapping    => INT(nDofsMapping,IK)      ,&
          nDOFOutput      => INT(nDOFOutput,IK)        ,&
          nTimePoints     => INT(nTimePoints,IK)        ,&
          offsetDOF       => INT(offsetDOF,IK)         )
  IF (UseTimeDepCoil) THEN
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
                            DataSetName = 'DG_Solution' , rank = 3                       , &
                            nValGlobal  = (/nVarOut     , nDofsMapping  , nTimePoints/)  , &
                            nVal        = (/nVarOut     , nDOFOutput    , nTimePoints/)  , &
                            offset      = (/0_IK        , offsetDOF     , 0_IK /)        , &
                            collective  = .TRUE.        , RealArray = U_N_3D_local)
  ELSE
    CALL GatheredWriteArray(FileName, create = .FALSE.                            , &
                            DataSetName = 'DG_Solution' , rank = 2                , &
                            nValGlobal  = (/nVarOut     , nDofsMapping/)          , &
                            nVal        = (/nVarOut     , nDOFOutput/)            , &
                            offset      = (/0_IK        , offsetDOF/)             , &
                            collective  = .TRUE.        , RealArray = U_N_2D_local)
  END IF ! UseTimeDepCoil
  SDEALLOCATE(U_N_2D_local)
  SDEALLOCATE(U_N_3D_local)
END ASSOCIATE

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
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5,copy_userblock
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: N_BG, NodeType, BGDataSize, NMax, PREF_VDM
USE MOD_Restart_Vars       ,ONLY: RestartTime
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
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
INTEGER                        :: nVal, iElem, Nloc
REAL                           :: outputArray(1:BGDataSize,0:NMax,0:NMax,0:NMax,1:nElems)
REAL                           :: StartT,EndT
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE BG-FIELD Analytic solution TO HDF5 FILE...'
GETTIME(StartT)

DO iElem = 1, INT(PP_nElems)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  IF(Nloc.EQ.Nmax)THEN
    outputArray(:,:,:,:,iElem) = N_BG(iElem)%BGFieldAnalytic(:,:,:,:)
  ELSE
    CALL ChangeBasis3D(BGDataSize,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
        N_BG(iElem)%BGFieldAnalytic(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc ) , &
                        outputArray(1:BGDataSize , 0:NMax , 0:NMax , 0:NMax   , iElem))
  END IF ! Nloc.Eq.Nmax
END DO ! iElem = 1, nElems

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
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=Nmax)
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
  N            => INT(NMax,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName='DG_Solution', rank=5 , &
                      nValGlobal=(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK        , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective=.false., RealArray=outputArray)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBGFieldAnalyticToHDF5
#endif /*!(PP_TimeDiscMethod==700)*/


END MODULE MOD_HDF5_Output_Fields