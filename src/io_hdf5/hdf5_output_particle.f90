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

MODULE MOD_HDF5_Output_Particles
#if defined(PARTICLES)
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteNodeSourceExtToHDF5
  MODULE PROCEDURE WriteNodeSourceExtToHDF5
END INTERFACE

INTERFACE WriteParticleToHDF5
  MODULE PROCEDURE WriteParticleToHDF5
END INTERFACE

INTERFACE WriteBoundaryParticleToHDF5
  MODULE PROCEDURE WriteBoundaryParticleToHDF5
END INTERFACE

INTERFACE WriteLostParticlesToHDF5
  MODULE PROCEDURE WriteLostParticlesToHDF5
END INTERFACE

INTERFACE WriteAdaptiveInfoToHDF5
  MODULE PROCEDURE WriteAdaptiveInfoToHDF5
END INTERFACE

INTERFACE WriteAdaptiveWallTempToHDF5
  MODULE PROCEDURE WriteAdaptiveWallTempToHDF5
END INTERFACE

INTERFACE WriteVibProbInfoToHDF5
  MODULE PROCEDURE WriteVibProbInfoToHDF5
END INTERFACE

INTERFACE WriteClonesToHDF5
  MODULE PROCEDURE WriteClonesToHDF5
END INTERFACE

#if USE_HDG
PUBLIC :: AddBRElectronFluidToPartSource
#endif /*USE_HDG*/
PUBLIC :: WriteNodeSourceExtToHDF5
PUBLIC :: WriteParticleToHDF5
PUBLIC :: WriteBoundaryParticleToHDF5
PUBLIC :: WriteLostParticlesToHDF5
PUBLIC :: WriteAdaptiveInfoToHDF5
PUBLIC :: WriteAdaptiveWallTempToHDF5
PUBLIC :: WriteVibProbInfoToHDF5
PUBLIC :: WriteClonesToHDF5
PUBLIC :: WriteElectroMagneticPICFieldToHDF5
PUBLIC :: WriteEmissionVariablesToHDF5
PUBLIC :: FillParticleData
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteNodeSourceExtToHDF5(OutputTime)
!===================================================================================================================================
! Write NodeSourceExt (external charge density) field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars    ,ONLY: NodeSourceExtGlobal
USE MOD_Mesh_Vars          ,ONLY: MeshFile,nGlobalElems,offsetElem,Vdm_EQ_N
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt,NodeVolume,DoDeposition
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,nUniqueGlobalNodes
USE MOD_TimeDisc_Vars      ,ONLY: iter
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVISU
USE MOD_Interpolation      ,ONLY: GetVandermonde
#if USE_MPI
USE MOD_PICDepo            ,ONLY: ExchangeNodeSourceExtTmp
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: N_variables=1
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255)             :: FileName,DataSetName
INTEGER                        :: iElem,i,iMax
REAL                           :: NodeSourceExtEqui(1:N_variables,0:1,0:1,0:1),sNodeVol(1:8)
INTEGER                        :: NodeID(1:8)
!===================================================================================================================================
! create global Eps field for parallel output of Eps distribution
ALLOCATE(NodeSourceExtGlobal(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='NodeSourceExt'
NodeSourceExtGlobal=0.

! Allocate and determine Vandermonde mapping from equidistant (visu) to NodeType node set
IF(.NOT.ALLOCATED(Vdm_EQ_N))THEN
  ALLOCATE(Vdm_EQ_N(0:PP_N,0:1))
  CALL GetVandermonde(1, NodeTypeVISU, PP_N, NodeType, Vdm_EQ_N, modal=.FALSE.)
END IF ! .NOT.ALLOCATED(Vdm_EQ_N)

! Skip MPI communication in the first step as nothing has been deposited yet
IF(iter.NE.0)THEN

#if USE_MPI
! Communicate the NodeSourceExtTmp values of the last boundary interaction before the state is written to .h5
! Only call when deposition is active (otherwise this routine only writes the old array from the restart file to keep the data)
IF(DoDeposition) CALL ExchangeNodeSourceExtTmp()
#endif /*USE_MPI*/

end if ! iter.NE.0

! Loop over all elements and store charge density values in equidistantly distributed nodes of PP_N=1
sNodeVol(1:8) = 1. ! default
! This array is not allocated when DoDeposition=F, however, the previously calculated surface charge might still be required in
! the future, when DoDeposition is activated again. Therefore, read the old data and store in the new state file.
IF(.NOT.ALLOCATED(NodeSourceExt)) THEN
  ALLOCATE(NodeSourceExt(1:nUniqueGlobalNodes))
  NodeSourceExt = 0.
END IF
DO iElem=1,PP_nElems
  ! Copy values to equidistant distribution
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem+offsetElem)))
  IF(DoDeposition) sNodeVol(1:8) = 1./NodeVolume(NodeID(1:8)) ! only required when actual deposition is performed
  NodeSourceExtEqui(1,0,0,0) = NodeSourceExt(NodeID(1))*sNodeVol(1)
  NodeSourceExtEqui(1,1,0,0) = NodeSourceExt(NodeID(2))*sNodeVol(2)
  NodeSourceExtEqui(1,1,1,0) = NodeSourceExt(NodeID(3))*sNodeVol(3)
  NodeSourceExtEqui(1,0,1,0) = NodeSourceExt(NodeID(4))*sNodeVol(4)
  NodeSourceExtEqui(1,0,0,1) = NodeSourceExt(NodeID(5))*sNodeVol(5)
  NodeSourceExtEqui(1,1,0,1) = NodeSourceExt(NodeID(6))*sNodeVol(6)
  NodeSourceExtEqui(1,1,1,1) = NodeSourceExt(NodeID(7))*sNodeVol(7)
  NodeSourceExtEqui(1,0,1,1) = NodeSourceExt(NodeID(8))*sNodeVol(8)

  ! Map equidistant distribution to G/GL (current node type)
  CALL ChangeBasis3D(1, 1, PP_N, Vdm_EQ_N, NodeSourceExtEqui(:,:,:,:),NodeSourceExtGlobal(:,:,:,:,iElem))
END DO!iElem

! Write data twice to .h5 file
! 1. to _State_.h5 file (or restart)
! 2. to separate file (for visu)
#if USE_DEBUG
iMax=2 ! write to state and to a separate file (for debugging)
#else
iMax=1 ! write to state file
#endif /*USE_DEBUG*/
DO i = 1, iMax
  IF(i.EQ.1)THEN
    ! Write field to _State_.h5 file (or restart)
    FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))//'.h5'
    DataSetName='DG_SourceExt'
  ELSE
    ! Generate skeleton for the file with all relevant data on a single processor (MPIRoot)
    ! Write field to separate file for debugging purposes
    FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_NodeSourceExtGlobal',OutputTime))//'.h5'
    IF(MPIRoot) CALL GenerateFileSkeleton('NodeSourceExtGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
      CALL WriteAttributeToHDF5(File_ID,'VarNamesNodeSourceExtGlobal',N_variables,StrArray=StrVarNames)
      CALL CloseDataFile()
    END IF ! MPIRoot
    DataSetName='DG_Solution'
  END IF ! i.EQ.2

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
        DataSetName=TRIM(DataSetName) , rank=5 , &
        nValGlobal =(/N_variables     , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
        nVal       =(/N_variables     , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
        offset     =(/       0_IK     , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
        collective =.TRUE.            , RealArray=NodeSourceExtGlobal)
  END ASSOCIATE
END DO ! i = 1, 2

SDEALLOCATE(NodeSourceExtGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteNodeSourceExtToHDF5


#if USE_HDG
SUBROUTINE AddBRElectronFluidToPartSource()
!===================================================================================================================================
! Add BR electron fluid density to PartSource for output to state.h5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals            ,ONLY: abort
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_PreProc
USE MOD_HDG_Vars           ,ONLY: ElemToBRRegion,RegionElectronRef
USE MOD_DG_Vars            ,ONLY: U
USE MOD_PICDepo_Vars       ,ONLY: PartSource
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,RegionID
INTEGER :: i,j,k
REAL    :: source_e
!===================================================================================================================================
! Loop over all elements and all DOF and add the contribution of the BR electron density to PartSource
DO iElem=1,nElems
  ! BR electron fluid region
  RegionID=ElemToBRRegion(iElem)
  IF (RegionID.GT.0) THEN
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if ((USE_HDG) && (PP_nVar==1))
      source_e = U(1,i,j,k,iElem)-RegionElectronRef(2,RegionID)
#else
      CALL abort(__STAMP__,' CalculateBRElectronsPerCell only implemented for electrostatic HDG!')
#endif
      IF (source_e .LT. 0.) THEN
        source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
            * EXP( (source_e) / RegionElectronRef(3,RegionID) )
      ELSE
        source_e = RegionElectronRef(1,RegionID) &         !--- linearized boltzmann relation at positive exponent
            * (1. + ((source_e) / RegionElectronRef(3,RegionID)) )
      END IF
      PartSource(4,i,j,k,iElem) = PartSource(4,i,j,k,iElem) - source_e
    END DO; END DO; END DO
  END IF
END DO ! iElem=1,PP_nElems

END SUBROUTINE AddBRElectronFluidToPartSource
#endif /*USE_HDG*/


SUBROUTINE WriteParticleToHDF5(FileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems, offsetElem
USE MOD_Particle_Vars          ,ONLY: VarTimeStep
USE MOD_Particle_Vars          ,ONLY: PartInt,PartData,PartDataSize,locnPart,offsetnPart,PartIntSize,PartDataVarNames
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY: UseDSMC, DSMC
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Vars          ,ONLY: VibQuantData,ElecDistriData,AD_Data,MaxQuantNum,MaxElecQuant
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
LOGICAL                        :: reSwitch
!===================================================================================================================================

IF (MPIRoot) THEN
  ALLOCATE(StrVarNames(PartIntSize))
  StrVarNames(1) = 'FirstPartID'
  StrVarNames(2) = 'LastPartID'

  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesPartInt',PartIntSize,StrArray=StrVarNames)
  ! Write the beginning of the emission for restart from fluid solution
  CALL CloseDataFile()

  DEALLOCATE(StrVarNames)
END IF

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch      = .TRUE.
  gatheredWrite = .FALSE.
END IF

!-----------------------------------------------------
! 1. Basic particle properties
!-----------------------------------------------------
! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems          => INT(nGlobalElems,IK)          ,&
      nVar                  => INT(PartIntSize,IK)           ,&
      PP_nElems             => INT(PP_nElems,IK)             ,&
      offsetElem            => INT(offsetElem,IK)            ,&
      PartDataSize          => INT(PartDataSize,IK)          )

  CALL GatheredWriteArray(FileName                    , create = .FALSE.            , &
                          DataSetName     = 'PartInt' , rank   = 2                  , &
                          nValGlobal      = (/nVar    , nGlobalElems/)              , &
                          nVal            = (/nVar    , PP_nElems   /)              , &
                          offset          = (/0_IK    , offsetElem  /)              , &
                          collective      = .TRUE.    , IntegerArray = PartInt)

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',INT(PartDataSize,4),StrArray=PartDataVarNames)
    CALL CloseDataFile()
  END IF

  IF(nGlobalNbrOfParticles(3).EQ.0_IK)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName = 'PartData'     , rank = 2                 , &
                            nValGlobal  = (/PartDataSize , nGlobalNbrOfParticles(3) /) , &
                            nVal        = (/PartDataSize , locnPart   /)            , &
                            offset      = (/0_IK         , offsetnPart/)            , &
                            collective  = .FALSE.        , RealArray = PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !nGlobalNbrOfParticles(3).EQ.0_IK
#if USE_MPI
  CALL DistributedWriteArray(FileName                                                  , &
                             DataSetName  = 'PartData'      , rank= 2                  , &
                             nValGlobal   = (/PartDataSize  , nGlobalNbrOfParticles(3) /) , &
                             nVal         = (/PartDataSize  , locnPart   /)            , &
                             offset       = (/0_IK          , offsetnPart/)            , &
                             collective   = UseCollectiveIO , offSetDim= 2             , &
                             communicator = MPI_COMM_PICLAS    , RealArray= PartData)
  ! Output of the element-wise time step as a separate container in state file
  IF(VarTimeStep%UseDistribution) THEN
    CALL DistributedWriteArray(FileName                                      , &
                              DataSetName  = 'ElemTimeStep'  , rank = 2      , &
                              nValGlobal   = (/nGlobalElems  , 1_IK/)        , &
                              nVal         = (/PP_nElems     , 1_IK/)        , &
                              offset       = (/offsetElem    , 0_IK/)        , &
                              collective   = UseCollectiveIO , offSetDim = 1 , &
                              communicator = MPI_COMM_PICLAS    , RealArray = VarTimeStep%ElemFac)
  END IF
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'     , rank = 2                 , &
                        nValGlobal  = (/PartDataSize , nGlobalNbrOfParticles(3) /) , &
                        nVal        = (/PartDataSize , locnPart   /)            , &
                        offset      = (/0_IK         , offsetnPart/)            , &
                        collective  = .FALSE.        , RealArray = PartData)
  ! Output of the element-wise time step as a separate container in state file
  IF(VarTimeStep%UseDistribution) THEN
    CALL WriteArrayToHDF5(DataSetName = 'ElemTimeStep' , rank=2 , &
                          nValGlobal  = (/nGlobalElems , 1_IK/) , &
                          nVal        = (/PP_nElems    , 1_IK/) , &
                          offset      = (/offsetElem   , 0_IK/) , &
                          collective  = .FALSE.        , RealArray=VarTimeStep%ElemFac)
  END IF
  CALL CloseDataFile()
#endif /*USE_MPI*/
END ASSOCIATE

!-----------------------------------------------------
! 2. Polyatomic
!-----------------------------------------------------
IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (MaxQuantNum           => INT(MaxQuantNum,IK))
    IF(nGlobalNbrOfParticles(3).EQ.0_IK)THEN ! zero particles present: write empty dummy container to .h5 file, required for (auto-)restart
      IF(MPIRoot)THEN ! only root writes the container
        CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
        CALL WriteArrayToHDF5(DataSetName = 'VibQuantData' , rank = 2                 , &
                              nValGlobal  = (/MaxQuantNum  , nGlobalNbrOfParticles(3) /) , &
                              nVal        = (/MaxQuantNum  , locnPart   /)            , &
                              offset      = (/0_IK         , offsetnPart/)            , &
                              collective  = .FALSE.        , IntegerArray_i4 = VibQuantData)
        CALL CloseDataFile()
      END IF !MPIRoot
    END IF !nGlobalNbrOfParticles(3).EQ.0_IK
#if USE_MPI
    CALL DistributedWriteArray(FileName                                                 , &
                              DataSetName  = 'VibQuantData'  , rank = 2                 , &
                              nValGlobal   = (/MaxQuantNum   , nGlobalNbrOfParticles(3) /) , &
                              nVal         = (/MaxQuantNum   , locnPart    /)           , &
                              offset       = (/0_IK          , offsetnPart /)           , &
                              collective   = UseCollectiveIO , offSetDim = 2            , &
                              communicator = MPI_COMM_PICLAS    , IntegerArray_i4 = VibQuantData)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArrayToHDF5(DataSetName = 'VibQuantData' , rank = 2                 , &
                          nValGlobal  = (/ MaxQuantNum , nGlobalNbrOfParticles(3) /) , &
                          nVal        = (/ MaxQuantNum , locnPart     /)          , &
                          offset      = (/ 0_IK        , offsetnPart  /)          , &
                          collective  = .FALSE.        , IntegerArray_i4 = VibQuantData)
    CALL CloseDataFile()
#endif /*USE_MPI*/
  END ASSOCIATE

#if USE_LOADBALANCE
  IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
    DEALLOCATE(VibQuantData)
#if USE_LOADBALANCE
  END IF
#endif /*USE_LOADBALANCE*/
END IF ! useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)

!-----------------------------------------------------
! 3. Electronic Excitation
!-----------------------------------------------------
IF (useDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (MaxElecQuant          => INT(MaxElecQuant,IK))
    IF(nGlobalNbrOfParticles(3).EQ.0_IK)THEN ! zero particles present: write empty dummy container to .h5 file, required for (auto-)restart
      IF(MPIRoot)THEN ! only root writes the container
        CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
        CALL WriteArrayToHDF5(DataSetName = 'ElecDistriData' , rank = 2                 , &
                              nValGlobal  = (/MaxElecQuant   , nGlobalNbrOfParticles(3) /) , &
                              nVal        = (/MaxElecQuant   , locnPart   /)            , &
                              offset      = (/0_IK           , offsetnPart/)            , &
                              collective  = .FALSE.          , RealArray = ElecDistriData)
        CALL CloseDataFile()
      END IF !MPIRoot
    END IF !nGlobalNbrOfParticles(3).EQ.0_IK
#if USE_MPI
    CALL DistributedWriteArray(FileName                                                  , &
                              DataSetName  = 'ElecDistriData' , rank = 2                 , &
                              nValGlobal   = (/MaxElecQuant   , nGlobalNbrOfParticles(3) /) , &
                              nVal         = (/MaxElecQuant   , locnPart    /)           , &
                              offset       = (/0_IK           , offsetnPart /)           , &
                              collective   = UseCollectiveIO  , offSetDim = 2            , &
                              communicator = MPI_COMM_PICLAS     , RealArray = ElecDistriData)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArrayToHDF5(DataSetName = 'ElecDistriData' , rank = 2                    , &
                          nValGlobal  = (/ MaxElecQuant  , nGlobalNbrOfParticles(3)   /)  , &
                          nVal        = (/ MaxElecQuant  , locnPart     /)             , &
                          offset      = (/ 0_IK          , offsetnPart  /)             , &
                          collective  = .FALSE.          , RealArray = ElecDistriData)
    CALL CloseDataFile()
#endif /*USE_MPI*/
  END ASSOCIATE

#if USE_LOADBALANCE
  IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
    DEALLOCATE(ElecDistriData)
#if USE_LOADBALANCE
  END IF
#endif /*USE_LOADBALANCE*/
END IF

!-----------------------------------------------------
! 4. Ambipolar diffusion
!-----------------------------------------------------
IF (useDSMC.AND.DSMC%DoAmbipolarDiff) THEN
  IF(nGlobalNbrOfParticles(3).EQ.0_IK)THEN ! zero particles present: write empty dummy container to .h5 file, required for (auto-)restart
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName = 'ADVeloData' , rank = 2                 , &
                            nValGlobal  = (/3_IK       , nGlobalNbrOfParticles(3) /) , &
                            nVal        = (/3_IK       , locnPart   /)            , &
                            offset      = (/0_IK       , offsetnPart/)            , &
                            collective  = .FALSE.      , RealArray = AD_Data)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !nGlobalNbrOfParticles(3).EQ.0_IK
#if USE_MPI
  CALL DistributedWriteArray(FileName                                                 , &
                            DataSetName  = 'ADVeloData'    , rank = 2                 , &
                            nValGlobal   = (/3_IK          , nGlobalNbrOfParticles(3) /) , &
                            nVal         = (/3_IK          , locnPart    /)           , &
                            offset       = (/0_IK          , offsetnPart /)           , &
                            collective   = UseCollectiveIO , offSetDim = 2            , &
                            communicator = MPI_COMM_PICLAS    , RealArray = AD_Data)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'ADVeloData' , rank = 2                   , &
                        nValGlobal  = (/ 3_IK      , nGlobalNbrOfParticles(3)   /) , &
                        nVal        = (/ 3_IK      , locnPart     /)            , &
                        offset      = (/ 0_IK      , offsetnPart  /)            , &
                        collective  = .FALSE.      , RealArray = AD_Data)
  CALL CloseDataFile()
#endif /*USE_MPI*/

#if USE_LOADBALANCE
  IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
    DEALLOCATE(AD_Data)
#if USE_LOADBALANCE
  END IF
#endif /*USE_LOADBALANCE*/
END IF

! For LoadBalance, keep arrays allocated to restart from
#if defined(PARTICLES) && USE_LOADBALANCE
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
  DEALLOCATE(PartData)
  DEALLOCATE(PartInt)
#if defined(PARTICLES) && USE_LOADBALANCE
END IF
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

END SUBROUTINE WriteParticleToHDF5


SUBROUTINE WriteBoundaryParticleToHDF5(MeshFileName,OutputTime,PreviousTime)
!===================================================================================================================================
! Write data of impacting particles on specific boundary conditions of .h5 file (position, velocity, species ID, kinetic energy [eV],
! macro particle factor, time of impact, impact obliqueness angle)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ElementaryCharge
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_PreProc
USE MOD_Equation_Vars          ,ONLY: StrVarNames
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems, offsetElem
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,nVarPartStateBoundary
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart2
USE MOD_TimeDisc_Vars          ,ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: PreviousTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames2(:)
INTEGER                        :: nVar
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,globnPart(6)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER                        :: PartDataSizeLoc     !> Number of entries in PartData, not to be confused with the global variable
CHARACTER(LEN=255)             :: FileName,PreviousFileName
REAL                           :: PreviousTime_loc
INTEGER                        :: ALLOCSTAT,SpecID
!===================================================================================================================================
! Do not write to file on restart or fresh computation
IF(iter.EQ.0) RETURN
! set local variables for output and previous times
IF(PRESENT(PreviousTime))PreviousTime_loc = PreviousTime
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PartStateBoundary',OutputTime))//'.h5'

#if USE_HDG
#if PP_nVar==1
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateBoundary',4,StrVarNames,MeshFileName,OutputTime)
#elif PP_nVar==3
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateBoundary',3,StrVarNames,MeshFileName,OutputTime)
#else
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateBoundary',7,StrVarNames,MeshFileName,OutputTime)
#endif
#else
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateBoundary',PP_nVar,StrVarNames,MeshFileName,OutputTime)
#endif /*USE_HDG*/
! generate nextfile info in previous output file
IF(PRESENT(PreviousTime))THEN
  PreviousFileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PartStateBoundary',PreviousTime))//'.h5'
  IF(MPIRoot.AND.PreviousTime_loc.LT.OutputTime.AND.FILEEXISTS(PreviousFileName)) THEN
    CALL GenerateNextFileInfo('PartStateBoundary',OutputTime,PreviousTime_loc)
  END IF
END IF

! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

! 3xPos [m], 3xvelo [m/s], species [-]
PartDataSizeLoc = 7
! Kinetic energy [eV]
PartDataSizeLoc = PartDataSizeLoc + 1
! MPF [-]
PartDataSizeLoc = PartDataSizeLoc + 1
! time [s]
PartDataSizeLoc = PartDataSizeLoc + 1
! Impact obliqueness angle [degree]
PartDataSizeLoc = PartDataSizeLoc + 1
! iBC [-]
PartDataSizeLoc = PartDataSizeLoc + 1

! Set number of local particles
locnPart = INT(PartStateBoundaryVecLength,IK)

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('WriteBoundaryParticleToHDF5',offsetnPart,globnPart,locnPart,.FALSE.)

ALLOCATE(PartData(INT(PartDataSizeLoc,IK),offsetnPart+1_IK:offsetnPart+locnPart), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
    ,'Error in WriteBoundaryParticleToHDF5: Cannot allocate PartData array for writing boundary particle data to .h5!')

pcount=1
DO iPart=offsetnPart+1_IK,offsetnPart+locnPart
  ! Position and Velocity
  PartData(1,iPart)=PartStateBoundary(1,pcount)
  PartData(2,iPart)=PartStateBoundary(2,pcount)
  PartData(3,iPart)=PartStateBoundary(3,pcount)
  PartData(4,iPart)=PartStateBoundary(4,pcount)
  PartData(5,iPart)=PartStateBoundary(5,pcount)
  PartData(6,iPart)=PartStateBoundary(6,pcount)

  ! SpeciesID
  PartData(7,iPart)=PartStateBoundary(7,pcount)

  ! Kinetic energy [J->eV] (do not consider the MPF here! Call CalcEkinPart2 with MPF=1.0)
  ! Take ABS() from SpecID as is might be negative (for storing particles that are emitted from a surface)
  SpecID = INT(ABS(PartStateBoundary(7,pcount)))
  IF(SpecID.EQ.0) CALL abort(__STAMP__,'Error in WriteBoundaryParticleToHDF5: SpecID = PartStateBoundary(7,pcount) = 0')
  ! Check if photon is emitted during ray tracing
  IF(SpecID.EQ.999)THEN
    PartData(8,iPart)=0 ! insert photon energy here
  ELSE
    PartData(8,iPart)=CalcEkinPart2(PartStateBoundary(4:6,pcount),SpecID,1.0) / ElementaryCharge
  END IF ! SpecID.EQ.999

  ! MPF: Macro particle factor
  PartData(9,iPart)=PartStateBoundary(8,pcount)

  ! Simulation time [s]
  PartData(10,iPart)=PartStateBoundary(9,pcount)

  ! Impact obliqueness angle [degree]
  PartData(11,iPart)=PartStateBoundary(10,pcount)

  ! iBC [-]
  PartData(12,iPart)=PartStateBoundary(11,pcount)

  pcount = pcount +1
END DO ! iPart=offsetnPart+1_IK,offsetnPart+locnPart

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch=.TRUE.
  gatheredWrite=.FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems => INT(nGlobalElems,IK) ,&
      nVar         => INT(nVar,IK)         ,&
      PP_nElems    => INT(PP_nElems,IK)    ,&
      offsetElem   => INT(offsetElem,IK)   ,&
      PartDataSizeLoc => INT(PartDataSizeLoc,IK) )

  ALLOCATE(StrVarNames2(PartDataSizeLoc))
  StrVarNames2(1)  = 'ParticlePositionX'
  StrVarNames2(2)  = 'ParticlePositionY'
  StrVarNames2(3)  = 'ParticlePositionZ'
  StrVarNames2(4)  = 'VelocityX'
  StrVarNames2(5)  = 'VelocityY'
  StrVarNames2(6)  = 'VelocityZ'
  StrVarNames2(7)  = 'Species'
  StrVarNames2(8)  = 'KineticEnergy_eV'
  StrVarNames2(9)  = 'MacroParticleFactor'
  StrVarNames2(10) = 'Time'
  StrVarNames2(11) = 'ImpactObliquenessAngle'
  StrVarNames2(12) = 'iBC'

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',INT(PartDataSizeLoc,4),StrArray=StrVarNames2)
    CALL CloseDataFile()
  END IF

  IF(globnPart(3).EQ.0_IK)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName = 'PartData'        , rank = 2       , &
                            nValGlobal  = (/ PartDataSizeLoc, globnPart(3)/) , &
                            nVal        = (/ PartDataSizeLoc, locnPart    /) , &
                            offset      = (/ 0_IK           , offsetnPart /) , &
                            collective  = .FALSE.           , RealArray = PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !globnPart(3) .EQ.0_IK
#if USE_MPI
  CALL DistributedWriteArray(FileName                                                 , &
                             DataSetName  = 'PartData'        , rank = 2              , &
                             nValGlobal   = (/ PartDataSizeLoc, globnPart(3)/)        , &
                             nVal         = (/ PartDataSizeLoc, locnPart    /)        , &
                             offset       = (/ 0_IK           , offsetnPart /)        , &
                             collective   = UseCollectiveIO   , offSetDim = 2         , &
                             communicator = MPI_COMM_PICLAS      , RealArray = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'        , rank = 2              , &
                        nValGlobal  = (/ PartDataSizeLoc, globnPart(3)/)        , &
                        nVal        = (/ PartDataSizeLoc, locnPart    /)        , &
                        offset      = (/ 0_IK           , offsetnPart /)        , &
                        collective  = .FALSE.           , RealArray = PartData)
  CALL CloseDataFile()
#endif /*USE_MPI*/

END ASSOCIATE
! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

DEALLOCATE(StrVarNames2)
DEALLOCATE(PartData)

! Nullify and reset boundary parts container after write out
PartStateBoundaryVecLength = 0

! Re-allocate PartStateBoundary for a small number of particles and double the array size each time the
! maximum is reached
DEALLOCATE(PartStateBoundary)
ALLOCATE(PartStateBoundary(1:nVarPartStateBoundary,1:10))
PartStateBoundary=0.

END SUBROUTINE WriteBoundaryParticleToHDF5


SUBROUTINE WriteLostParticlesToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
! Write data of lost particles to .h5 file (position, velocity, species ID, MPF, time of loss, element ID and particle ID
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems, offsetElem
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLost,PartLostDataSize,PartStateLostVecLength,NbrOfLostParticles
USE MOD_Particle_Tracking_Vars ,ONLY: TotalNbrOfMissingParticlesSum
USE MOD_Equation_Vars          ,ONLY: StrVarNames
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames2(:)
INTEGER                        :: nVar
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,globnPart(6)
REAL,ALLOCATABLE               :: PartData(:,:)
CHARACTER(LEN=255)             :: FileName
INTEGER                        :: ALLOCSTAT
!===================================================================================================================================
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PartStateLost',OutputTime))//'.h5'

#if USE_HDG
#if PP_nVar==1
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateLost',4,StrVarNames,MeshFileName,OutputTime)
#elif PP_nVar==3
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateLost',3,StrVarNames,MeshFileName,OutputTime)
#else
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateLost',7,StrVarNames,MeshFileName,OutputTime)
#endif
#else
IF(MPIRoot) CALL GenerateFileSkeleton('PartStateLost',PP_nVar,StrVarNames,MeshFileName,OutputTime)
#endif /*USE_HDG*/

! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

! Set number of local particles
locnPart = INT(PartStateLostVecLength,IK)

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('WriteLostParticlesToHDF5',offsetnPart,globnPart,locnPart,.FALSE.)

ALLOCATE(PartData(INT(PartLostDataSize,IK),offsetnPart+1_IK:offsetnPart+locnPart), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'ERROR in hdf5_output.f90: Cannot allocate PartData array for writing lost particle data to .h5!')

pcount=1
DO iPart=offsetnPart+1_IK,offsetnPart+locnPart
  ! Position and Velocity
  PartData(1,iPart)=PartStateLost(1,pcount) ! LastPartPos-X
  PartData(2,iPart)=PartStateLost(2,pcount) ! LastPartPos-Y
  PartData(3,iPart)=PartStateLost(3,pcount) ! LastPartPos-Z
  PartData(4,iPart)=PartStateLost(4,pcount)
  PartData(5,iPart)=PartStateLost(5,pcount)
  PartData(6,iPart)=PartStateLost(6,pcount)

  ! SpeciesID
  PartData(7,iPart)=PartStateLost(7,pcount)

  ! MPF: Macro particle factor
  PartData(8,iPart)=PartStateLost(8,pcount)

  ! Simulation time [s]
  PartData(9,iPart)=PartStateLost(9,pcount)

  ! ElemID
  PartData(10,iPart)=PartStateLost(10,pcount)

  ! PartID
  PartData(11,iPart)=PartStateLost(11,pcount)

  ! PartPos (PartState(1:3))
  PartData(12,iPart)=PartStateLost(12,pcount)
  PartData(13,iPart)=PartStateLost(13,pcount)
  PartData(14,iPart)=PartStateLost(14,pcount)

  ! myrank
  PartData(15,iPart)=PartStateLost(15,pcount)

  ! MissingType
  PartData(16,iPart)=PartStateLost(16,pcount)

  pcount = pcount +1
END DO ! iPart=offsetnPart+1_IK,offsetnPart+locnPart

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch=.TRUE.
  gatheredWrite=.FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems     => INT(nGlobalElems,IK) ,&
      nVar             => INT(nVar,IK)         ,&
      PP_nElems        => INT(PP_nElems,IK)    ,&
      offsetElem       => INT(offsetElem,IK)   ,&
      PartLostDataSize => INT(PartLostDataSize,IK) )

  ALLOCATE(StrVarNames2(PartLostDataSize))
  StrVarNames2(1)  = 'LastPartPosX'
  StrVarNames2(2)  = 'LastPartPosY'
  StrVarNames2(3)  = 'LastPartPosZ'
  StrVarNames2(4)  = 'VelocityX'
  StrVarNames2(5)  = 'VelocityY'
  StrVarNames2(6)  = 'VelocityZ'
  StrVarNames2(7)  = 'Species'
  StrVarNames2(8)  = 'MacroParticleFactor'
  StrVarNames2(9)  = 'Time'
  StrVarNames2(10) = 'ElemID'
  StrVarNames2(11) = 'PartID'
  StrVarNames2(12)  = 'ParticlePositionX'
  StrVarNames2(13)  = 'ParticlePositionY'
  StrVarNames2(14)  = 'ParticlePositionZ'
  StrVarNames2(15)  = 'MyRank'
  StrVarNames2(16)  = 'MissingType'

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',INT(PartLostDataSize,4),StrArray=StrVarNames2)
    CALL CloseDataFile()
  END IF

  IF(globnPart(3).EQ.0_IK)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName = 'PartData'          , rank = 2       , &
                            nValGlobal  = (/ PartLostDataSize , globnPart(3)/) , &
                            nVal        = (/ PartLostDataSize , locnPart    /) , &
                            offset      = (/ 0_IK             , offsetnPart /) , &
                            collective  = .FALSE.             , RealArray = PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !globnPart(3).EQ.0_IK
#if USE_MPI
  CALL DistributedWriteArray(FileName                                                   , &
                             DataSetName  = 'PartData'          , rank = 2              , &
                             nValGlobal   = (/ PartLostDataSize , globnPart(3)/)        , &
                             nVal         = (/ PartLostDataSize , locnPart    /)        , &
                             offset       = (/ 0_IK             , offsetnPart /)        , &
                             collective   = UseCollectiveIO     , offSetDim = 2         , &
                             communicator = MPI_COMM_PICLAS        , RealArray = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'          , rank = 2              , &
                        nValGlobal  = (/ PartLostDataSize , globnPart(3)/)        , &
                        nVal        = (/ PartLostDataSize , locnPart    /)        , &
                        offset      = (/ 0_IK             , offsetnPart /)        , &
                        collective  = .FALSE.             , RealArray = PartData)
  CALL CloseDataFile()
#endif /*USE_MPI*/

END ASSOCIATE
! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

DEALLOCATE(StrVarNames2)
DEALLOCATE(PartData)

! Nullify and reset lost parts container after write out
PartStateLostVecLength  = 0
NbrOfLostParticles      = 0 ! only reset local counter but not the global counter (all procs)
TotalNbrOfMissingParticlesSum = 0 ! reset missing particle counter (only required during restart) after writing to .h5

! Re-allocate PartStateLost for a small number of particles and double the array size each time the maximum is reached
DEALLOCATE(PartStateLost)
ALLOCATE(PartStateLost(1:PartLostDataSize,1:10))
PartStateLost=0.

END SUBROUTINE WriteLostParticlesToHDF5


SUBROUTINE WriteAdaptiveInfoToHDF5(FileName)
!===================================================================================================================================
!> Subroutine that generates the adaptive boundary info and writes it out into State-File
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems, nElems
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_Particle_Sampling_Vars ,ONLY: AdaptBCMacroVal,AdaptBCSampleElemNum,AdaptBCMapSampleToElem,AdaptiveData,AdaptBCTruncAverage
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: UseAdaptiveType4
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255)             :: SpecID
INTEGER                        :: nVar, nVarTotal
INTEGER                        :: ElemID,iVar,iSpec,iSF,SampleElemID
!===================================================================================================================================

nVar = 7
iVar = 1
nVarTotal = nVar*nSpecies
ALLOCATE(StrVarNames(nVar*nSpecies))
DO iSpec=1,nSpecies
  WRITE(SpecID,'(I3.3)') iSpec
  StrVarNames(iVar)   = 'Spec'//TRIM(SpecID)//'-VeloX'
  StrVarNames(iVar+1) = 'Spec'//TRIM(SpecID)//'-VeloY'
  StrVarNames(iVar+2) = 'Spec'//TRIM(SpecID)//'-VeloZ'
  StrVarNames(iVar+3) = 'Spec'//TRIM(SpecID)//'-Density'
  StrVarNames(iVar+4) = 'Spec'//TRIM(SpecID)//'-PumpVeloPerArea'
  StrVarNames(iVar+5) = 'Spec'//TRIM(SpecID)//'-Pressure'
  StrVarNames(iVar+6) = 'Spec'//TRIM(SpecID)//'-PumpIntegralError'
  iVar = iVar + nVar
END DO

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesAdaptive',nVarTotal,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

iVar = 1
AdaptiveData = 0.
UseAdaptiveType4 = .FALSE.

DO iSpec = 1, nSpecies
  DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
    IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) UseAdaptiveType4 = .TRUE.
  END DO
  DO SampleElemID = 1,AdaptBCSampleElemNum
    ElemID = AdaptBCMapSampleToElem(SampleElemID)
    AdaptiveData(iVar:iVar-1+nVar,ElemID) = AdaptBCMacroVal(1:7,SampleElemID,iSpec)
  END DO
  iVar = iVar + nVar
END DO

WRITE(H5_Name,'(A)') 'AdaptiveInfo'
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      nElems          => INT(nElems,IK)          ,&
      nVarTotal       => INT(nVarTotal,IK)       ,&
      nSpecies        => INT(nSpecies,IK)        ,&
      offsetElem      => INT(offsetElem,IK)      )
  CALL WriteArrayToHDF5(DataSetName = H5_Name     , rank = 2                  , &
                        nValGlobal  = (/nVarTotal , nGlobalElems/)            , &
                        nVal        = (/nVarTotal , nElems      /)            , &
                        offset      = (/0_IK      , offsetElem  /)            , &
                        collective  = .FALSE.     , RealArray = AdaptiveData)
END ASSOCIATE
CALL CloseDataFile()
SDEALLOCATE(StrVarNames)

IF(AdaptBCTruncAverage) CALL WriteAdaptiveRunningAverageToHDF5(FileName)
IF(UseAdaptiveType4) CALL WriteAdaptBCPartNumOutToHDF5(FileName)

END SUBROUTINE WriteAdaptiveInfoToHDF5


SUBROUTINE WriteAdaptiveRunningAverageToHDF5(FileName)
!===================================================================================================================================
!> Output of the running average required for the sampling at the adaptive boundary conditions. Required to keep the average
!> during a load balance step.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Timedisc_Vars          ,ONLY: iter
USE MOD_Restart_Vars           ,ONLY: DoRestart
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_Particle_Sampling_Vars ,ONLY: AdaptBCAverage, AdaptBCSampleElemNum, AdaptBCMapSampleToElem, AdaptBCSampIter
USE MOD_Particle_Sampling_Vars ,ONLY: AdaptBCSampleElemNumGlobal, offSetElemAdaptBCSample, AdaptBCSampIterReadIn
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nVar,ElemID,SampleElemID
INTEGER, ALLOCATABLE           :: AdaptBCAverageIndex(:)
!===================================================================================================================================

IF(.NOT.DoRestart.AND.iter.EQ.0) RETURN

nVar = 8
ALLOCATE(AdaptBCAverageIndex(1:AdaptBCSampleElemNumGlobal))

DO SampleElemID = 1,AdaptBCSampleElemNum
  ElemID = AdaptBCMapSampleToElem(SampleElemID)
  AdaptBCAverageIndex(SampleElemID) = ElemID + offsetElem
END DO

! Store the position in the array for early restarts and the used weighting factor
IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'AdaptBCWeightingFactor',nSpecies,RealArray=Species(1:nSpecies)%MacroParticleFactor)
  IF(INT(iter,4)+AdaptBCSampIterReadIn.LT.AdaptBCSampIter) THEN
    CALL WriteAttributeToHDF5(File_ID,'AdaptBCSampIter',1,IntegerScalar=INT(iter,4)+AdaptBCSampIterReadIn)
  ELSE
    CALL WriteAttributeToHDF5(File_ID,'AdaptBCSampIter',1,IntegerScalar=AdaptBCSampIter)
  END IF
  CALL CloseDataFile()
END IF

CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      AdaptBCSampleElemNumGlobal    => INT(AdaptBCSampleElemNumGlobal,IK)    ,&
      AdaptBCSampleElemNum          => INT(AdaptBCSampleElemNum,IK)          ,&
      nVar                          => INT(nVar,IK)            ,&
      nSpecies                      => INT(nSpecies,IK)        ,&
      AdaptBCSampIter               => INT(AdaptBCSampIter,IK) ,&
      offSetElemAdaptBCSample       => INT(offSetElemAdaptBCSample,IK)      )
  CALL WriteArrayToHDF5(DataSetName = 'AdaptiveRunningAverage' , rank = 4                   , &
                        nValGlobal  = (/nVar  , AdaptBCSampIter  , AdaptBCSampleElemNumGlobal, nSpecies/) , &
                        nVal        = (/nVar  , AdaptBCSampIter  , AdaptBCSampleElemNum      , nSpecies/) , &
                        offset      = (/0_IK  , 0_IK             , offSetElemAdaptBCSample   , 0_IK/) , &
                        collective  = .FALSE. , RealArray = AdaptBCAverage)
  CALL WriteArrayToHDF5(DataSetName = 'AdaptiveRunningAverageIndex' , rank = 1                   , &
                        nValGlobal  = (/AdaptBCSampleElemNumGlobal/) , &
                        nVal        = (/AdaptBCSampleElemNum      /) , &
                        offset      = (/offSetElemAdaptBCSample  /) , &
                        collective  = .FALSE. , IntegerArray_i4 = AdaptBCAverageIndex)
END ASSOCIATE
CALL CloseDataFile()
DEALLOCATE(AdaptBCAverageIndex)

END SUBROUTINE WriteAdaptiveRunningAverageToHDF5


SUBROUTINE WriteAdaptBCPartNumOutToHDF5(FileName)
!===================================================================================================================================
!> Write the number of particles that left the domain in the previous time step to allow the continuation of the simulation
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Timedisc_Vars           ,ONLY: iter
USE MOD_Restart_Vars            ,ONLY: DoRestart
USE MOD_Particle_Vars           ,ONLY: nSpecies, Species
#if USE_MPI
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCPartNumOut
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nSurfacefluxBCs
INTEGER, ALLOCATABLE            :: AdaptBCPartNumOutTemp(:,:)
!===================================================================================================================================

IF(.NOT.DoRestart.AND.iter.EQ.0) RETURN

nSurfacefluxBCs = MAXVAL(Species(:)%nSurfacefluxBCs)

! Sum-up values across the processors for output in temporary variable
#if USE_MPI
IF(MPIRoot)THEN
  ALLOCATE(AdaptBCPartNumOutTemp(1:nSpecies,1:nSurfacefluxBCs))
  CALL MPI_REDUCE(AdaptBCPartNumOut,AdaptBCPartNumOutTemp,nSpecies*nSurfacefluxBCs,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
ELSE
  CALL MPI_REDUCE(AdaptBCPartNumOut,MPI_IN_PLACE         ,nSpecies*nSurfacefluxBCs,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
END IF
#endif

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nSpecies          => INT(nSpecies,IK)            ,&
            nSurfacefluxBCs   => INT(nSurfacefluxBCs,IK)        )
    CALL WriteArrayToHDF5(DataSetName = 'AdaptBCPartNumOut' , rank = 2, &
                          nValGlobal  = (/nSpecies,nSurfacefluxBCs/),   &
                          nVal        = (/nSpecies,nSurfacefluxBCs/),   &
                          offset      = (/0_IK,0_IK/), &
                          collective  = .FALSE. , IntegerArray_i4 = AdaptBCPartNumOutTemp)
  END ASSOCIATE
  CALL CloseDataFile()
END IF

END SUBROUTINE WriteAdaptBCPartNumOutToHDF5


SUBROUTINE WriteAdaptiveWallTempToHDF5(FileName)
!===================================================================================================================================
!> Output of the adaptive cell-local wall temperature and the corresponding global side index
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, nGlobalSurfSides, nComputeNodeSurfSides, offsetComputeNodeSurfSide
USE MOD_Particle_Boundary_Vars    ,ONLY: BoundaryWallTemp, SurfSide2GlobalSide
#if USE_MPI
USE MOD_MPI_Shared_Vars                ,ONLY: MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: H5_Name, H5_Name2
!===================================================================================================================================

#if USE_MPI
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nGlobalSurfSides      .EQ.0) RETURN
#endif

WRITE(H5_Name,'(A)') 'AdaptiveBoundaryWallTemp'
WRITE(H5_Name2,'(A)') 'AdaptiveBoundaryGlobalSideIndx'

#if USE_MPI
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nSurfSample          => INT(nSurfSample,IK)               , &
      nGlobalSides         => INT(nGlobalSurfSides,IK)           , &
      nLocalSides          => INT(nComputeNodeSurfSides,IK)     , &
      offsetSurfSide       => INT(offsetComputeNodeSurfSide,IK))
  CALL WriteArrayToHDF5(DataSetName = H5_Name , rank = 3                   , &
                        nValGlobal  = (/nSurfSample  , nSurfSample  , nGlobalSides/) , &
                        nVal        = (/nSurfSample  , nSurfSample  , nLocalSides      /) , &
                        offset      = (/0_IK  , 0_IK      , offsetSurfSide  /) , &
                        collective  = .FALSE.  , RealArray = BoundaryWallTemp(:,:,1:nComputeNodeSurfSides))
  CALL WriteArrayToHDF5(DataSetName = H5_Name2 , rank = 1                  , &
                        nValGlobal  = (/nGlobalSides/) , &
                        nVal        = (/nLocalSides/) , &
                        offset      = (/offsetSurfSide  /) , &
                        collective  = .FALSE.  , IntegerArray_i4 = SurfSide2GlobalSide(SURF_SIDEID,1:nComputeNodeSurfSides))
END ASSOCIATE
CALL CloseDataFile()

END SUBROUTINE WriteAdaptiveWallTempToHDF5


SUBROUTINE WriteVibProbInfoToHDF5(FileName)
!===================================================================================================================================
!> Subroutine that generates the adaptive boundary info and writes it out into State-File
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems, nElems
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: VarVibRelaxProb, CollisMode, DSMC
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
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255)             :: SpecID
INTEGER                        :: iSpec
!===================================================================================================================================
IF(CollisMode.GT.1) THEN
  IF(DSMC%VibRelaxProb.GE.2.0) THEN
    ALLOCATE(StrVarNames(nSpecies))
    DO iSpec=1,nSpecies
      WRITE(SpecID,'(I3.3)') iSpec
      StrVarNames(iSpec)   = 'Spec'//TRIM(SpecID)//'-VibProbRelax'
    END DO

    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttributeToHDF5(File_ID,'VarNamesVibProbInfo',nSpecies,StrArray=StrVarNames)
      CALL CloseDataFile()
    END IF

    WRITE(H5_Name,'(A)') 'VibProbInfo'
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nGlobalElems    => INT(nGlobalElems,IK)    ,&
          nElems          => INT(nElems,IK)          ,&
          nSpecies        => INT(nSpecies,IK)        ,&
          offsetElem      => INT(offsetElem,IK)      )
      CALL WriteArrayToHDF5(DataSetName = H5_Name        , rank = 2    , &
                            nValGlobal  = (/nGlobalElems , nSpecies /) , &
                            nVal        = (/nElems       , nSpecies /) , &
                            offset      = (/offsetElem   , 0_IK     /) , &
                            collective  = .FALSE.        , RealArray = VarVibRelaxProb%ProbVibAv)
    END ASSOCIATE
    CALL CloseDataFile()
    SDEALLOCATE(StrVarNames)
  ELSE ! DSMC%VibRelaxProb < 2.0
#if USE_MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif
    CALL WriteAttributeToHDF5(File_ID,'VibProbConstInfo',1,RealScalar=DSMC%VibRelaxProb)
    CALL CloseDataFile()
  END IF
ELSE ! CollisMode <= 1
  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VibProbConstInfo',1,RealScalar=0.)
    CALL CloseDataFile()
  END IF
END IF

END SUBROUTINE WriteVibProbInfoToHDF5


SUBROUTINE WriteClonesToHDF5(FileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars ,ONLY: ManualTimeStep
USE MOD_DSMC_Vars     ,ONLY: UseDSMC, CollisMode, DSMC, PolyatomMolDSMC, SpecDSMC
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting, ClonedParticles
USE MOD_PARTICLE_Vars ,ONLY: nSpecies, usevMPF, Species, PartDataSize
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
!INTEGER(HID_T)                 :: Dset_ID
!INTEGER                        :: nVal
INTEGER                        :: pcount, iDelay, iElem_glob, iPos
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,globnPart(6)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER, ALLOCATABLE           :: VibQuantData(:,:)
REAL, ALLOCATABLE              :: ElecDistriData(:,:), AD_Data(:,:)
INTEGER                        :: PartDataSizeLoc       !number of entries in each line of PartData
INTEGER                        :: MaxQuantNum, iPolyatMole, iSpec, tempDelay, MaxElecQuant
!-----------------------------------------------------------------------------------------------------------------------------------
! Additional output of clone delay and global element ID
PartDataSizeLoc = PartDataSize + 2

IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  MaxQuantNum = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

IF (useDSMC.AND.(DSMC%ElectronicModel.EQ.2)) THEN
  MaxElecQuant = 0
  DO iSpec = 1, nSpecies
    IF (.NOT.((SpecDSMC(iSpec)%InterID.EQ.4).OR.SpecDSMC(iSpec)%FullyIonized)) THEN
      IF (SpecDSMC(iSpec)%MaxElecQuant.GT.MaxElecQuant) MaxElecQuant = SpecDSMC(iSpec)%MaxElecQuant
    END IF
  END DO
END IF

locnPart =   0

SELECT CASE(RadialWeighting%CloneMode)
CASE(1)
  tempDelay = RadialWeighting%CloneInputDelay - 1
CASE(2)
  tempDelay = RadialWeighting%CloneInputDelay
CASE DEFAULT
  CALL abort(__STAMP__, 'RadialWeighting: CloneMode is not supported!')
END SELECT

DO pcount = 0,tempDelay
  locnPart = locnPart + RadialWeighting%ClonePartNum(pcount)
END DO

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('WriteClonesToHDF5',offsetnPart,globnPart,locnPart,.FALSE.)

ALLOCATE(PartData(PartDataSizeLoc,offsetnPart+1:offsetnPart+locnPart))

IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  ALLOCATE(VibQuantData(MaxQuantNum,offsetnPart+1:offsetnPart+locnPart))
  VibQuantData = 0
  !+1 is real number of necessary vib quants for the particle
END IF
IF (useDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
  ALLOCATE(ElecDistriData(MaxElecQuant,offsetnPart+1_IK:offsetnPart+locnPart))
  ElecDistriData = 0
  !+1 is real number of necessary vib quants for the particle
END IF
IF (useDSMC.AND.DSMC%DoAmbipolarDiff)  THEN
  ALLOCATE(AD_Data(3,offsetnPart+1_IK:offsetnPart+locnPart))
  AD_Data = 0
  !+1 is real number of necessary vib quants for the particle
END IF
iPart=offsetnPart
DO iDelay=0,tempDelay
  DO pcount = 1, RadialWeighting%ClonePartNum(iDelay)
    iElem_glob = ClonedParticles(pcount,iDelay)%Element
    iPart = iPart + 1
    PartData(1:6,iPart)=ClonedParticles(pcount,iDelay)%PartState(1:6)
    PartData(7,iPart)=REAL(ClonedParticles(pcount,iDelay)%Species)
    PartData(8,iPart)=REAL(iElem_glob)
    PartData(9,iPart)=REAL(iDelay)
    iPos = 9
    IF (useDSMC) THEN
      ! Internal energy modelling: vibrational + rotational
      IF(CollisMode.GT.1) THEN
        PartData(1+iPos:2+iPos,iPart) = ClonedParticles(pcount,iDelay)%PartStateIntEn(1:2)
        iPos = iPos + 2
        ! Electronic energy modelling
        IF(DSMC%ElectronicModel.GT.0) THEN
          PartData(1+iPos,iPart) = ClonedParticles(pcount,iDelay)%PartStateIntEn(3)
          iPos = iPos + 1
        END IF
      END IF
    END IF
    ! Variable particle weighting
    IF(usevMPF) THEN
      PartData(1+iPos,iPart) = ClonedParticles(pcount,iDelay)%WeightingFactor
      iPos = iPos + 1
    END IF
    IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
      IF (SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%SpecToPolyArray
        VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,iPart) = &
          ClonedParticles(pcount,iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
      ELSE
          VibQuantData(:,iPart) = 0
      END IF
    END IF
    IF (useDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
      IF (.NOT.((SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%InterID.EQ.4) &
          .OR.SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%FullyIonized)) THEN
          ElecDistriData(1:SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%MaxElecQuant,iPart) = &
            ClonedParticles(pcount,iDelay)%DistriFunc(1:SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%MaxElecQuant)
      END IF
    END IF
    IF (useDSMC.AND.DSMC%DoAmbipolarDiff)  THEN
      IF (Species(ClonedParticles(pcount,iDelay)%Species)%ChargeIC.GT.0.0) THEN
        AD_Data(1:3,iPart) = ClonedParticles(pcount,iDelay)%AmbiPolVelo(1:3)
      END IF
    END IF
  END DO
END DO

ALLOCATE(StrVarNames(PartDataSizeLoc))
StrVarNames(1)='ParticlePositionX'
StrVarNames(2)='ParticlePositionY'
StrVarNames(3)='ParticlePositionZ'
StrVarNames(4)='VelocityX'
StrVarNames(5)='VelocityY'
StrVarNames(6)='VelocityZ'
StrVarNames(7)='Species'
StrVarNames(8)='Element'
StrVarNames(9)='CloneDelay'

iPos = 9
! DSMC-specific variables
IF(useDSMC)THEN
  IF(CollisMode.GT.1) THEN
    ! Internal energy modelling: vibrational + rotational
    StrVarNames(1+iPos)='Vibrational'
    StrVarNames(2+iPos)='Rotational'
    iPos=iPos+2
    ! Electronic energy modelling
    IF(DSMC%ElectronicModel.GT.0) THEN
      StrVarNames(1+iPos)='Electronic'
      iPos=iPos+1
    END IF
  END IF
END IF
! Variable particle weighting
IF (usevMPF) THEN
  StrVarNames(1+iPos)='MPF'
  iPos=iPos+1
END IF

#if USE_MPI
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

ASSOCIATE (&
      offsetnPart     => INT(offsetnPart,IK)   ,&
      MaxQuantNum     => INT(MaxQuantNum,IK)   ,&
      MaxElecQuant    => INT(MaxElecQuant,IK)   ,&
      PartDataSizeLoc => INT(PartDataSizeLoc,IK)  )
CALL WriteArrayToHDF5(DataSetName='CloneData'     , rank=2          , &
                      nValGlobal=(/PartDataSizeLoc, globnPart(3)/)  , &
                      nVal=      (/PartDataSizeLoc, locnPart    /)  , &
                      offset=    (/0_IK           , offsetnPart /)  , &
                      collective=.FALSE.          , RealArray=PartData)
IF (useDSMC) THEN
  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    CALL WriteArrayToHDF5(DataSetName='CloneVibQuantData' , rank=2              , &
                          nValGlobal=(/MaxQuantNum        , globnPart(3)   /)   , &
                          nVal=      (/MaxQuantNum        , locnPart       /)   , &
                          offset=    (/0_IK               , offsetnPart    /)   , &
                          collective=.FALSE.              , IntegerArray_i4=VibQuantData)
    DEALLOCATE(VibQuantData)
  END IF
  IF (DSMC%ElectronicModel.EQ.2) THEN
    CALL WriteArrayToHDF5(DataSetName='CloneElecDistriData' , rank=2              , &
                          nValGlobal=(/MaxElecQuant       , globnPart(3)   /)   , &
                          nVal=      (/MaxElecQuant        , locnPart      /)   , &
                          offset=    (/0_IK               , offsetnPart    /)   , &
                          collective=.FALSE.              , RealArray=ElecDistriData)
    DEALLOCATE(ElecDistriData)
  END IF
  IF (DSMC%DoAmbipolarDiff) THEN
    CALL WriteArrayToHDF5(DataSetName='CloneADVeloData'   , rank=2              , &
                          nValGlobal=(/3_IK               , globnPart(3)   /)   , &
                          nVal=      (/3_IK               , locnPart       /)   , &
                          offset=    (/0_IK               , offsetnPart    /)   , &
                          collective=.FALSE.              , RealArray=AD_Data)
    DEALLOCATE(AD_Data)
  END IF

END IF
END ASSOCIATE

CALL CloseDataFile()

! Output of clone species variables as attribute to the dataset
IF(MPIRoot) THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesParticleClones',PartDataSizeLoc,StrArray=StrVarNames,DatasetName='CloneData')
  CALL WriteAttributeToHDF5(File_ID,'ManualTimeStep',1,RealScalar=ManualTimeStep,DatasetName='CloneData')
  CALL WriteAttributeToHDF5(File_ID,'WeightingFactor',1,RealScalar=Species(1)%MacroParticleFactor,DatasetName='CloneData')
  CALL WriteAttributeToHDF5(File_ID,'RadialWeightingFactor',1,RealScalar=RadialWeighting%PartScaleFactor,DatasetName='CloneData')
  CALL CloseDataFile()
END IF

DEALLOCATE(StrVarNames)
DEALLOCATE(PartData)

END SUBROUTINE WriteClonesToHDF5


!===================================================================================================================================
!> Store the magnetic filed acting on particles at each DOF for all elements to .h5
!===================================================================================================================================
SUBROUTINE WriteElectroMagneticPICFieldToHDF5()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems, nElems,MeshFile,Elem_xGP
USE MOD_Output_Vars            ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars     ,ONLY: NodeType
USE MOD_PICInterpolation_tools ,ONLY: GetExternalFieldAtParticle,GetEMField
USE MOD_Interpolation_Vars     ,ONLY: xGP
USE MOD_Restart_Vars           ,ONLY: RestartTime
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
INTEGER,PARAMETER              :: outputVars=6
REAL,ALLOCATABLE               :: outputArray(:,:,:,:,:)
REAL                           :: StartT,EndT
INTEGER                        :: iElem,i,j,k
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE PIC EM-FIELD TO HDF5 FILE...'
GETTIME(StartT)

ALLOCATE(outputArray(1:outputVars,0:PP_N,0:PP_N,0:PP_N,1:nElems))
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        ASSOCIATE( x => Elem_xGP(1,i,j,k,iElem), y => Elem_xGP(2,i,j,k,iElem), z => Elem_xGP(3,i,j,k,iElem))
          ! Superposition of the external and calculated electromagnetic field
          !   GetExternalFieldAtParticle : Get the 1 of 4 external fields (analytic, variable, etc.) at position x,y,z
          !                   GetEMField : Evaluate the electro-(magnetic) field using the reference position and return the field
          outputArray(1:6,i,j,k,iElem) = GetExternalFieldAtParticle((/x,y,z/)) + GetEMField(iElem,(/xGP(i),xGP(j),xGP(k)/))
        END ASSOCIATE
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem=1,PP_nElems


! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:outputVars))
StrVarNames(1)='ElectricFieldX'
StrVarNames(2)='ElectricFieldY'
StrVarNames(3)='ElectricFieldZ'
StrVarNames(4)='MagneticFieldX'
StrVarNames(5)='MagneticFieldY'
StrVarNames(6)='MagneticFieldZ'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_PIC-EMField.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BField',File_ID) ! File_Type='BField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=PP_N)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',outputVars,StrArray=StrVarNames)
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
  outputVars   => INT(outputVars,IK)   ,&
  N            => INT(PP_N,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName = 'DG_Solution' , rank = 5                                    , &
                      nValGlobal  = (/outputVars  , N+1_IK   , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal        = (/outputVars  , N+1_IK   , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset      = (/0_IK        , 0_IK     , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective  = .FALSE.       , RealArray = outputArray)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(outputArray)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

END SUBROUTINE WriteElectroMagneticPICFieldToHDF5


!===================================================================================================================================
!> Write particle emission variables from state.h5
!> E.g. arrays containing information that have to be restored after restart (not necessarily required for automatic load balance
!> restarts, but maybe required for some)
!> Synchronize the read-in variables across all procs within the emission communicator (for the specific Species and Init) if
!> required
!===================================================================================================================================
SUBROUTINE WriteEmissionVariablesToHDF5(FileName)
! MODULES
#if USE_MPI
USE mpi
#endif /*USE_MPI*/
!USE MOD_io_HDF5
USE MOD_Globals
!USE MOD_PreProc
USE MOD_Particle_Vars     ,ONLY: Species,nSpecies
USE MOD_Particle_Vars     ,ONLY: NeutralizationBalanceGlobal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec,iInit ! ,InitGroup
CHARACTER(LEN=50) :: InitName
INTEGER(KIND=IK)  :: NeutralizationBalanceTmp(1:1) ! This is a dummy array of size 1 !
!===================================================================================================================================
! Only root writes the data
IF(.NOT.MPIRoot) RETURN

! Loop over all species and inits
DO iSpec=1,nSpecies
  DO iInit = 1, Species(iSpec)%NumberOfInits
    SELECT CASE(Species(iSpec)%Init(iInit)%ParticleEmissionType)
     CASE(9) ! '2D_landmark_neutralization'
       ! Re-load the value because the emission communicator can change during load balance restarts: MPIRoot is always part of this
       ! specific communicator

       NeutralizationBalanceTmp(1) = NeutralizationBalanceGlobal

       WRITE(InitName,'(A,I0,A,I0)') 'Spec',iSpec,'Init',iInit
       CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
       ! Associate construct for integer KIND=8 possibility
       ASSOCIATE (&
             nGlobalEntries => INT(1,IK)  ,&
             nEntries       => INT(1,IK)  ,&
             offsetEntries  => INT(0,IK)  )
         CALL WriteArrayToHDF5(DataSetName = TRIM(InitName) , rank = 1 , &
                               nValGlobal  = (/nGlobalEntries/) , &
                               nVal        = (/nEntries      /) , &
                               offset      = (/offsetEntries /) , &
                               collective  = .FALSE. , IntegerArray = NeutralizationBalanceTmp)
       END ASSOCIATE
       CALL CloseDataFile()

     END SELECT
  END DO  ! iInit
END DO  ! iSpec=1,nSpecies

END SUBROUTINE WriteEmissionVariablesToHDF5


!===================================================================================================================================
!> Calculate the particle offset and global number of particles across all processors
!> In this routine the number are calculated using integer KIND=8, but are returned with integer KIND=ICC in order to test if using
!> integer KIND=8 is required for total number of particles, particle boundary state, lost particles or clones
!===================================================================================================================================
SUBROUTINE GetOffsetAndGlobalNumberOfParts(CallingRoutine,offsetnPart,globnPart,locnPart,GetMinMaxNbrOfParticles)
! MODULES
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: CallingRoutine
INTEGER(KIND=IK),INTENT(IN)  :: locnPart
LOGICAL,INTENT(IN)           :: GetMinMaxNbrOfParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=IK),INTENT(OUT) :: offsetnPart
INTEGER(KIND=IK),INTENT(OUT) :: globnPart(6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER(KIND=8)              :: locnPart8,locnPart8Recv,globnPart8 ! always integer KIND=8
INTEGER(KIND=IK)             :: SimNumSpecMin,SimNumSpecMax
#endif
!===================================================================================================================================
#if USE_MPI
locnPart8     = INT(locnPart,8)
locnPart8Recv = 0_IK
CALL MPI_EXSCAN(locnPart8,locnPart8Recv,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_PICLAS,iError)
offsetnPart = INT(locnPart8Recv,KIND=IK)
! Last proc calculates the global number and broadcasts it
IF(myrank.EQ.nProcessors-1) locnPart8=locnPart8Recv+locnPart8
CALL MPI_BCAST(locnPart8,1,MPI_INTEGER8,nProcessors-1,MPI_COMM_PICLAS,iError)
!global numbers
globnPart8=locnPart8
GlobalNbrOfParticlesUpdated = .TRUE.
LOGWRITE(*,*) TRIM(CallingRoutine)//'offsetnPart,locnPart,globnPart8',offsetnPart,locnPart,globnPart8

! Sanity check: Add up all particles with integer KIND=8 and compare
IF(MPIRoot)THEN
  ! Check if offsetnPart is kind=8 is the number of particles is larger than integer KIND=4
  IF(globnPart8.GT.INT(HUGE(offsetnPart),8)) THEN
    WRITE(UNIT_stdOut,'(A,I0)') '\n\n\nTotal number of particles  : ',globnPart8
    WRITE(UNIT_stdOut,'(A,I0)')       'Maximum number of particles: ',HUGE(offsetnPart)
    CALL abort(__STAMP__,TRIM(CallingRoutine)//' has encountered more than integer KIND=4 particles. Activate PICLAS_INTKIND8!')
  END IF
END IF ! MPIRoot

! Get min/max number of particles
SimNumSpecMin = 0
SimNumSpecMax = 0
IF(GetMinMaxNbrOfParticles)THEN
  IF (MPIRoot) THEN
    CALL MPI_REDUCE(locnPart , SimNumSpecMin , 1 , MPI_INTEGER_INT_KIND , MPI_MIN , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(locnPart , SimNumSpecMax , 1 , MPI_INTEGER_INT_KIND , MPI_MAX , 0 , MPI_COMM_PICLAS , IERROR)
  ELSE
    CALL MPI_REDUCE(locnPart , 0             , 1 , MPI_INTEGER_INT_KIND , MPI_MIN , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(locnPart , 0             , 1 , MPI_INTEGER_INT_KIND , MPI_MAX , 0 , MPI_COMM_PICLAS , IERROR)
  END IF
END IF ! GetMinMaxNbrOfParticles

! Cast to Kind=IK before returning the number
globnPart(1) = INT(SimNumSpecMin , KIND = IK)
globnPart(2) = INT(SimNumSpecMax , KIND = IK)
globnPart(3) = INT(globnPart8    , KIND = IK)

#else
offsetnPart=0_IK
globnPart(1:3)=INT(locnPart,KIND=IK)
#endif /*USE_MPI*/

! Get extrema over the complete simulation only during WriteParticleToHDF5
IF(GetMinMaxNbrOfParticles)THEN
  globnPart(4) = MIN(globnPart(1),globnPart(4))
  globnPart(5) = MAX(globnPart(2),globnPart(5))
  globnPart(6) = MAX(globnPart(3),globnPart(6))
END IF ! GetMinMaxNbrOfParticles

END SUBROUTINE GetOffsetAndGlobalNumberOfParts


SUBROUTINE FillParticleData()
!===================================================================================================================================
! Fills the particle data arrays required for loadbalance/output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Vars          ,ONLY: PartDataSize
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartSpecies,nSpecies
USE MOD_Particle_Vars          ,ONLY: PartInt,PartData
USE MOD_Particle_Vars          ,ONLY: locnPart,offsetnPart
USE MOD_Particle_Vars          ,ONLY: VibQuantData,ElecDistriData,AD_Data,MaxQuantNum,MaxElecQuant
USE MOD_Particle_Vars          ,ONLY: PartState, PartSpecies, PartMPF, usevMPF, nSpecies, Species
USE MOD_Particle_Vars          ,ONLY: UseRotRefFrame, PartVeloRotRef
USE MOD_DSMC_Vars              ,ONLY: UseDSMC, CollisMode,PartStateIntEn, DSMC, PolyatomMolDSMC, SpecDSMC, VibQuantsPar
USE MOD_DSMC_Vars              ,ONLY: ElectronicDistriPart, AmbipolElecVelo
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars          ,ONLY: velocityAtTime, velocityOutputAtTime
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: PartIntSize=2
INTEGER                        :: pcount, iPos
INTEGER(KIND=IK)               :: iPart
LOGICAL                        :: withDSMC=.FALSE.
INTEGER                        :: iPolyatMole, iSpec
INTEGER                        :: ALLOCSTAT
INTEGER                        :: iElem_glob, iElem_loc
!===================================================================================================================================

! Required default values for KIND=IK
MaxQuantNum=-1

withDSMC = useDSMC

IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  MaxQuantNum = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

IF (useDSMC.AND.(DSMC%ElectronicModel.EQ.2)) THEN
  MaxElecQuant = 0
  DO iSpec = 1, nSpecies
    IF (.NOT.((SpecDSMC(iSpec)%InterID.EQ.4).OR.SpecDSMC(iSpec)%FullyIonized)) THEN
      IF (SpecDSMC(iSpec)%MaxElecQuant.GT.MaxElecQuant) MaxElecQuant = SpecDSMC(iSpec)%MaxElecQuant
    END IF
  END DO
END IF

locnPart =   0_IK
DO pcount = 1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(pcount)) THEN
    locnPart = locnPart + 1_IK
  END IF
END DO

! Communicate the total number and offset of particles
CALL GetOffsetAndGlobalNumberOfParts('WriteParticleToHDF5',offsetnPart,nGlobalNbrOfParticles,locnPart,.TRUE.)

#if USE_LOADBALANCE
! Arrays might still be allocated from previous loadbalance step
SDEALLOCATE(PartInt)
SDEALLOCATE(PartData)
#endif /*USE_LOADBALANCE*/
ALLOCATE(PartInt(     PartIntSize,      offsetElem+1   :offsetElem+PP_nElems))
ALLOCATE(PartData(INT(PartDataSize,IK),offsetnPart+1_IK:offsetnPart+locnPart), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in hdf5_output.f90: Cannot allocate PartData array for writing particle data to .h5!')
!!! Kleiner Hack von JN (Teil 1/2):

IF (.NOT.(useDSMC.OR.usevMPF)) THEN
  ALLOCATE(PEM%pStart(1:PP_nElems)           , &
           PEM%pNumber(1:PP_nElems)          , &
           PEM%pNext(1:PDM%maxParticleNumber), &
           PEM%pEnd(1:PP_nElems) )
  useDSMC=.TRUE.
END IF
CALL UpdateNextFreePosition()
!!! Ende kleiner Hack von JN (Teil 1/2)
iPart=offsetnPart
DO iElem_loc=1,PP_nElems
  iElem_glob = iElem_loc + offsetElem
  PartInt(1,iElem_glob)=iPart
  IF (ALLOCATED(PEM%pNumber)) THEN
    nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
    PartInt(2,iElem_glob) = PartInt(1,iElem_glob) + INT(PEM%pNumber(iElem_loc),IK)
    pcount = PEM%pStart(iElem_loc)
    DO iPart=PartInt(1,iElem_glob)+1_IK,PartInt(2,iElem_glob)
      PartData(1:3,iPart)=PartState(1:3,pcount)
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
      IF (velocityOutputAtTime) THEN
        PartData(4:6,iPart)=velocityAtTime(1:3,pcount)
      ELSE
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
      PartData(4:6,iPart)=PartState(4:6,pcount)
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
      END IF
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
      PartData(7,iPart)=REAL(PartSpecies(pcount))
      ! Sanity check: output of particles with species ID zero is prohibited
      IF(PartData(7,iPart).LE.0) CALL abort(__STAMP__,&
          'Found particle for output to .h5 with species ID zero, which indicates a corrupted simulation.')
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(pcount.EQ.PARTOUT)THEN
          PartData(7,iPart)=-PartData(7,iPart)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      iPos = 7
      ! Rotational frame of reference
      IF(UseRotRefFrame) THEN
        PartData(1+iPos:3+iPos,iPart) = PartVeloRotRef(1:3,pcount)
        iPos = iPos + 3
      END IF
      ! DSMC-specific variable
      IF (withDSMC) THEN
        ! Internal energy modelling: vibrational + rotational
        IF(CollisMode.GT.1) THEN
          PartData(1+iPos:2+iPos,iPart) = PartStateIntEn(1:2,pcount)
          iPos = iPos + 2
          ! Electronic energy modelling
          IF(DSMC%ElectronicModel.GT.0) THEN
            PartData(1+iPos,iPart) = PartStateIntEn(3,pcount)
            iPos = iPos + 1
          END IF
        END IF
      END IF
      ! Variable particle weighting
      IF(usevMPF) THEN
        PartData(1+iPos,iPart) = PartMPF(pcount)
        iPos = iPos + 1
      END IF
      pcount = PEM%pNext(pcount)
    END DO
    iPart = PartInt(2,iElem_glob)
  ELSE
    CALL abort(__STAMP__, " Particle HDF5-Output method not supported! PEM%pNumber not associated")
  END IF
  PartInt(2,iElem_glob)=iPart
END DO

!-----------------------------------------------------
! 2. Polyatomic
!-----------------------------------------------------
IF(withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0))THEN
#if USE_LOADBALANCE
  SDEALLOCATE(VibQuantData)
#endif /*USE_LOADBALANCE*/
  ALLOCATE(VibQuantData(MaxQuantNum,offsetnPart+1_IK:offsetnPart+locnPart))
  VibQuantData = 0
  !+1 is real number of necessary vib quants for the particle
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(1,iElem_glob)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
      PartInt(2,iElem_glob) = PartInt(1,iElem_glob) + INT(PEM%pNumber(iElem_loc),IK)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(1,iElem_glob)+1_IK,PartInt(2,iElem_glob)
        IF (SpecDSMC(PartSpecies(pcount))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(PartSpecies(pcount))%SpecToPolyArray
          VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,iPart) = &
            VibQuantsPar(pcount)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
        ELSE
           VibQuantData(:,iPart) = 0
        END IF
        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(2,iElem_glob)
    ELSE
      CALL abort(&
      __STAMP__&
      , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(2,iElem_glob)=iPart
  END DO
END IF ! withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)

!-----------------------------------------------------
! 3. Electronic Excitation
!-----------------------------------------------------
IF(withDSMC.AND.(DSMC%ElectronicModel.EQ.2))THEN
#if USE_LOADBALANCE
  SDEALLOCATE(ElecDistriData)
#endif /*USE_LOADBALANCE*/
  ALLOCATE(ElecDistriData(MaxElecQuant,offsetnPart+1_IK:offsetnPart+locnPart))
  ElecDistriData = 0
  !+1 is real number of necessary vib quants for the particle
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(1,iElem_glob)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
      PartInt(2,iElem_glob) = PartInt(1,iElem_glob) + INT(PEM%pNumber(iElem_loc),IK)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(1,iElem_glob)+1_IK,PartInt(2,iElem_glob)
        IF (.NOT.((SpecDSMC(PartSpecies(pcount))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(pcount))%FullyIonized)) THEN
          ElecDistriData(1:SpecDSMC(PartSpecies(pcount))%MaxElecQuant,iPart) = &
            ElectronicDistriPart(pcount)%DistriFunc(1:SpecDSMC(PartSpecies(pcount))%MaxElecQuant)
        ELSE
           ElecDistriData(:,iPart) = 0
        END IF
        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(2,iElem_glob)
    ELSE
      CALL abort(&
      __STAMP__&
      , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(2,iElem_glob)=iPart
  END DO
END IF ! withDSMC.AND.(DSMC%ElectronicModel.EQ.2)

!-----------------------------------------------------
! 4. Ambipolar diffusion
!-----------------------------------------------------
IF(withDSMC.AND.DSMC%DoAmbipolarDiff)THEN
#if USE_LOADBALANCE
  SDEALLOCATE(AD_Data)
#endif /*USE_LOADBALANCE*/
  ALLOCATE(AD_Data(3,offsetnPart+1_IK:offsetnPart+locnPart))
  AD_Data = 0.0
  !+1 is real number of necessary vib quants for the particle
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(1,iElem_glob)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
      PartInt(2,iElem_glob) = PartInt(1,iElem_glob) + INT(PEM%pNumber(iElem_loc),IK)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(1,iElem_glob)+1_IK,PartInt(2,iElem_glob)
        IF (Species(PartSpecies(pcount))%ChargeIC.GT.0.0) THEN
          AD_Data(1:3,iPart) = AmbipolElecVelo(pcount)%ElecVelo(1:3)
        ELSE
          AD_Data(1:3,iPart) = 0
        END IF
        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(2,iElem_glob)
    ELSE
      CALL abort(__STAMP__, " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(2,iElem_glob)=iPart
  END DO

END IF ! withDSMC.AND.DSMC%DoAmbipolarDiff


!!! Kleiner Hack von JN (Teil 2/2):
useDSMC=withDSMC
IF (.NOT.(useDSMC.OR.usevMPF)) THEN
  DEALLOCATE(PEM%pStart , &
             PEM%pNumber, &
             PEM%pNext  , &
             PEM%pEnd   )
END IF
!!! Ende kleiner Hack von JN (Teil 2/2)

END SUBROUTINE FillParticleData

#endif /*defined(PARTICLES)*/
END MODULE MOD_HDF5_Output_Particles
