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

MODULE MOD_HDF5_Output_Particles
#if defined(PARTICLES)
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
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt,NodeVolume,NodeSourceExtTmp
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: nUniqueGlobalNodes
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExtTmpLoc, NodeMapping, NodeSourceExtTmp_Shared_Win,NodeSourceExt_Shared_Win
USE MOD_MPI_Shared         ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_LEADERS_SHARED, MPI_COMM_SHARED, myComputeNodeRank, myLeaderGroupRank
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors, nLeaderGroupProcs
#endif  /*USE_MPI*/
USE MOD_TimeDisc_Vars      ,ONLY: iter
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
REAL                           :: NodeSourceExtEqui(1:N_variables,0:1,0:1,0:1)
INTEGER                        :: NodeID(1:8)
#if USE_MPI
INTEGER                        :: iProc
INTEGER                        :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER                        :: MessageSize
#endif /*USE_MPI*/
INTEGER                        :: firstNode, lastNode, iNode
!===================================================================================================================================
! create global Eps field for parallel output of Eps distribution
ALLOCATE(NodeSourceExtGlobal(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='NodeSourceExt'
NodeSourceExtGlobal=0.

! Skip MPI communication in the first step as nothing has been deposited yet
IF(iter.NE.0)THEN

  ! Communicate the NodeSourceExtTmp values of the last boundary interaction before the state is written to .h5
#if USE_MPI
  MessageSize = nUniqueGlobalNodes
  CALL MPI_REDUCE(NodeSourceExtTmpLoc(:) ,NodeSourceExtTmp(:), &
      MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  ! Reset local surface charge
  NodeSourceExtTmpLoc = 0.
  CALL BARRIER_AND_SYNC(NodeSourceExtTmp_Shared_Win,MPI_COMM_SHARED)
  IF ((myComputeNodeRank.EQ.0).AND.(nLeaderGroupProcs.GT.1)) THEN
    DO iProc = 0, nLeaderGroupProcs - 1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        CALL MPI_IRECV( NodeMapping(iProc)%RecvNodeSourceExt(:)        &
            , NodeMapping(iProc)%nRecvUniqueNodes           &
            , MPI_DOUBLE_PRECISION                                        &
            , iProc                                                       &
            , 666                                                         &
            , MPI_COMM_LEADERS_SHARED                                       &
            , RecvRequest(iProc)                                          &
            , IERROR)
      END IF
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        DO iNode = 1, NodeMapping(iProc)%nSendUniqueNodes
          NodeMapping(iProc)%SendNodeSourceExt(iNode) = &
              NodeSourceExtTmp(NodeMapping(iProc)%SendNodeUniqueGlobalID(iNode))
        END DO
        CALL MPI_ISEND( NodeMapping(iProc)%SendNodeSourceExt(:)            &
            , NodeMapping(iProc)%nSendUniqueNodes           &
            , MPI_DOUBLE_PRECISION                                        &
            , iProc                                                       &
            , 666                                                         &
            , MPI_COMM_LEADERS_SHARED                                     &
            , SendRequest(iProc)                                          &
            , IERROR)
      END IF
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
        IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
        IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
    END DO

    DO iProc = 0, nLeaderGroupProcs - 1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        DO iNode = 1, NodeMapping(iProc)%nRecvUniqueNodes
          NodeSourceExtTmp(NodeMapping(iProc)%RecvNodeUniqueGlobalID(iNode)) = &
              NodeSourceExtTmp(NodeMapping(iProc)%RecvNodeUniqueGlobalID(iNode)) + &
              NodeMapping(iProc)%RecvNodeSourceExt(iNode)
        END DO
      END IF
    END DO
  END IF
  CALL BARRIER_AND_SYNC(NodeSourceExtTmp_Shared_Win,MPI_COMM_SHARED)
  firstNode = INT(REAL( myComputeNodeRank   *nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))+1
  lastNode  = INT(REAL((myComputeNodeRank+1)*nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))
#else
  firstNode = 1
  lastNode = nUniqueGlobalNodes
#endif

  ! Add NodeSourceExtTmp values of the last boundary interaction
  DO iNode=firstNode, lastNode
    NodeSourceExt(iNode) = NodeSourceExt(iNode) + NodeSourceExtTmp(iNode)
  END DO
#if USE_MPI
  CALL BARRIER_AND_SYNC(NodeSourceExt_Shared_Win,MPI_COMM_SHARED)
#else
  ! Reset local surface charge
  NodeSourceExtTmp = 0.
#endif

end if ! iter.NE.0


! Loop over all elements and store charge density values in equidistantly distributed nodes of PP_N=1
DO iElem=1,PP_nElems
  ! Copy values to equidistant distribution
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem+offsetElem)))
  NodeSourceExtEqui(1,0,0,0) = NodeSourceExt(NodeID(1))/NodeVolume(NodeID(1))
  NodeSourceExtEqui(1,1,0,0) = NodeSourceExt(NodeID(2))/NodeVolume(NodeID(2))
  NodeSourceExtEqui(1,1,1,0) = NodeSourceExt(NodeID(3))/NodeVolume(NodeID(3))
  NodeSourceExtEqui(1,0,1,0) = NodeSourceExt(NodeID(4))/NodeVolume(NodeID(4))
  NodeSourceExtEqui(1,0,0,1) = NodeSourceExt(NodeID(5))/NodeVolume(NodeID(5))
  NodeSourceExtEqui(1,1,0,1) = NodeSourceExt(NodeID(6))/NodeVolume(NodeID(6))
  NodeSourceExtEqui(1,1,1,1) = NodeSourceExt(NodeID(7))/NodeVolume(NodeID(7))
  NodeSourceExtEqui(1,0,1,1) = NodeSourceExt(NodeID(8))/NodeVolume(NodeID(8))

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
    !FutureTime=0.0 ! not required
    FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_NodeSourceExtGlobal',OutputTime))//'.h5'
    IF(MPIRoot) CALL GenerateFileSkeleton('NodeSourceExtGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)!,FutureTime)
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesNodeSourceExtGlobal',N_variables,StrArray=StrVarNames)
    CALL CloseDataFile()
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
USE MOD_Mesh_Vars          ,ONLY: offsetElem
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
USE MOD_Particle_Vars          ,ONLY: PDM, PEM, PartState, PartSpecies, PartMPF, usevMPF, nSpecies, VarTimeStep, Species
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY: UseDSMC, CollisMode,PartStateIntEn, DSMC, PolyatomMolDSMC, SpecDSMC, VibQuantsPar
USE MOD_DSMC_Vars              ,ONLY: ElectronicDistriPart, AmbipolElecVelo
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars          ,ONLY: velocityAtTime, velocityOutputAtTime
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem
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
INTEGER                        :: nVar
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
LOGICAL                        :: withDSMC=.FALSE.
INTEGER                        :: iElem_glob, iElem_loc
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER, ALLOCATABLE           :: VibQuantData(:,:)
REAL, ALLOCATABLE              :: ElecDistriData(:,:), AD_Data(:,:)
INTEGER,PARAMETER              :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER                        :: MaxQuantNum, iPolyatMole, iSpec, MaxElecQuant
! Integers of KIND=IK
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart
INTEGER(KIND=IK),ALLOCATABLE   :: PartInt(:,:)
INTEGER(KIND=IK)               :: locnPart_max
#if USE_MPI
INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
INTEGER(KIND=IK)               :: nParticles(0:nProcessors-1)
#endif
INTEGER                        :: ALLOCSTAT
!=============================================
! Required default values for KIND=IK
MaxQuantNum=-1
! Write properties -----------------------------------------------------------------------------------------------------------------
! Open dataset
!CALL H5DOPEN_F(File_ID,'DG_Solution',Dset_id,iError)

!!added for Evib, Erot writeout
withDSMC=useDSMC
IF (withDSMC) THEN
!IF (withDSMC) THEN
  IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0)) THEN !int ener + 3, vmpf +1
    PartDataSize=11
  ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. (DSMC%ElectronicModel.GT.0)) ) THEN !int ener + 2 and vmpf + 1
                                                                            ! or int energ +3 but no vmpf +1
    PartDataSize=10
  ELSE IF (CollisMode.GT.1) THEN
    PartDataSize=9 !int ener + 2
  ELSE IF (usevMPF) THEN
    PartDataSize=8 !+ 1 vmpf
  ELSE
    PartDataSize=7 !+ 0
  END IF
ELSE IF (usevMPF) THEN
  PartDataSize=8 !vmpf +1
ELSE
  PartDataSize=7
END IF

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  MaxQuantNum = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

IF (withDSMC.AND.(DSMC%ElectronicModel.EQ.2)) THEN
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

#if USE_MPI
sendbuf(1)=locnPart
recvbuf=0_IK
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart=recvbuf(1)
sendbuf(1)=recvbuf(1)+locnPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nGlobalNbrOfParticles=sendbuf(1)
GlobalNbrOfParticlesUpdated = .TRUE.
CALL MPI_GATHER(locnPart,1,MPI_INTEGER_INT_KIND,nParticles,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_WORLD,iError)
!IF (myRank.EQ.0) THEN
!  WRITE(*,*) 'PARTICLE-ELEMENT DISTRIBUTION'
!  WRITE(*,*) 'iProc, firstelemInd,   nElems,  locnPart,  totalnPart'
!  DO pcount=0,nProcessors-1
!    WRITE(*,'(I5,4I12)')pcount,offsetElemMPI(pcount),offsetElemMPI(pcount+1)-offsetElemMPI(pcount),&
!                       nParticles(pcount),SUM(nParticles(0:pcount))
!  END DO
!END IF
LOGWRITE(*,*)'offsetnPart,locnPart,nGlobalNbrOfParticles',offsetnPart,locnPart,nGlobalNbrOfParticles
CALL MPI_REDUCE(locnPart, locnPart_max, 1, MPI_INTEGER_INT_KIND, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
#else
offsetnPart=0_IK
nGlobalNbrOfParticles=locnPart
locnPart_max=locnPart
#endif
ALLOCATE(PartInt(offsetElem+1:offsetElem+PP_nElems,PartIntSize))
ALLOCATE(PartData(INT(PartDataSize,IK),offsetnPart+1_IK:offsetnPart+locnPart), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'ERROR in hdf5_output.f90: Cannot allocate PartData array for writing particle data to .h5!')
!!! Kleiner Hack von JN (Teil 1/2):

IF (.NOT.(useDSMC.OR.usevMPF)) THEN
  ALLOCATE(PEM%pStart(1:PP_nElems)           , &
           PEM%pNumber(1:PP_nElems)          , &
           PEM%pNext(1:PDM%maxParticleNumber), &
           PEM%pEnd(1:PP_nElems) )!            , &
           !PDM%nextUsedPosition(PDM%maxParticleNumber)  )
  useDSMC=.TRUE.
END IF
CALL UpdateNextFreePosition()
!!! Ende kleiner Hack von JN (Teil 1/2)
iPart=offsetnPart
DO iElem_loc=1,PP_nElems
  iElem_glob = iElem_loc + offsetElem
  PartInt(iElem_glob,1)=iPart
  IF (ALLOCATED(PEM%pNumber)) THEN
    nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
    PartInt(iElem_glob,2) = PartInt(iElem_glob,1) + INT(PEM%pNumber(iElem_loc),IK)
    pcount = PEM%pStart(iElem_loc)
    DO iPart=PartInt(iElem_glob,1)+1_IK,PartInt(iElem_glob,2)
      PartData(1,iPart)=PartState(1,pcount)
      PartData(2,iPart)=PartState(2,pcount)
      PartData(3,iPart)=PartState(3,pcount)
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
      IF (velocityOutputAtTime) THEN
        PartData(4,iPart)=velocityAtTime(1,pcount)
        PartData(5,iPart)=velocityAtTime(2,pcount)
        PartData(6,iPart)=velocityAtTime(3,pcount)
      ELSE
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
      PartData(4,iPart)=PartState(4,pcount)
      PartData(5,iPart)=PartState(5,pcount)
      PartData(6,iPart)=PartState(6,pcount)
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
      IF (withDSMC) THEN
      !IF (withDSMC) THEN
        IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
          PartData(8,iPart)=PartStateIntEn(1,pcount)
          PartData(9,iPart)=PartStateIntEn(2,pcount)
          PartData(10,iPart)=PartStateIntEn(3,pcount)
          PartData(11,iPart)=PartMPF(pcount)
        ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
          PartData(8,iPart)=PartStateIntEn(1,pcount)
          PartData(9,iPart)=PartStateIntEn(2,pcount)
          PartData(10,iPart)=PartMPF(pcount)
        ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
          PartData(8,iPart)=PartStateIntEn(1,pcount)
          PartData(9,iPart)=PartStateIntEn(2,pcount)
          PartData(10,iPart)=PartStateIntEn(3,pcount)
        ELSE IF (CollisMode.GT.1) THEN
          PartData(8,iPart)=PartStateIntEn(1,pcount)
          PartData(9,iPart)=PartStateIntEn(2,pcount)
        ELSE IF (usevMPF) THEN
          PartData(8,iPart)=PartMPF(pcount)
        END IF
      ELSE IF (usevMPF) THEN
          PartData(8,iPart)=PartMPF(pcount)
      END IF
      pcount = PEM%pNext(pcount)
    END DO
    iPart = PartInt(iElem_glob,2)
  ELSE
    CALL abort(&
    __STAMP__&
    , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
  END IF
  PartInt(iElem_glob,2)=iPart
END DO

nVar=2
ALLOCATE(StrVarNames(nVar))
StrVarNames(1)='FirstPartID'
StrVarNames(2)='LastPartID'
!CALL WriteAttributeToHDF5(File_ID,'VarNamesPartInt',2,StrArray=StrVarNames)

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesPartInt',nVar,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch=.TRUE.
  gatheredWrite=.FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems          => INT(nGlobalElems,IK)          ,&
      nVar                  => INT(nVar,IK)                  ,&
      PP_nElems             => INT(PP_nElems,IK)             ,&
      offsetElem            => INT(offsetElem,IK)            ,&
      PartDataSize          => INT(PartDataSize,IK)          )

  CALL GatheredWriteArray(FileName                         , create = .FALSE.            , &
                          DataSetName     = 'PartInt'      , rank   = 2                  , &
                          nValGlobal      = (/nGlobalElems , nVar/)                      , &
                          nVal            = (/PP_nElems    , nVar/)                      , &
                          offset          = (/offsetElem   , 0_IK/)                      , &
                          collective      = .TRUE.         , IntegerArray = PartInt)

  !  CALL WriteArrayToHDF5(DataSetName='PartInt', rank=2,&
  !                        nValGlobal=(/nGlobalElems, nVar/),&
  !                        nVal=      (/PP_nElems, nVar   /),&
  !                        offset=    (/offsetElem, 0_IK  /),&
  !                        collective=.TRUE., existing=.FALSE., IntegerArray=PartInt)
  !
  DEALLOCATE(StrVarNames)
  !CALL CloseDataFile()

  ALLOCATE(StrVarNames(PartDataSize))
  StrVarNames(1)='ParticlePositionX'
  StrVarNames(2)='ParticlePositionY'
  StrVarNames(3)='ParticlePositionZ'
  StrVarNames(4)='VelocityX'
  StrVarNames(5)='VelocityY'
  StrVarNames(6)='VelocityZ'
  StrVarNames(7)='Species'

  IF(withDSMC)THEN
    ! IF(withDSMC)THEN
    IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicModel.GT.0))THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
      StrVarNames(11)='MPF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='MPF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
    ELSE IF (CollisMode.GT.1) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
    ELSE IF (usevMPF) THEN
      StrVarNames( 8)='MPF'
    END IF
  ELSE IF (usevMPF) THEN
    StrVarNames( 8)='MPF'
  END IF

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',INT(PartDataSize,4),StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

  IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName='PartData'    , rank=2              , &
                            nValGlobal=(/PartDataSize , nGlobalNbrOfParticles /)       , &
                            nVal=      (/PartDataSize , locnPart   /)       , &
                            offset=    (/0_IK         , offsetnPart/)       , &
                            collective=.FALSE.        , RealArray=PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !locnPart_max.EQ.0
#if USE_MPI
  CALL DistributedWriteArray(FileName                      , &
                             DataSetName  = 'PartData'     , rank= 2        , &
                             nValGlobal   = (/PartDataSize , nGlobalNbrOfParticles /)  , &
                             nVal         = (/PartDataSize , locnPart   /)  , &
                             offset       = (/0_IK         , offsetnPart/)  , &
                             collective   = .FALSE.        , offSetDim= 2   , &
                             communicator = PartMPI%COMM   , RealArray= PartData)
  ! Output of the element-wise time step as a separate container in state file
  IF(VarTimeStep%UseDistribution) THEN
    CALL DistributedWriteArray(FileName , &
                              DataSetName = 'PartTimeStep'  , rank=2      , &
                              nValGlobal  = (/nGlobalElems  , 1_IK/)      , &
                              nVal        = (/PP_nElems     , 1_IK/)      , &
                              offset      = (/offsetElem    , 0_IK/)      , &
                              collective  =.FALSE.          , offSetDim=1 , &
                              communicator=PartMPI%COMM     , RealArray=VarTimeStep%ElemFac)
  END IF
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'     , rank = 2       , &
                        nValGlobal  = (/PartDataSize , nGlobalNbrOfParticles /)  , &
                        nVal        = (/PartDataSize , locnPart   /)  , &
                        offset      = (/0_IK         , offsetnPart/)  , &
                        collective  = .TRUE.         , RealArray = PartData)
    ! Output of the element-wise time step as a separate container in state file
  IF(VarTimeStep%UseDistribution) THEN
    CALL WriteArrayToHDF5(DataSetName = 'PartTimeStep'  , rank=2, &
                          nValGlobal  = (/nGlobalElems  , 1_IK/), &
                          nVal        = (/PP_nElems     , 1_IK/)   ,&
                          offset      = (/offsetElem    , 0_IK/)  ,&
                          collective  = .FALSE.         , RealArray=VarTimeStep%ElemFac)
  END IF
  CALL CloseDataFile()
#endif /*USE_MPI*/
END ASSOCIATE

DEALLOCATE(StrVarNames)
DEALLOCATE(PartData)

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  ALLOCATE(VibQuantData(MaxQuantNum,offsetnPart+1_IK:offsetnPart+locnPart))
  VibQuantData = 0
  !+1 is real number of necessary vib quants for the particle
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(iElem_glob,1)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
      PartInt(iElem_glob,2) = PartInt(iElem_glob,1) + INT(PEM%pNumber(iElem_loc),IK)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(iElem_glob,1)+1_IK,PartInt(iElem_glob,2)
        IF (SpecDSMC(PartSpecies(pcount))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(PartSpecies(pcount))%SpecToPolyArray
          VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,iPart) = &
            VibQuantsPar(pcount)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
        ELSE
           VibQuantData(:,iPart) = 0
        END IF
        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(iElem_glob,2)
    ELSE
      CALL abort(&
      __STAMP__&
      , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(iElem_glob,2)=iPart
  END DO

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (MaxQuantNum           => INT(MaxQuantNum,IK))
    IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file, required for (auto-)restart
      IF(MPIRoot)THEN ! only root writes the container
        CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
        CALL WriteArrayToHDF5(DataSetName='VibQuantData', rank=2              , &
                              nValGlobal=(/MaxQuantNum  , nGlobalNbrOfParticles /)       , &
                              nVal=      (/MaxQuantNum  , locnPart   /)       , &
                              offset=    (/0_IK         , offsetnPart/)       , &
                              collective=.FALSE.        , IntegerArray_i4=VibQuantData)
        CALL CloseDataFile()
      END IF !MPIRoot
    END IF !locnPart_max.EQ.0
#if USE_MPI
    CALL DistributedWriteArray(FileName , &
                              DataSetName ='VibQuantData', rank=2           , &
                              nValGlobal  =(/MaxQuantNum , nGlobalNbrOfParticles  /)   , &
                              nVal        =(/MaxQuantNum , locnPart    /)   , &
                              offset      =(/0_IK        , offsetnPart /)   , &
                              collective  =.FALSE.       , offSetDim=2      , &
                              communicator=PartMPI%COMM  , IntegerArray_i4=VibQuantData)
    DEALLOCATE(VibQuantData)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArrayToHDF5(DataSetName = 'VibQuantData' , rank = 2             , &
                          nValGlobal  = (/ MaxQuantNum , nGlobalNbrOfParticles   /)      , &
                          nVal        = (/ MaxQuantNum , locnPart     /)      , &
                          offset      = (/ 0_IK        , offsetnPart  /)      , &
                          collective  = .TRUE.         , IntegerArray_i4 = VibQuantData)
    DEALLOCATE(VibQuantData)
    CALL CloseDataFile()
#endif /*USE_MPI*/
  END ASSOCIATE
END IF

IF (withDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
  ALLOCATE(ElecDistriData(MaxElecQuant,offsetnPart+1_IK:offsetnPart+locnPart))
  ElecDistriData = 0
  !+1 is real number of necessary vib quants for the particle
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(iElem_glob,1)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
      PartInt(iElem_glob,2) = PartInt(iElem_glob,1) + INT(PEM%pNumber(iElem_loc),IK)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(iElem_glob,1)+1_IK,PartInt(iElem_glob,2)
        IF (.NOT.((SpecDSMC(PartSpecies(pcount))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(pcount))%FullyIonized)) THEN
          ElecDistriData(1:SpecDSMC(PartSpecies(pcount))%MaxElecQuant,iPart) = &
            ElectronicDistriPart(pcount)%DistriFunc(1:SpecDSMC(PartSpecies(pcount))%MaxElecQuant)
        ELSE
           ElecDistriData(:,iPart) = 0
        END IF
        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(iElem_glob,2)
    ELSE
      CALL abort(&
      __STAMP__&
      , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(iElem_glob,2)=iPart
  END DO

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (MaxElecQuant          => INT(MaxElecQuant,IK))
    IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file, required for (auto-)restart
      IF(MPIRoot)THEN ! only root writes the container
        CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
        CALL WriteArrayToHDF5(DataSetName='ElecDistriData', rank=2              , &
                              nValGlobal=(/MaxElecQuant   , nGlobalNbrOfParticles /)       , &
                              nVal=      (/MaxElecQuant   , locnPart   /)       , &
                              offset=    (/0_IK           , offsetnPart/)       , &
                              collective=.FALSE.          , RealArray=ElecDistriData)
        CALL CloseDataFile()
      END IF !MPIRoot
    END IF !locnPart_max.EQ.0
#if USE_MPI
    CALL DistributedWriteArray(FileName , &
                              DataSetName ='ElecDistriData', rank=2           , &
                              nValGlobal  =(/MaxElecQuant  , nGlobalNbrOfParticles  /)   , &
                              nVal        =(/MaxElecQuant  , locnPart    /)   , &
                              offset      =(/0_IK          , offsetnPart /)   , &
                              collective  =.FALSE.         , offSetDim=2      , &
                              communicator=PartMPI%COMM    , RealArray=ElecDistriData)
    DEALLOCATE(ElecDistriData)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArrayToHDF5(DataSetName = 'ElecDistriData' , rank = 2             , &
                          nValGlobal  = (/ MaxElecQuant  , nGlobalNbrOfParticles   /)      , &
                          nVal        = (/ MaxElecQuant  , locnPart     /)      , &
                          offset      = (/ 0_IK          , offsetnPart  /)      , &
                          collective  = .TRUE.           , RealArray = ElecDistriData)
    DEALLOCATE(ElecDistriData)
    CALL CloseDataFile()
#endif /*USE_MPI*/
  END ASSOCIATE
END IF

IF (withDSMC.AND.DSMC%DoAmbipolarDiff) THEN
  ALLOCATE(AD_Data(3,offsetnPart+1_IK:offsetnPart+locnPart))
  AD_Data = 0.0
  !+1 is real number of necessary vib quants for the particle
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(iElem_glob,1)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
      PartInt(iElem_glob,2) = PartInt(iElem_glob,1) + INT(PEM%pNumber(iElem_loc),IK)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(iElem_glob,1)+1_IK,PartInt(iElem_glob,2)
        IF (Species(PartSpecies(pcount))%ChargeIC.GT.0.0) THEN
          AD_Data(1:3,iPart) = AmbipolElecVelo(pcount)%ElecVelo(1:3)
        ELSE
          AD_Data(1:3,iPart) = 0
        END IF
        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(iElem_glob,2)
    ELSE
      CALL abort(&
      __STAMP__&
      , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(iElem_glob,2)=iPart
  END DO

    IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file, required for (auto-)restart
      IF(MPIRoot)THEN ! only root writes the container
        CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
        CALL WriteArrayToHDF5(DataSetName='ADVeloData' , rank=2              , &
                              nValGlobal=(/3_IK        , nGlobalNbrOfParticles /)       , &
                              nVal=      (/3_IK        , locnPart   /)       , &
                              offset=    (/0_IK        , offsetnPart/)       , &
                              collective=.FALSE.       , RealArray=AD_Data)
        CALL CloseDataFile()
      END IF !MPIRoot
    END IF !locnPart_max.EQ.0
#if USE_MPI
    CALL DistributedWriteArray(FileName , &
                              DataSetName ='ADVeloData'  , rank=2           , &
                              nValGlobal  =(/3_IK        , nGlobalNbrOfParticles  /)   , &
                              nVal        =(/3_IK        , locnPart    /)   , &
                              offset      =(/0_IK        , offsetnPart /)   , &
                              collective  =.FALSE.       , offSetDim=2      , &
                              communicator=PartMPI%COMM  , RealArray=AD_Data)
    DEALLOCATE(AD_Data)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArrayToHDF5(DataSetName = 'ADVeloData'   , rank = 2             , &
                          nValGlobal  = (/ 3_IK        , nGlobalNbrOfParticles   /)      , &
                          nVal        = (/ 3_IK        , locnPart     /)      , &
                          offset      = (/ 0_IK        , offsetnPart  /)      , &
                          collective  = .TRUE.         , RealArray = AD_Data)
    DEALLOCATE(AD_Data)
    CALL CloseDataFile()
#endif /*USE_MPI*/
END IF

DEALLOCATE(PartInt)
! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

!!! Kleiner Hack von JN (Teil 2/2):
useDSMC=withDSMC
IF (.NOT.(useDSMC.OR.usevMPF)) THEN
  DEALLOCATE(PEM%pStart , &
             PEM%pNumber, &
             PEM%pNext  , &
             PEM%pEnd   )!, &
             !PDM%nextUsedPosition  )
END IF
!!! Ende kleiner Hack von JN (Teil 2/2)


END SUBROUTINE WriteParticleToHDF5


SUBROUTINE WriteBoundaryParticleToHDF5(MeshFileName,OutputTime,PreviousTime)
!===================================================================================================================================
! Write data of impacting particles on specific boundary conditions of .h5 file (position, velocity, species ID, kinetic energy [eV],
! macro particle factor, time of impact, impact obliqueness angle)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ElementaryCharge
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems, offsetElem
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,nVarPartStateBoundary
USE MOD_Equation_Vars          ,ONLY: StrVarNames
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart2
USE MOD_TimeDisc_Vars          ,ONLY: iter
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
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
#if USE_MPI
INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
INTEGER(KIND=IK)               :: nParticles(0:nProcessors-1)
#endif
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,nPart_glob
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER(KIND=IK)               :: locnPart_max
CHARACTER(LEN=255)             :: FileName,PreviousFileName
REAL                           :: PreviousTime_loc
INTEGER                        :: ALLOCSTAT
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
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

! 3xPos [m], 3xvelo [m/s], species [-]
PartDataSize = 7
! Kinetic energy [eV]
PartDataSize = PartDataSize + 1
! MPF [-]
PartDataSize = PartDataSize + 1
! time [s]
PartDataSize = PartDataSize + 1
! Impact obliqueness angle [degree]
PartDataSize = PartDataSize + 1
! iBC [-]
PartDataSize = PartDataSize + 1

! Set number of local particles
locnPart = INT(PartStateBoundaryVecLength,IK)

#if USE_MPI
sendbuf(1)  = locnPart
recvbuf     = 0_IK
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart = recvbuf(1)
sendbuf(1)  = recvbuf(1)+locnPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nPart_glob  = sendbuf(1)
CALL MPI_GATHER(locnPart,1,MPI_INTEGER_INT_KIND,nParticles,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_WORLD,iError)
LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
CALL MPI_REDUCE(locnPart, locnPart_max, 1, MPI_INTEGER_INT_KIND, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
#else
offsetnPart  = 0_IK
nPart_glob   = locnPart
locnPart_max = locnPart
#endif
ALLOCATE(PartData(INT(PartDataSize,IK),offsetnPart+1_IK:offsetnPart+locnPart), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'ERROR in hdf5_output.f90: Cannot allocate PartData array for writing boundary particle data to .h5!')

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
  PartData(8,iPart)=CalcEkinPart2(PartStateBoundary(4:6,pcount),INT(ABS(PartStateBoundary(7,pcount))),1.0) / ElementaryCharge

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
      PartDataSize => INT(PartDataSize,IK) )

  ALLOCATE(StrVarNames2(PartDataSize))
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
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',INT(PartDataSize,4),StrArray=StrVarNames2)
    CALL CloseDataFile()
  END IF

  IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName='PartData'     , rank=2              , &
                            nValGlobal=(/ PartDataSize , nPart_glob  /)      , &
                            nVal=      (/ PartDataSize , locnPart    /)      , &
                            offset=    (/ 0_IK         , offsetnPart /)      , &
                            collective=.FALSE.         , RealArray=PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !locnPart_max.EQ.0
#if USE_MPI
  CALL DistributedWriteArray(FileName                       , &
                             DataSetName  = 'PartData'      , rank = 2              , &
                             nValGlobal   = (/ PartDataSize , nPart_glob  /)        , &
                             nVal         = (/ PartDataSize , locnPart    /)        , &
                             offset       = (/ 0_IK         , offsetnPart /)        , &
                             collective   = .FALSE.         , offSetDim = 2         , &
                             communicator = PartMPI%COMM    , RealArray = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'      , rank = 2              , &
                        nValGlobal  = (/ PartDataSize , nPart_glob  /)        , &
                        nVal        = (/ PartDataSize , locnPart    /)        , &
                        offset      = (/ 0_IK         , offsetnPart /)        , &
                        collective  = .TRUE.          , RealArray = PartData)
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
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
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
#if USE_MPI
INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
INTEGER(KIND=IK)               :: nParticles(0:nProcessors-1)
#endif
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,nPart_glob
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER(KIND=IK)               :: locnPart_max
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
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

! Set number of local particles
locnPart = INT(PartStateLostVecLength,IK)

#if USE_MPI
sendbuf(1)  = locnPart
recvbuf     = 0_IK
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart = recvbuf(1)
sendbuf(1)  = recvbuf(1)+locnPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nPart_glob  = sendbuf(1)
CALL MPI_GATHER(locnPart,1,MPI_INTEGER_INT_KIND,nParticles,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_WORLD,iError)
LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
CALL MPI_REDUCE(locnPart, locnPart_max, 1, MPI_INTEGER_INT_KIND, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
#else
offsetnPart  = 0_IK
nPart_glob   = locnPart
locnPart_max = locnPart
#endif

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

  IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName='PartData'         , rank=2              , &
                            nValGlobal=(/ PartLostDataSize , nPart_glob  /)      , &
                            nVal=      (/ PartLostDataSize , locnPart    /)      , &
                            offset=    (/ 0_IK             , offsetnPart /)      , &
                            collective=.FALSE.             , RealArray=PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !locnPart_max.EQ.0
#if USE_MPI
  CALL DistributedWriteArray(FileName                           , &
                             DataSetName  = 'PartData'          , rank = 2              , &
                             nValGlobal   = (/ PartLostDataSize , nPart_glob  /)        , &
                             nVal         = (/ PartLostDataSize , locnPart    /)        , &
                             offset       = (/ 0_IK             , offsetnPart /)        , &
                             collective   = .FALSE.             , offSetDim = 2         , &
                             communicator = PartMPI%COMM        , RealArray = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'      , rank = 2              , &
                        nValGlobal  = (/ PartLostDataSize , nPart_glob  /)        , &
                        nVal        = (/ PartLostDataSize , locnPart    /)        , &
                        offset      = (/ 0_IK             , offsetnPart /)        , &
                        collective  = .TRUE.              , RealArray = PartData)
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
USE MOD_Particle_Vars          ,ONLY: nSpecies
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
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255)             :: SpecID
INTEGER                        :: nVar, nVarTotal
INTEGER                        :: ElemID,iVar,iSpec,SampleElemID
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
  StrVarNames(iVar+5) = 'Spec'//TRIM(SpecID)//'-PumpPressure'
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
DO iSpec = 1, nSpecies
  DO SampleElemID = 1,AdaptBCSampleElemNum
    ElemID = AdaptBCMapSampleToElem(SampleElemID)
    AdaptiveData(iVar:iVar-1+nVar,ElemID) = AdaptBCMacroVal(1:7,SampleElemID,iSpec)
  END DO
  iVar = iVar + nVar
END DO

WRITE(H5_Name,'(A)') 'AdaptiveInfo'
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      nElems          => INT(nElems,IK)          ,&
      nVarTotal       => INT(nVarTotal,IK)       ,&
      nSpecies        => INT(nSpecies,IK)        ,&
      offsetElem      => INT(offsetElem,IK)      )
  CALL WriteArrayToHDF5(DataSetName = H5_Name, rank = 2            , &
                        nValGlobal  = (/nVarTotal , nGlobalElems/) , &
                        nVal        = (/nVarTotal , nElems      /) , &
                        offset      = (/0_IK      , offsetElem  /) , &
                        collective  = .false., RealArray = AdaptiveData)
END ASSOCIATE
CALL CloseDataFile()
SDEALLOCATE(StrVarNames)

IF(AdaptBCTruncAverage) CALL WriteAdaptiveRunningAverageToHDF5(FileName)

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
USE MOD_Particle_Vars          ,ONLY: nSpecies
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

! Store the position in the array for early restarts
IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  IF(INT(iter,4)+AdaptBCSampIterReadIn.LT.AdaptBCSampIter) THEN
    CALL WriteAttributeToHDF5(File_ID,'AdaptBCSampIter',1,IntegerScalar=INT(iter,4)+AdaptBCSampIterReadIn)
  ELSE
    CALL WriteAttributeToHDF5(File_ID,'AdaptBCSampIter',1,IntegerScalar=AdaptBCSampIter)
  END IF
  CALL CloseDataFile()
END IF

CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
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
                        collective  = .false. , RealArray = AdaptBCAverage)
  CALL WriteArrayToHDF5(DataSetName = 'AdaptiveRunningAverageIndex' , rank = 1                   , &
                        nValGlobal  = (/AdaptBCSampleElemNumGlobal/) , &
                        nVal        = (/AdaptBCSampleElemNum      /) , &
                        offset      = (/offSetElemAdaptBCSample  /) , &
                        collective  = .false. , IntegerArray_i4 = AdaptBCAverageIndex)
END ASSOCIATE
CALL CloseDataFile()
DEALLOCATE(AdaptBCAverageIndex)

END SUBROUTINE WriteAdaptiveRunningAverageToHDF5


SUBROUTINE WriteAdaptiveWallTempToHDF5(FileName)
!===================================================================================================================================
!> Output of the adaptive cell-local wall temperature and the corresponding global side index
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, nSurfTotalSides, nComputeNodeSurfSides, offsetComputeNodeSurfSide
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
IF (nSurfTotalSides      .EQ.0) RETURN
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
      nGlobalSides         => INT(nSurfTotalSides,IK)           , &
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
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nGlobalElems    => INT(nGlobalElems,IK)    ,&
          nElems          => INT(nElems,IK)          ,&
          nSpecies        => INT(nSpecies,IK)        ,&
          offsetElem      => INT(offsetElem,IK)      )
      CALL WriteArrayToHDF5(DataSetName = H5_Name , rank = 2        , &
                            nValGlobal  = (/nGlobalElems, nSpecies/) , &
                            nVal        = (/nElems      , nSpecies/) , &
                            offset      = (/offsetElem  ,0_IK     /) , &
                            collective  = .false. , RealArray = VarVibRelaxProb%ProbVibAv)
    END ASSOCIATE
    CALL CloseDataFile()
    SDEALLOCATE(StrVarNames)
  ELSE ! DSMC%VibRelaxProb < 2.0
#if USE_MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
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
USE MOD_DSMC_Vars     ,ONLY: UseDSMC, CollisMode, DSMC, PolyatomMolDSMC, SpecDSMC
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting, ClonedParticles
USE MOD_PARTICLE_Vars ,ONLY: nSpecies, usevMPF, Species
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
#if USE_MPI
INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
#endif
INTEGER                        :: pcount, iDelay, iElem_glob
LOGICAL                        :: withDSMC=.FALSE.
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,nPart_glob
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER, ALLOCATABLE           :: VibQuantData(:,:)
REAL, ALLOCATABLE              :: ElecDistriData(:,:), AD_Data(:,:)
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER                        :: MaxQuantNum, iPolyatMole, iSpec, tempDelay, MaxElecQuant
!-----------------------------------------------------------------------------------------------------------------------------------
!!added for Evib, Erot writeout
withDSMC=useDSMC
IF (withDSMC) THEN
  IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0)) THEN !int ener + 3, vmpf +1
    PartDataSize=13
  ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. (DSMC%ElectronicModel.GT.0)) ) THEN !int ener + 2 and vmpf + 1
                                                                            ! or int energ +3 but no vmpf +1
    PartDataSize=12
  ELSE IF (CollisMode.GT.1) THEN
    PartDataSize=11 !int ener + 2
  ELSE IF (usevMPF) THEN
    PartDataSize=10 !+ 1 vmpf
  ELSE
    PartDataSize=9 !+ 0
  END IF
ELSE IF (usevMPF) THEN
  PartDataSize=10 !vmpf +1
ELSE
  PartDataSize=9
END IF

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  MaxQuantNum = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

IF (withDSMC.AND.(DSMC%ElectronicModel.EQ.2)) THEN
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
  CALL abort(__STAMP__,&
              'RadialWeighting: CloneMode is not supported!')
END SELECT

DO pcount = 0,tempDelay
    locnPart = locnPart + RadialWeighting%ClonePartNum(pcount)
END DO

#if USE_MPI
sendbuf(1)=locnPart
recvbuf=0
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart=recvbuf(1)
sendbuf(1)=recvbuf(1)+locnPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nPart_glob=sendbuf(1)
#else
offsetnPart=0
nPart_glob=locnPart
#endif
ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  ALLOCATE(VibQuantData(MaxQuantNum,offsetnPart+1:offsetnPart+locnPart))
  VibQuantData = 0
  !+1 is real number of necessary vib quants for the particle
END IF
IF (withDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
  ALLOCATE(ElecDistriData(MaxElecQuant,offsetnPart+1_IK:offsetnPart+locnPart))
  ElecDistriData = 0
  !+1 is real number of necessary vib quants for the particle
END IF
IF (withDSMC.AND.DSMC%DoAmbipolarDiff)  THEN
  ALLOCATE(AD_Data(3,offsetnPart+1_IK:offsetnPart+locnPart))
  AD_Data = 0
  !+1 is real number of necessary vib quants for the particle
END IF
iPart=offsetnPart
DO iDelay=0,tempDelay
  DO pcount = 1, RadialWeighting%ClonePartNum(iDelay)
    iElem_glob = ClonedParticles(pcount,iDelay)%Element
    iPart = iPart + 1
    PartData(1,iPart)=ClonedParticles(pcount,iDelay)%PartState(1)
    PartData(2,iPart)=ClonedParticles(pcount,iDelay)%PartState(2)
    PartData(3,iPart)=ClonedParticles(pcount,iDelay)%PartState(3)
    PartData(4,iPart)=ClonedParticles(pcount,iDelay)%PartState(4)
    PartData(5,iPart)=ClonedParticles(pcount,iDelay)%PartState(5)
    PartData(6,iPart)=ClonedParticles(pcount,iDelay)%PartState(6)
    PartData(7,iPart)=REAL(ClonedParticles(pcount,iDelay)%Species)
    PartData(8,iPart)=REAL(iElem_glob)
    PartData(9,iPart)=REAL(iDelay)
    IF (withDSMC) THEN
      IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
        PartData(10,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(11,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
        PartData(12,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(3)
        PartData(13,iPart)=ClonedParticles(pcount,iDelay)%WeightingFactor
      ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
        PartData(10,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(11,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
        PartData(12,iPart)=ClonedParticles(pcount,iDelay)%WeightingFactor
      ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
        PartData(10,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(11,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
        PartData(12,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(3)
      ELSE IF (CollisMode.GT.1) THEN
        PartData(10,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(11,iPart)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
      ELSE IF (usevMPF) THEN
        PartData(10,iPart)=ClonedParticles(pcount,iDelay)%WeightingFactor
      END IF
    ELSE IF (usevMPF) THEN
        PartData(10,iPart)=ClonedParticles(pcount,iDelay)%WeightingFactor
    END IF
    IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
      IF (SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%SpecToPolyArray
        VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,iPart) = &
          ClonedParticles(pcount,iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
      ELSE
          VibQuantData(:,iPart) = 0
      END IF
    END IF
    IF (withDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
      IF (.NOT.((SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%InterID.EQ.4) &
          .OR.SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%FullyIonized)) THEN
          ElecDistriData(1:SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%MaxElecQuant,iPart) = &
            ClonedParticles(pcount,iDelay)%DistriFunc(1:SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%MaxElecQuant)
      END IF
    END IF
    IF (withDSMC.AND.DSMC%DoAmbipolarDiff)  THEN
      IF (Species(ClonedParticles(pcount,iDelay)%Species)%ChargeIC.GT.0.0) THEN
        AD_Data(1:3,iPart) = ClonedParticles(pcount,iDelay)%AmbiPolVelo(1:3)
      END IF
    END IF
  END DO
END DO

ALLOCATE(StrVarNames(PartDataSize))
StrVarNames(1)='ParticlePositionX'
StrVarNames(2)='ParticlePositionY'
StrVarNames(3)='ParticlePositionZ'
StrVarNames(4)='VelocityX'
StrVarNames(5)='VelocityY'
StrVarNames(6)='VelocityZ'
StrVarNames(7)='Species'
StrVarNames(8)='Element'
StrVarNames(9)='CloneDelay'

IF(withDSMC)THEN
  IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicModel.GT.0))THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
    StrVarNames(12)='Electronic'
    StrVarNames(13)='MPF'
  ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
    StrVarNames(12)='MPF'
  ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
    StrVarNames(12)='Electronic'
  ELSE IF (CollisMode.GT.1) THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
  ELSE IF (usevMPF) THEN
    StrVarNames(10)='MPF'
  END IF
END IF

#if USE_MPI
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif
CALL WriteAttributeToHDF5(File_ID,'VarNamesParticleClones',PartDataSize,StrArray=StrVarNames)

ASSOCIATE (&
      nPart_glob      => INT(nPart_glob,IK)    ,&
      offsetnPart     => INT(offsetnPart,IK)   ,&
      MaxQuantNum     => INT(MaxQuantNum,IK)   ,&
      MaxElecQuant    => INT(MaxElecQuant,IK)   ,&
      PartDataSize    => INT(PartDataSize,IK)  )
CALL WriteArrayToHDF5(DataSetName='CloneData'   , rank=2         , &
                      nValGlobal=(/PartDataSize , nPart_glob /)  , &
                      nVal=      (/PartDataSize , locnPart   /)  , &
                      offset=    (/0_IK         , offsetnPart/)  , &
                      collective=.FALSE.        , RealArray=PartData)
IF (withDSMC) THEN
  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    CALL WriteArrayToHDF5(DataSetName='CloneVibQuantData' , rank=2              , &
                          nValGlobal=(/MaxQuantNum        , nPart_glob     /)   , &
                          nVal=      (/MaxQuantNum        , locnPart       /)   , &
                          offset=    (/0_IK               , offsetnPart    /)   , &
                          collective=.FALSE.              , IntegerArray_i4=VibQuantData)
    DEALLOCATE(VibQuantData)
  END IF
  IF (DSMC%ElectronicModel.EQ.2) THEN
    CALL WriteArrayToHDF5(DataSetName='CloneElecDistriData' , rank=2              , &
                          nValGlobal=(/MaxElecQuant       , nPart_glob     /)   , &
                          nVal=      (/MaxElecQuant        , locnPart       /)   , &
                          offset=    (/0_IK               , offsetnPart    /)   , &
                          collective=.FALSE.              , RealArray=ElecDistriData)
    DEALLOCATE(ElecDistriData)
  END IF
  IF (DSMC%DoAmbipolarDiff) THEN
    CALL WriteArrayToHDF5(DataSetName='CloneADVeloData'   , rank=2              , &
                          nValGlobal=(/3_IK               , nPart_glob     /)   , &
                          nVal=      (/3_IK               , locnPart       /)   , &
                          offset=    (/0_IK               , offsetnPart    /)   , &
                          collective=.FALSE.              , RealArray=AD_Data)
    DEALLOCATE(AD_Data)
  END IF
END IF
END ASSOCIATE

CALL CloseDataFile()

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
#if USE_MPI
REAL                           :: StartT,EndT
#endif /*USE_MPI*/
INTEGER                        :: iElem,i,j,k
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE PIC EM-FIELD TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#endif /*USE_MPI*/

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
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  outputVars   => INT(outputVars,IK)   ,&
  N            => INT(PP_N,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName='DG_Solution'    , rank=5 , &
                      nValGlobal=(/outputVars , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/outputVars , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
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
USE MOD_Particle_MPI_Vars ,ONLY: PartMPI
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
IF(.NOT.PartMPI%MPIRoot) RETURN

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
                               collective  = .false. , IntegerArray = NeutralizationBalanceTmp)
       END ASSOCIATE
       CALL CloseDataFile()

     END SELECT
  END DO  ! iInit
END DO  ! iSpec=1,nSpecies

END SUBROUTINE WriteEmissionVariablesToHDF5


#endif /*defined(PARTICLES)*/
END MODULE MOD_HDF5_Output_Particles
