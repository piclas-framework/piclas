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

MODULE MOD_HDF5_Input_Particles
#if defined(PARTICLES)
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ReadNodeSourceExtFromHDF5
PUBLIC :: ReadEmissionVariablesFromHDF5
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Read NodeSourceExt from h5 file, which is stored as DG solution type field 'DG_SourceExt'.
!> Map this solution to equidistant-node polynomial (NodeTypeVISU with N=1) and then map the solution to the global nodes
!> 'NodeSourceExt'.
!===================================================================================================================================
SUBROUTINE ReadNodeSourceExtFromHDF5()
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Dielectric_Vars        ,ONLY: DoDielectric
USE MOD_HDF5_Input             ,ONLY: ReadArray
USE MOD_HDF5_Input             ,ONLY: File_ID,DatasetExists
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeVISU,NodeType
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Mesh_Vars              ,ONLY: Vdm_N_EQ,offsetElem
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID,GetGlobalElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,nUniqueGlobalNodes!,NodeToElemMapping,NodeToElemInfo
USE MOD_PICDepo_Vars           ,ONLY: NodeSourceExt,NodeVolume,DoDeposition
USE MOD_Restart_Vars           ,ONLY: N_Restart
!#if USE_MPI
!USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
!USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
!USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
!#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: U_local(1,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems)
LOGICAL                            :: DG_SourceExtExists
REAL                               :: NodeSourceExtEqui(1,0:1,0:1,0:1),NodeVol(1:8)
INTEGER(KIND=IK)                   :: OffsetElemTmp,PP_nElemsTmp,N_RestartTmp
INTEGER                            :: iElem!,CNElemID
INTEGER                            :: NodeID(1:8)!,firstNode,lastNode,firstGlobalElemID(1:8),iNode
!===================================================================================================================================
IF(.NOT.DoDielectric) RETURN

! This array is not allocated when DoDeposition=F, however, the previously calculated surface charge might still be required in the
! future, when DoDeposition is activated again. Therefore, read the old data and store in the new state file.
IF(.NOT.DoDeposition) THEN
  ALLOCATE(NodeSourceExt(1:nUniqueGlobalNodes))
  NodeSourceExt = 0.
END IF

! Temp. vars for integer KIND=8 possibility
OffsetElemTmp = INT(OffsetElem,IK)
PP_nElemsTmp  = INT(PP_nElems,IK)
N_RestartTmp  = INT(N_Restart,IK)

CALL DatasetExists(File_ID,'DG_SourceExt',DG_SourceExtExists)

IF(DG_SourceExtExists)THEN

  ! We need to interpolate the solution to the new computational grid
  CALL ReadArray('DG_SourceExt',5,(/1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
      OffsetElemTmp,5,RealArray=U_local)

  ! Allocate and determine Vandermonde mapping from NodeType to equidistant (visu) node set
  ALLOCATE(Vdm_N_EQ(0:1,0:N_Restart))
  CALL GetVandermonde(N_Restart, NodeType, 1, NodeTypeVISU, Vdm_N_EQ, modal=.FALSE.)

!  #if USE_MPI
!    firstNode = INT(REAL( myComputeNodeRank   )*REAL(nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))+1
!    lastNode  = INT(REAL((myComputeNodeRank+1))*REAL(nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))
!  #else
!    firstNode = 1
!    lastNode = nUniqueGlobalNodes
!  #endif
  NodeVol = 1. ! default
  DO iElem =1, PP_nElems
    ! Map G/GL (current node type) to equidistant distribution
    CALL ChangeBasis3D(1, N_Restart, 1, Vdm_N_EQ, U_local(:,:,:,:,iElem),NodeSourceExtEqui(:,:,:,:))

    ! Map the solution to the global nodes 'NodeSourceExt' and apply the volumes (charge density -> charge)
    ! Map non-unique to unique node ID
    NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem+offsetElem)))
    !DO iNode = 1, 8
    !  firstGlobalElemID(iNode) = GetGlobalElemID(NodeToElemInfo(NodeToElemMapping(1,NodeID(iNode)) + 1))
    !END DO ! I = 1, 8

    ! method 1: Only change the nodes which are assigned to the proc
    !IF(firstNode.LE.NodeID(1).AND.NodeID(1).LE.lastNode) NodeSourceExt(NodeID(1)) = NodeSourceExtEqui(1,0,0,0) * NodeVolume(NodeID(1))
    !IF(firstNode.LE.NodeID(2).AND.NodeID(2).LE.lastNode) NodeSourceExt(NodeID(2)) = NodeSourceExtEqui(1,1,0,0) * NodeVolume(NodeID(2))
    !IF(firstNode.LE.NodeID(3).AND.NodeID(3).LE.lastNode) NodeSourceExt(NodeID(3)) = NodeSourceExtEqui(1,1,1,0) * NodeVolume(NodeID(3))
    !IF(firstNode.LE.NodeID(4).AND.NodeID(4).LE.lastNode) NodeSourceExt(NodeID(4)) = NodeSourceExtEqui(1,0,1,0) * NodeVolume(NodeID(4))
    !IF(firstNode.LE.NodeID(5).AND.NodeID(5).LE.lastNode) NodeSourceExt(NodeID(5)) = NodeSourceExtEqui(1,0,0,1) * NodeVolume(NodeID(5))
    !IF(firstNode.LE.NodeID(6).AND.NodeID(6).LE.lastNode) NodeSourceExt(NodeID(6)) = NodeSourceExtEqui(1,1,0,1) * NodeVolume(NodeID(6))
    !IF(firstNode.LE.NodeID(7).AND.NodeID(7).LE.lastNode) NodeSourceExt(NodeID(7)) = NodeSourceExtEqui(1,1,1,1) * NodeVolume(NodeID(7))
    !IF(firstNode.LE.NodeID(8).AND.NodeID(8).LE.lastNode) NodeSourceExt(NodeID(8)) = NodeSourceExtEqui(1,0,1,1) * NodeVolume(NodeID(8))

    ! method 2: change any node
    ! this can lead to a race condition when two or more procs write to the same position in the
    ! shared memory array NodeSourceExt
    ! Fix: try MPI_FETCH_AND_OP(...)
    !NodeSourceExt(NodeID(1)) = NodeSourceExtEqui(1,0,0,0) * NodeVolume(NodeID(1))
    !NodeSourceExt(NodeID(2)) = NodeSourceExtEqui(1,1,0,0) * NodeVolume(NodeID(2))
    !NodeSourceExt(NodeID(3)) = NodeSourceExtEqui(1,1,1,0) * NodeVolume(NodeID(3))
    !NodeSourceExt(NodeID(4)) = NodeSourceExtEqui(1,0,1,0) * NodeVolume(NodeID(4))
    !NodeSourceExt(NodeID(5)) = NodeSourceExtEqui(1,0,0,1) * NodeVolume(NodeID(5))
    !NodeSourceExt(NodeID(6)) = NodeSourceExtEqui(1,1,0,1) * NodeVolume(NodeID(6))
    !NodeSourceExt(NodeID(7)) = NodeSourceExtEqui(1,1,1,1) * NodeVolume(NodeID(7))
    !NodeSourceExt(NodeID(8)) = NodeSourceExtEqui(1,0,1,1) * NodeVolume(NodeID(8))

    ! method 3: check node ID
    ! Compare global Elem ID and only write the data for the first element
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(1)) NodeSourceExt(NodeID(1)) = NodeSourceExtEqui(1,0,0,0) * NodeVolume(NodeID(1))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(2)) NodeSourceExt(NodeID(2)) = NodeSourceExtEqui(1,1,0,0) * NodeVolume(NodeID(2))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(3)) NodeSourceExt(NodeID(3)) = NodeSourceExtEqui(1,1,1,0) * NodeVolume(NodeID(3))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(4)) NodeSourceExt(NodeID(4)) = NodeSourceExtEqui(1,0,1,0) * NodeVolume(NodeID(4))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(5)) NodeSourceExt(NodeID(5)) = NodeSourceExtEqui(1,0,0,1) * NodeVolume(NodeID(5))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(6)) NodeSourceExt(NodeID(6)) = NodeSourceExtEqui(1,1,0,1) * NodeVolume(NodeID(6))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(7)) NodeSourceExt(NodeID(7)) = NodeSourceExtEqui(1,1,1,1) * NodeVolume(NodeID(7))
    !IF(iElem+OffsetElem.EQ.firstGlobalElemID(8)) NodeSourceExt(NodeID(8)) = NodeSourceExtEqui(1,0,1,1) * NodeVolume(NodeID(8))

    ! method 4: all nodes are local (removed shared memory)
    IF(DoDeposition) NodeVol = NodeVolume(NodeID(1:8)) ! this is only required when deposition is activated
    NodeSourceExt(NodeID(1)) = NodeSourceExtEqui(1,0,0,0) * NodeVol(1)
    NodeSourceExt(NodeID(2)) = NodeSourceExtEqui(1,1,0,0) * NodeVol(2)
    NodeSourceExt(NodeID(3)) = NodeSourceExtEqui(1,1,1,0) * NodeVol(3)
    NodeSourceExt(NodeID(4)) = NodeSourceExtEqui(1,0,1,0) * NodeVol(4)
    NodeSourceExt(NodeID(5)) = NodeSourceExtEqui(1,0,0,1) * NodeVol(5)
    NodeSourceExt(NodeID(6)) = NodeSourceExtEqui(1,1,0,1) * NodeVol(6)
    NodeSourceExt(NodeID(7)) = NodeSourceExtEqui(1,1,1,1) * NodeVol(7)
    NodeSourceExt(NodeID(8)) = NodeSourceExtEqui(1,0,1,1) * NodeVol(8)
  END DO

!#if USE_MPI
!  CALL BARRIER_AND_SYNC(NodeSourceExt_Shared_Win,MPI_COMM_SHARED)
!#endif /*USE_MPI*/
END IF ! DG_SourceExtExists

END SUBROUTINE ReadNodeSourceExtFromHDF5


!===================================================================================================================================
!> Read particle emission variables from state.h5
!> E.g. arrays containing information that have to be restored after restart (not necessarily required for automatic load balance
!> restarts, but maybe required for some)
!> Synchronize the read-in variables across all procs within the emission communicator (for the specific Species and Init) if
!> required
!===================================================================================================================================
SUBROUTINE ReadEmissionVariablesFromHDF5()
! MODULES
#if USE_MPI
USE mpi
#endif /*USE_MPI*/
!USE MOD_io_HDF5
USE MOD_Globals
!USE MOD_PreProc
USE MOD_Particle_Vars     ,ONLY: Species,nSpecies
USE MOD_Particle_Vars     ,ONLY: NeutralizationBalanceGlobal
!USE MOD_Particle_Vars     ,ONLY: NeutralizationBalance
USE MOD_HDF5_Input        ,ONLY: ReadArray,DatasetExists
USE MOD_Restart_Vars      ,ONLY: RestartFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec,iInit ! ,InitGroup
LOGICAL           :: DataExists
CHARACTER(LEN=50) :: InitName
INTEGER(KIND=IK)  :: NeutralizationBalanceTmp(1:1) ! This is a dummy array of size 1 !
!===================================================================================================================================
! Loop over all species and inits
DO iSpec=1,nSpecies
  DO iInit = 1, Species(iSpec)%NumberOfInits
    SELECT CASE(Species(iSpec)%Init(iInit)%ParticleEmissionType)
     CASE(9) ! '2D_landmark_neutralization'
       ! Re-load the value because the emission communicator can change during load balance restarts: MPIRoot is always part of this
       ! specific communicator

       ! Only the root reads the data and replaces his value, which will be communicated via the all-reduce (he also does the
       ! initial output)
       IF(MPIRoot)THEN

         IF(.NOT.FILEEXISTS(RestartFile)) &
             CALL abort(__STAMP__,'Error in ReadEmissionVariablesFromHDF5() becuase RestartFile does not exist: '//TRIM(RestartFile))

         CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
         WRITE(InitName,'(A,I0,A,I0)') 'Spec',iSpec,'Init',iInit
         CALL DatasetExists(File_ID,TRIM(InitName),DataExists)
         IF(DataExists)THEN
           CALL ReadArray(TRIM(InitName),1,(/1_IK/),0_IK,1,IntegerArray=NeutralizationBalanceTmp)
         ELSE
           !CALL abort(__STAMP__,"Read array ["//TRIM(InitName)//"] from restart file ["//TRIM(RestartFile)//"] failed.")
           WRITE (*,*) "Read array ["//TRIM(InitName)//"] from restart file ["//TRIM(RestartFile)//"] failed. "//&
             "Setting NeutralizationBalanceGlobal =0"
           NeutralizationBalanceTmp = 0
         END IF ! DataExists
         CALL CloseDataFile()

         NeutralizationBalanceGlobal = INT(NeutralizationBalanceTmp(1),4)
         !NeutralizationBalance       = NeutralizationBalanceGlobal
       END IF

! Only broadcast the information when MPI is used
!#if USE_MPI
!       ! Communicate number of particles with all procs in the same init group to the global root for output
!       InitGroup=Species(iSpec)%Init(iInit)%InitCOMM
!       ! Only processors which are part of group take part in the communication
!        CALL MPI_BCAST(Box_X, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_PICLAS, iError)
!#endif /*USE_MPI*/

     END SELECT
  END DO  ! iInit
END DO  ! iSpec=1,nSpecies

END SUBROUTINE ReadEmissionVariablesFromHDF5

#endif /*defined(PARTICLES)*/
END MODULE MOD_HDF5_Input_Particles
