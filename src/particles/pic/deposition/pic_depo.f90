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

MODULE MOD_PICDepo
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
 IMPLICIT NONE
 PRIVATE
!===================================================================================================================================
INTERFACE Deposition
  MODULE PROCEDURE Deposition
END INTERFACE

INTERFACE FinalizeDeposition
  MODULE PROCEDURE FinalizeDeposition
END INTERFACE

PUBLIC:: Deposition, InitializeDeposition, FinalizeDeposition, DefineParametersPICDeposition
#if USE_MPI
PUBLIC :: ExchangeNodeSourceExtTmp
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for PIC Deposition
!==================================================================================================================================
SUBROUTINE DefineParametersPICDeposition()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools      ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%CreateStringOption( 'PIC-TimeAverageFile'      , 'Read charge density from .h5 file and save to PartSource\n'//&
                                                           'WARNING: Currently not correctly implemented for shared memory', 'none')
CALL prms%CreateLogicalOption('PIC-RelaxDeposition'      , 'Relaxation of current PartSource with RelaxFac\n'//&
                                                           'into PartSourceOld', '.FALSE.')
CALL prms%CreateRealOption(   'PIC-RelaxFac'             , 'Relaxation factor of current PartSource with RelaxFac\n'//&
                                                           'into PartSourceOld', '0.001')

CALL prms%CreateRealOption(   'PIC-shapefunction-radius'             , 'Radius of shape function'   , '1.')
CALL prms%CreateIntOption(    'PIC-shapefunction-alpha'              , 'Exponent of shape function' , '2')
CALL prms%CreateIntOption(    'PIC-shapefunction-dimension'          , '1D, 2D or 3D shape function', '3')
CALL prms%CreateIntOption(    'PIC-shapefunction-direction'          , &
    'Only required for PIC-shapefunction-dimension 1 or 2: Shape function direction for 1D (the direction in which the charge '//&
    'will be distributed) and 2D (the direction in which the charge will be constant)', '1')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-3D-deposition' ,'Deposit the charge over volume (3D)\n'//&
                                                                   ' or over a line (1D)/area(2D)\n'//&
                                                                   '1D shape function: volume or line\n'//&
                                                                   '2D shape function: volume or area', '.TRUE.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'Minimum shape function radius (for cylindrical and spherical)', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'Scaling factor of shape function radius '//&
                                                             '(for cylindrical and spherical)', '0.')
CALL prms%CreateRealOption(     'PIC-shapefunction-adaptive-DOF'  ,'Average number of DOF in shape function radius (assuming a '//&
    'Cartesian grid with equal elements). Only implemented for PIC-Deposition-Type = shape_function_adaptive (2). The maximum '//&
    'number of DOF is limited by the polynomial degree and depends on PIC-shapefunction-dimension=1, 2 or 3\n'//&
   '1D: 2*(N+1)\n'//&
   '2D: Pi*(N+1)^2\n'//&
   '3D: (4/3)*Pi*(N+1)^3\n')
CALL prms%CreateLogicalOption('PIC-shapefunction-adaptive-smoothing', 'Enable smooth transition of element-dependent radius when'//&
                                                                      ' using shape_function_adaptive.', '.FALSE.')

END SUBROUTINE DefineParametersPICDeposition


!===================================================================================================================================
!> Initialize the deposition variables first
!===================================================================================================================================
SUBROUTINE InitializeDeposition()
! MODULES
USE MOD_Globals
USE MOD_Basis                  ,ONLY: BarycentricWeights,InitializeVandermonde
USE MOD_Basis                  ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars     ,ONLY: xGP,wBary,NodeType,NodeTypeVISU
USE MOD_Mesh_Vars              ,ONLY: nElems,sJ,Vdm_EQ_N
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueGlobalNodes
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Tools          ,ONLY: CalcCellLocNodeVolumes,ReadTimeAverage
USE MOD_PICInterpolation_Vars  ,ONLY: InterpolationType
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_MPI_vars               ,ONLY: offsetElemMPI
USE MOD_MPI_Shared_Vars        ,ONLY: ComputeNodeRootRank
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,NodeToElemInfo,NodeToElemMapping
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID, GetCNElemID
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: nLeaderGroupProcs, nProcessors_Global
USE MOD_MPI_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE          :: xGP_tmp(:),wGP_tmp(:)
INTEGER                   :: ALLOCSTAT, iElem, i, j, k, kk, ll, mm, iNode
REAL                      :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N), DetJac(1,0:1,0:1,0:1)
REAL, ALLOCATABLE         :: Vdm_tmp(:,:)
CHARACTER(255)            :: TimeAverageFile
#if USE_MPI
INTEGER                   :: UniqueNodeID
INTEGER                   :: jElem, NonUniqueNodeID
INTEGER                   :: SendNodeCount, GlobalElemRank, iProc
INTEGER                   :: TestElemID, GlobalElemRankOrig, iRank
LOGICAL,ALLOCATABLE       :: NodeDepoMapping(:,:), DoNodeMapping(:), SendNode(:), IsDepoNode(:)
LOGICAL                   :: bordersMyrank
! Non-symmetric particle exchange
INTEGER                   :: SendRequestNonSymDepo(0:nProcessors_Global-1)      , RecvRequestNonSymDepo(0:nProcessors_Global-1)
INTEGER                   :: nSendUniqueNodesNonSymDepo(0:nProcessors_Global-1) , nRecvUniqueNodesNonSymDepo(0:nProcessors_Global-1)
INTEGER                   :: GlobalRankToNodeSendDepoRank(0:nProcessors_Global-1)
#endif
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE DEPOSITION...'

IF(.NOT.DoDeposition) THEN
  ! fill deposition type with empty string
  DepositionType='NONE'
  OutputSource=.FALSE.
  RelaxDeposition=.FALSE.
  RETURN
END IF

!--- Allocate arrays for charge density collection and initialize
ALLOCATE(PartSource(1:4,0:PP_N,0:PP_N,0:PP_N,nElems))
PartSource = 0.0
PartSourceConstExists=.FALSE.

!--- check if relaxation of current PartSource with RelaxFac into PartSourceOld
RelaxDeposition = GETLOGICAL('PIC-RelaxDeposition','F')
IF (RelaxDeposition) THEN
  RelaxFac     = GETREAL('PIC-RelaxFac','0.001')
#if ((USE_HDG) && (PP_nVar==1))
  ALLOCATE(PartSourceOld(1,1:2,0:PP_N,0:PP_N,0:PP_N,nElems),STAT=ALLOCSTAT)
#else
  ALLOCATE(PartSourceOld(1:4,1:2,0:PP_N,0:PP_N,0:PP_N,nElems),STAT=ALLOCSTAT)
#endif
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__,'ERROR in pic_depo.f90: Cannot allocate PartSourceOld!')
  END IF
  PartSourceOld=0.
  OutputSource = .TRUE.
ELSE
  OutputSource = GETLOGICAL('PIC-OutputSource','F')
END IF

!--- check if charge density is computed from TimeAverageFile
TimeAverageFile = GETSTR('PIC-TimeAverageFile','none')
IF (TRIM(TimeAverageFile).NE.'none') THEN
  CALL abort(__STAMP__,'This feature is currently not working! PartSource must be correctly handled in shared memory context.')
  CALL ReadTimeAverage(TimeAverageFile)
  IF (.NOT.RelaxDeposition) THEN
  !-- switch off deposition: use only the read PartSource
    DoDeposition=.FALSE.
    DepositionType='constant'
    RETURN
  ELSE
  !-- use read PartSource as initialValue for relaxation
  !-- CAUTION: will be overwritten by DG_Source if present in restart-file!
    DO iElem = 1, nElems
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
#if ((USE_HDG) && (PP_nVar==1))
            PartSourceOld(1,1,mm,ll,kk,iElem) = PartSource(4,mm,ll,kk,iElem)
            PartSourceOld(1,2,mm,ll,kk,iElem) = PartSource(4,mm,ll,kk,iElem)
#else
            PartSourceOld(1:4,1,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem)
            PartSourceOld(1:4,2,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem)
#endif
          END DO !mm
        END DO !ll
      END DO !kk
    END DO !iElem
  END IF
END IF

!--- init DepositionType-specific vars
SELECT CASE(TRIM(DepositionType))
CASE('cell_volweight')
  ALLOCATE(CellVolWeightFac(0:PP_N),wGP_tmp(0:PP_N) , xGP_tmp(0:PP_N))
  ALLOCATE(CellVolWeight_Volumes(0:1,0:1,0:1,nElems))
  CellVolWeightFac(0:PP_N) = xGP(0:PP_N)
  CellVolWeightFac(0:PP_N) = (CellVolWeightFac(0:PP_N)+1.0)/2.0
  CALL LegendreGaussNodesAndWeights(1,xGP_tmp,wGP_tmp)
  ALLOCATE( Vdm_tmp(0:1,0:PP_N))
  CALL InitializeVandermonde(PP_N,1,wBary,xGP,xGP_tmp,Vdm_tmp)
  DO iElem=1, nElems
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          DetLocal(1,i,j,k)=1./sJ(i,j,k,iElem)
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
    CALL ChangeBasis3D(1,PP_N, 1,Vdm_tmp, DetLocal(:,:,:,:),DetJac(:,:,:,:))
    DO k=0,1
      DO j=0,1
        DO i=0,1
          CellVolWeight_Volumes(i,j,k,iElem) = DetJac(1,i,j,k)*wGP_tmp(i)*wGP_tmp(j)*wGP_tmp(k)
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
  END DO
  DEALLOCATE(Vdm_tmp)
  DEALLOCATE(wGP_tmp, xGP_tmp)
CASE('cell_volweight_mean')
#if USE_MPI
  ALLOCATE(RecvRequestCN(0:nLeaderGroupProcs-1), SendRequestCN(0:nLeaderGroupProcs-1))
#endif
  IF ((TRIM(InterpolationType).NE.'cell_volweight')) THEN
    ALLOCATE(CellVolWeightFac(0:PP_N))
    CellVolWeightFac(0:PP_N) = xGP(0:PP_N)
    CellVolWeightFac(0:PP_N) = (CellVolWeightFac(0:PP_N)+1.0)/2.0
  END IF

  ! Initialize sub-cell volumes around nodes
  CALL CalcCellLocNodeVolumes()
  ALLOCATE(NodeSource(1:4,1:nUniqueGlobalNodes))
  NodeSource=0.0
  IF(DoDielectricSurfaceCharge)THEN
    ALLOCATE(NodeSourceExt(1:nUniqueGlobalNodes))
    NodeSourceExt    = 0.
  END IF ! DoDielectricSurfaceCharge
#if USE_MPI
  IF(DoDielectricSurfaceCharge)THEN
    ALLOCATE(NodeSourceExtTmp(1:nUniqueGlobalNodes))
    NodeSourceExtTmp = 0.
  END IF ! DoDielectricSurfaceCharge

  ALLOCATE(DoNodeMapping(0:nProcessors_Global-1),SendNode(1:nUniqueGlobalNodes))
  DoNodeMapping = .FALSE.
  SendNode = .FALSE.
  DO iElem = 1,nComputeNodeTotalElems
    IF (FlagShapeElem(iElem)) THEN
      bordersMyrank = .FALSE.
      ! Loop all local nodes
      TestElemID = GetGlobalElemID(iElem)
      GlobalElemRankOrig = ElemInfo_Shared(ELEM_RANK,TestElemID)
      IF (DoHaloDepo.AND.(GlobalElemRankOrig.NE.myRank)) DoNodeMapping(GlobalElemRankOrig) = .TRUE.

      DO iNode = 1, 8
        NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
        UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
        ! Loop 1D array [offset + 1 : offset + NbrOfElems]
        ! (all CN elements that are connected to the local nodes)
        DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
          TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
          ! check if element for this side is on the current compute-node. Alternative version to the check above
          IF (DoHaloDepo) THEN
            SendNode(UniqueNodeID) = .TRUE.
            IF (GlobalElemRank.NE.myRank) DoNodeMapping(GlobalElemRank) = .TRUE.
          ELSE
            IF (GlobalElemRank.EQ.myRank) THEN
              bordersMyrank = .TRUE.
              SendNode(UniqueNodeID) = .TRUE.
            END IF
          END IF
        END DO
        IF (.NOT.DoHaloDepo.AND.bordersMyrank) DoNodeMapping(GlobalElemRankOrig) = .TRUE.
      END DO
    END IF
  END DO

  nDepoNodes = 0
  ALLOCATE(IsDepoNode(1:nUniqueGlobalNodes))
  IsDepoNode = .FALSE.
  DO iElem =1, nElems
    TestElemID = GetCNElemID(iElem + offsetElem)
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,TestElemID)
      UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
      IsDepoNode(UniqueNodeID) = .TRUE.
    END DO
  END DO
  nDepoNodes = COUNT(IsDepoNode)
  nDepoNodesTotal = nDepoNodes
  DO iNode=1, nUniqueGlobalNodes
    IF (.NOT.IsDepoNode(iNode).AND.SendNode(iNode)) THEN
      nDepoNodesTotal = nDepoNodesTotal + 1
    END IF
  END DO

  ALLOCATE(DepoNodetoGlobalNode(1:nDepoNodesTotal))
  nDepoNodesTotal = 0
  DO iNode=1, nUniqueGlobalNodes
    IF (IsDepoNode(iNode)) THEN
      nDepoNodesTotal = nDepoNodesTotal + 1
      DepoNodetoGlobalNode(nDepoNodesTotal) = iNode
    END IF
  END DO
  DO iNode=1, nUniqueGlobalNodes
    IF (.NOT.IsDepoNode(iNode).AND.SendNode(iNode)) THEN
      nDepoNodesTotal = nDepoNodesTotal + 1
      DepoNodetoGlobalNode(nDepoNodesTotal) = iNode
    END IF
  END DO

  GlobalRankToNodeSendDepoRank = -1
  nNodeSendExchangeProcs = COUNT(DoNodeMapping)
  ALLOCATE(NodeSendDepoRankToGlobalRank(1:nNodeSendExchangeProcs))
  NodeSendDepoRankToGlobalRank = 0
  nNodeSendExchangeProcs = 0
  DO iRank= 0, nProcessors_Global-1
    IF (iRank.EQ.myRank) CYCLE
    IF (DoNodeMapping(iRank)) THEN
      nNodeSendExchangeProcs = nNodeSendExchangeProcs + 1
      GlobalRankToNodeSendDepoRank(iRank) = nNodeSendExchangeProcs
      NodeSendDepoRankToGlobalRank(nNodeSendExchangeProcs) = iRank
    END IF
  END DO
  ALLOCATE(NodeDepoMapping(1:nNodeSendExchangeProcs, 1:nUniqueGlobalNodes))
  NodeDepoMapping = .FALSE.

  DO iNode = 1, nUniqueGlobalNodes
    IF (SendNode(iNode)) THEN
      DO jElem = NodeToElemMapping(1,iNode) + 1, NodeToElemMapping(1,iNode) + NodeToElemMapping(2,iNode)
          TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
          ! check if element for this side is on the current compute-node. Alternative version to the check above
          IF (GlobalElemRank.NE.myRank) THEN
            iRank = GlobalRankToNodeSendDepoRank(GlobalElemRank)
            IF (iRank.LT.1) CALL ABORT(__STAMP__,'Found not connected Rank!', IERROR)
            NodeDepoMapping(iRank, iNode) = .TRUE.
          END IF
        END DO
    END IF
  END DO

  ! Get number of send nodes for each proc: Size of each message for each proc for deposition
  nSendUniqueNodesNonSymDepo         = 0
  nRecvUniqueNodesNonSymDepo(myrank) = 0
  ALLOCATE(NodeMappingSend(1:nNodeSendExchangeProcs))
  DO iProc = 1, nNodeSendExchangeProcs
    NodeMappingSend(iProc)%nSendUniqueNodes = 0
    DO iNode = 1, nUniqueGlobalNodes
      IF (NodeDepoMapping(iProc,iNode)) NodeMappingSend(iProc)%nSendUniqueNodes = NodeMappingSend(iProc)%nSendUniqueNodes + 1
    END DO
    ! local to global array
    nSendUniqueNodesNonSymDepo(NodeSendDepoRankToGlobalRank(iProc)) = NodeMappingSend(iProc)%nSendUniqueNodes
  END DO

  ! Open receive buffer for non-symmetric exchange identification
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    CALL MPI_IRECV( nRecvUniqueNodesNonSymDepo(iProc)  &
                  , 1                            &
                  , MPI_INTEGER                  &
                  , iProc                        &
                  , 1999                         &
                  , MPI_COMM_WORLD               &
                  , RecvRequestNonSymDepo(iProc)           &
                  , IERROR)
  END DO

  ! Send each proc the number of nodes that can be reached by deposition
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    CALL MPI_ISEND( nSendUniqueNodesNonSymDepo(iProc) &
                  , 1                                 &
                  , MPI_INTEGER                       &
                  , iProc                             &
                  , 1999                              &
                  , MPI_COMM_WORLD                    &
                  , SendRequestNonSymDepo(iProc)      &
                  , IERROR)
  END DO

  ! Finish communication
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO

  nNodeRecvExchangeProcs = COUNT(nRecvUniqueNodesNonSymDepo.GT.0)
  ALLOCATE(NodeMappingRecv(1:nNodeRecvExchangeProcs))
  ALLOCATE(NodeRecvDepoRankToGlobalRank(1:nNodeRecvExchangeProcs))
  NodeRecvDepoRankToGlobalRank = 0
  nNodeRecvExchangeProcs = 0
  DO iRank= 0, nProcessors_Global-1
    IF (iRank.EQ.myRank) CYCLE
    IF (nRecvUniqueNodesNonSymDepo(iRank).GT.0) THEN
      nNodeRecvExchangeProcs = nNodeRecvExchangeProcs + 1
      ! Store global rank of iRecvRank
      NodeRecvDepoRankToGlobalRank(nNodeRecvExchangeProcs) = iRank
      ! Store number of nodes of iRecvRank
      NodeMappingRecv(nNodeRecvExchangeProcs)%nRecvUniqueNodes = nRecvUniqueNodesNonSymDepo(iRank)
    END IF
  END DO

  ! Open receive buffer
  ALLOCATE(RecvRequest(1:nNodeRecvExchangeProcs))
  DO iProc = 1, nNodeRecvExchangeProcs
    ALLOCATE(NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
    ALLOCATE(NodeMappingRecv(iProc)%RecvNodeSourceCharge(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
    ALLOCATE(NodeMappingRecv(iProc)%RecvNodeSourceCurrent(1:3,1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
    IF(DoDielectricSurfaceCharge) ALLOCATE(NodeMappingRecv(iProc)%RecvNodeSourceExt(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
    CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID                   &
                  , NodeMappingRecv(iProc)%nRecvUniqueNodes                         &
                  , MPI_INTEGER                                                 &
                  , NodeRecvDepoRankToGlobalRank(iProc)                         &
                  , 666                                                         &
                  , MPI_COMM_WORLD                                              &
                  , RecvRequest(iProc)                                          &
                  , IERROR)
  END DO

  ! Open send buffer
  ALLOCATE(SendRequest(1:nNodeSendExchangeProcs))
  DO iProc = 1, nNodeSendExchangeProcs
    ALLOCATE(NodeMappingSend(iProc)%SendNodeUniqueGlobalID(1:NodeMappingSend(iProc)%nSendUniqueNodes))
    NodeMappingSend(iProc)%SendNodeUniqueGlobalID=-1
    ALLOCATE(NodeMappingSend(iProc)%SendNodeSourceCharge(1:NodeMappingSend(iProc)%nSendUniqueNodes))
    NodeMappingSend(iProc)%SendNodeSourceCharge=0.
    ALLOCATE(NodeMappingSend(iProc)%SendNodeSourceCurrent(1:3,1:NodeMappingSend(iProc)%nSendUniqueNodes))
    NodeMappingSend(iProc)%SendNodeSourceCurrent=0.
    IF(DoDielectricSurfaceCharge) ALLOCATE(NodeMappingSend(iProc)%SendNodeSourceExt(1:NodeMappingSend(iProc)%nSendUniqueNodes))
    SendNodeCount = 0
    DO iNode = 1, nUniqueGlobalNodes
      IF (NodeDepoMapping(iProc,iNode)) THEN
        SendNodeCount = SendNodeCount + 1
        NodeMappingSend(iProc)%SendNodeUniqueGlobalID(SendNodeCount) = iNode
      END IF
    END DO
    CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeUniqueGlobalID                   &
                  , NodeMappingSend(iProc)%nSendUniqueNodes                         &
                  , MPI_INTEGER                                                 &
                  , NodeSendDepoRankToGlobalRank(iProc)                         &
                  , 666                                                         &
                  , MPI_COMM_WORLD                                              &
                  , SendRequest(iProc)                                          &
                  , IERROR)
  END DO

  ! Finish send
  DO iProc = 1, nNodeSendExchangeProcs
    CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO

  ! Finish receive
  DO iProc = 1, nNodeRecvExchangeProcs
    CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO
#else
  nDepoNodes      = nUniqueGlobalNodes
  nDepoNodesTotal = nDepoNodes
  ALLOCATE(DepoNodetoGlobalNode(1:nDepoNodesTotal))
  DO iNode=1, nUniqueGlobalNodes
    DepoNodetoGlobalNode(iNode) = iNode
  END DO
#endif /*USE_MPI*/

  IF(DoDielectricSurfaceCharge)THEN
    ! Allocate and determine Vandermonde mapping from equidistant (visu) to NodeType node set
    ALLOCATE(Vdm_EQ_N(0:PP_N,0:1))
    CALL GetVandermonde(1, NodeTypeVISU, PP_N, NodeType, Vdm_EQ_N, modal=.FALSE.)
  END IF ! DoDielectricSurfaceCharge

CASE('shape_function', 'shape_function_cc', 'shape_function_adaptive')
#if USE_MPI
  ALLOCATE(RecvRequest(nShapeExchangeProcs),SendRequest(nShapeExchangeProcs))
#endif
  ! --- Set shape function radius in each cell when using adaptive shape function
  IF(TRIM(DepositionType).EQ.'shape_function_adaptive') CALL InitShapeFunctionAdaptive()

  ! --- Set integration factor only for uncorrected shape function methods
  SELECT CASE(TRIM(DepositionType))
  CASE('shape_function_cc', 'shape_function_adaptive')
    w_sf  = 1.0 ! set dummy value
  END SELECT

  ! --- Set periodic case matrix for shape function deposition (virtual displacement of particles in the periodic directions)
  CALL InitPeriodicSFCaseMatrix()

  ! --- Set element flag for cycling already completed elements
!  ALLOCATE(PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N,1:nElems))
#if USE_MPI
  ALLOCATE(ChargeSFDone(1:nComputeNodeTotalElems))
  ALLOCATE(nDepoDOFPerProc(0:nComputeNodeProcessors-1),nDepoOffsetProc(0:nComputeNodeProcessors-1))
  nDepoDOFPerProc =0; nDepoOffsetProc = 0
  DO iProc = 0, nComputeNodeProcessors-1
    nDepoDOFPerProc(iProc) = (offsetElemMPI(iProc + 1 + ComputeNodeRootRank) - offsetElemMPI(iProc + ComputeNodeRootRank))*(PP_N+1)**3*4
  END DO
  DO iProc = 0, nComputeNodeProcessors-1
    nDepoOffsetProc(iProc) = (offsetElemMPI(iProc + ComputeNodeRootRank) - offsetElemMPI(ComputeNodeRootRank))*(PP_N+1)**3*4
  END DO
  IF(myComputeNodeRank.EQ.0) THEN
   ALLOCATE(PartSourceGlob(1:4,0:PP_N,0:PP_N,0:PP_N,1:nComputeNodeTotalElems))
  END IF
#else
  ALLOCATE(ChargeSFDone(1:nElems))
#endif /*USE_MPI*/

CASE DEFAULT
  CALL abort(__STAMP__,'Unknown DepositionType in pic_depo.f90')
END SELECT

IF (PartSourceConstExists) THEN
  ALLOCATE(PartSourceConst(1:4,0:PP_N,0:PP_N,0:PP_N,nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in pic_depo.f90: Cannot allocate PartSourceConst!')
  PartSourceConst=0.
END IF

ALLOCATE(PartSourceTmp(    1:4,0:PP_N,0:PP_N,0:PP_N))

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition


!===================================================================================================================================
!> Calculate the shape function radius for each element depending on the neighbouring element sizes and the own element size
!===================================================================================================================================
SUBROUTINE InitShapeFunctionAdaptive()
! MODULES
USE MOD_Preproc
USE MOD_Globals                     ,ONLY: UNIT_stdOut,abort
USE MOD_PICDepo_Vars                ,ONLY: SFAdaptiveDOF,SFAdaptiveSmoothing,SFElemr2_Shared,dim_sf,dimFactorSF
USE MOD_ReadInTools                 ,ONLY: GETREAL,GETLOGICAL
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,BoundsOfElem_Shared
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID
USE MOD_Globals_Vars                ,ONLY: PI
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemMidPoint_Shared,ElemToElemMapping,ElemToElemInfo
USE MOD_Mesh_Tools                  ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars          ,ONLY: NodeCoords_Shared
USE MOD_PICDepo_Shapefunction_Tools ,ONLY: SFNorm
#if USE_MPI
USE MOD_Globals                     ,ONLY: IERROR
USE MOD_PICDepo_Vars                ,ONLY: SFElemr2_Shared_Win
USE MOD_Globals                     ,ONLY: MPIRoot
USE MOD_MPI_Shared_Vars             ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank,MPI_COMM_SHARED
USE MOD_MPI_Shared
#else
USE MOD_Mesh_Vars                   ,ONLY: nElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: UniqueNodeID,NonUniqueNodeID,iNode,NeighUniqueNodeID
REAL                           :: SFDepoScaling
LOGICAL                        :: ElemDone
INTEGER                        :: ppp,globElemID
REAL                           :: r_sf_tmp,SFAdaptiveDOFDefault,DOFMax
INTEGER                        :: iCNElem,firstElem,lastElem,jNode,NbElemID,NeighNonUniqueNodeID
CHARACTER(32)                  :: hilf2,hilf3
#if USE_MPI
#endif /*USE_MPI*/
REAL                           :: CharacteristicLength,BoundingBoxVolume
!===================================================================================================================================
! Set the number of DOF/SF
! Check which shape function dimension is used and set default value
SELECT CASE(dim_sf)
CASE(1)
  SFAdaptiveDOFDefault=2.0*(1.+1.)
  DOFMax = 2.0*(REAL(PP_N)+1.) ! Max. DOF per element in 1D
  hilf2 = '2*(N+1)'            ! Max. DOF per element in 1D for abort message
CASE(2)
  SFAdaptiveDOFDefault=PI*(1.+1.)**2
  DOFMax = PI*(REAL(PP_N)+1.)**2 ! Max. DOF per element in 2D
  hilf2 = 'PI*(N+1)**2'          ! Max. DOF per element in 1D for abort message
CASE(3)
  SFAdaptiveDOFDefault=(4./3.)*PI*(1.+1.)**3
  DOFMax = (4./3.)*PI*(REAL(PP_N)+1.)**3 ! Max. DOF per element in 2D
  hilf2 = '(4/3)*PI*(N+1)**3'            ! Max. DOF per element in 1D for abort message
END SELECT
WRITE(UNIT=hilf3,FMT='(G0)') SFAdaptiveDOFDefault
SFAdaptiveDOF = GETREAL('PIC-shapefunction-adaptive-DOF',TRIM(hilf3))

IF(SFAdaptiveDOF.GT.DOFMax)THEN
  SWRITE(UNIT_StdOut,'(A,F10.2)') "         PIC-shapefunction-adaptive-DOF =", SFAdaptiveDOF
  SWRITE(UNIT_StdOut,'(A,A19,A,F10.2)') " Maximum allowed is ",TRIM(hilf2)," =", DOFMax
  SWRITE(UNIT_StdOut,*) "Reduce the number of DOF/SF in order to have no DOF outside of the deposition range (neighbour elems)"
  SWRITE(UNIT_StdOut,*) "Set a value lower or equal to than the maximum for a given polynomial degree N\n"
  SWRITE(UNIT_StdOut,*) "              N:     1      2      3      4      5       6       7"
  SWRITE(UNIT_StdOut,*) "  ----------------------------------------------------------------"
  SWRITE(UNIT_StdOut,*) "           | 1D:     4      6      8     10     12      14      16"
  SWRITE(UNIT_StdOut,*) "  Max. DOF | 2D:    12     28     50     78    113     153     201"
  SWRITE(UNIT_StdOut,*) "           | 3D:    33    113    268    523    904    1436    2144"
  SWRITE(UNIT_StdOut,*) "  ----------------------------------------------------------------"
  CALL abort(__STAMP__,'PIC-shapefunction-adaptive-DOF > '//TRIM(hilf2)//' is not allowed')
ELSE
  ! Check which shape function dimension is used
  SELECT CASE(dim_sf)
  CASE(1)
    SFDepoScaling  = SFAdaptiveDOF/2.0
  CASE(2)
    SFDepoScaling  = SQRT(SFAdaptiveDOF/PI)
  CASE(3)
    SFDepoScaling  = (3.*SFAdaptiveDOF/(4.*PI))**(1./3.)
  END SELECT
END IF

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

CALL Allocate_Shared((/2,nComputeNodeTotalElems/),SFElemr2_Shared_Win,SFElemr2_Shared)
CALL MPI_WIN_LOCK_ALL(0,SFElemr2_Shared_Win,IERROR)
#else
ALLOCATE(SFElemr2_Shared(1:2,1:nElems))
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif
  SFElemr2_Shared = HUGE(1.)
#if USE_MPI
END IF
CALL BARRIER_AND_SYNC(SFElemr2_Shared_Win,MPI_COMM_SHARED)
#endif
DO iCNElem = firstElem,lastElem
  ElemDone = .FALSE.

  DO ppp = 1,ElemToElemMapping(2,iCNElem)
    ! Get neighbour global element ID
    globElemID = GetGlobalElemID(ElemToElemInfo(ElemToElemMapping(1,iCNElem)+ppp))
    NbElemID = GetCNElemID(globElemID)
    ! Loop neighbour nodes
    Nodeloop: DO jNode = 1, 8
      NeighNonUniqueNodeID = ElemNodeID_Shared(jNode,NbElemID)
      NeighUniqueNodeID = NodeInfo_Shared(NeighNonUniqueNodeID)
      ! Loop my nodes
      DO iNode = 1, 8
        NonUniqueNodeID = ElemNodeID_Shared(iNode,iCNElem)
        UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
        IF (UniqueNodeID.EQ.NeighUniqueNodeID) CYCLE Nodeloop ! Skip coinciding nodes of my and my neighbours element
      END DO
      ElemDone =.TRUE.

      ! Measure distance from my corner nodes to neighbour elem corner nodes
      DO iNode = 1, 8
        NonUniqueNodeID = ElemNodeID_Shared(iNode,iCNElem)

        ! Only measure distances in the dimension in which the nodes to not coincide (i.e. they are not projected onto each other
        ! in 1D or 2D deposition)
        ASSOCIATE( v1 => NodeCoords_Shared(1:3,NonUniqueNodeID)      ,&
                   v2 => NodeCoords_Shared(1:3,NeighNonUniqueNodeID) )
          IF(SFMeasureDistance(v1,v2)) THEN
            r_sf_tmp = SFNorm(v1-v2)
            IF (r_sf_tmp.LT.SFElemr2_Shared(1,iCNElem)) SFElemr2_Shared(1,iCNElem) = r_sf_tmp
          END IF ! SFMeasureDistance(v1,v2)
        END ASSOCIATE
      END DO ! iNode = 1, 8

    END DO Nodeloop
  END DO

  ! Sanity check if no neighbours are present
  IF (.NOT.ElemDone) THEN
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iCNElem)
      r_sf_tmp = SFNorm(ElemMidPoint_Shared(1:3,iCNElem)-NodeCoords_Shared(1:3,NonUniqueNodeID))
      IF (r_sf_tmp.LT.SFElemr2_Shared(1,iCNElem)) SFElemr2_Shared(1,iCNElem) = r_sf_tmp
    END DO
  END IF

  ! Because ElemVolume_Shared(CNElemID) is not available for halo elements, the bounding box volume is used as an approximate
  ! value for the element volume from which the characteristic length of the element is calculated
  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,GetGlobalElemID(iCNElem)) ) ! 1-2: Min, Max value; 1-3: x,y,z
    BoundingBoxVolume = (Bounds(2,1)-Bounds(1,1)) * (Bounds(2,2)-Bounds(1,2)) * (Bounds(2,3)-Bounds(1,3))
  END ASSOCIATE
  ! Check which shape function dimension is used
  SELECT CASE(dim_sf)
  CASE(1)
    !CharacteristicLength = ElemVolume_Shared(iCNElem) / dimFactorSF
    CharacteristicLength = BoundingBoxVolume / dimFactorSF
  CASE(2)
    !CharacteristicLength = SQRT(ElemVolume_Shared(iCNElem) / dimFactorSF)
    CharacteristicLength = SQRT(BoundingBoxVolume / dimFactorSF)
  CASE(3)
    !CharacteristicLength = ElemCharLength_Shared(iCNElem)
    CharacteristicLength = BoundingBoxVolume**(1./3.)
  END SELECT

  ! Check characteristic length of cell (or when using SFAdaptiveSmoothing)
  IF(CharacteristicLength.LT.SFElemr2_Shared(1,iCNElem).OR.SFAdaptiveSmoothing)THEN
    SFElemr2_Shared(1,iCNElem) = (SFElemr2_Shared(1,iCNElem) + CharacteristicLength)/2.0
  END IF

  ! Scale the radius so that it reaches at most the neighbouring cells but no further (all neighbours of the 8 corner nodes)
  SFElemr2_Shared(1,iCNElem) = SFElemr2_Shared(1,iCNElem) * SFDepoScaling / (PP_N+1.)
  SFElemr2_Shared(2,iCNElem) = SFElemr2_Shared(1,iCNElem)**2

  ! Sanity checks
  IF(SFElemr2_Shared(1,iCNElem).LE.0.0)      CALL abort(__STAMP__,'Shape function radius <= zero!')
  IF(SFElemr2_Shared(1,iCNElem).GE.HUGE(1.)) CALL abort(__STAMP__,'Shape function radius >= HUGE(1.)!')
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(SFElemr2_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE InitShapeFunctionAdaptive


!===================================================================================================================================
!> Fill PeriodicSFCaseMatrix when using shape function deposition in combination with periodic boundaries
!===================================================================================================================================
SUBROUTINE InitPeriodicSFCaseMatrix()
! MODULES
USE MOD_Globals            ,ONLY: UNIT_StdOut
#if USE_MPI
USE MOD_Globals            ,ONLY: MPIRoot
#endif /*USE_MPI*/
USE MOD_Particle_Mesh_Vars ,ONLY: PeriodicSFCaseMatrix,NbrOfPeriodicSFCases
USE MOD_PICDepo_Vars       ,ONLY: dim_sf,dim_periodic_vec1,dim_periodic_vec2,dim_sf_dir1,dim_sf_dir2
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_ReadInTools        ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: I,J
!===================================================================================================================================
IF (GEO%nPeriodicVectors.LE.0) THEN

  ! Set defaults and return in non-periodic case
  NbrOfPeriodicSFCases = 1
  ALLOCATE(PeriodicSFCaseMatrix(1:1,1:3))
  PeriodicSFCaseMatrix(:,:) = 0

ELSE

  ! Build case matrix:
  ! Particles may move in more periodic directions than their charge is deposited, e.g., fully periodic in combination with
  ! 1D shape function
  NbrOfPeriodicSFCases = 3**dim_sf

  ALLOCATE(PeriodicSFCaseMatrix(1:NbrOfPeriodicSFCases,1:3))
  PeriodicSFCaseMatrix(:,:) = 0
  IF (dim_sf.EQ.1) THEN
    PeriodicSFCaseMatrix(1,1) = 1
    PeriodicSFCaseMatrix(3,1) = -1
  END IF
  IF (dim_sf.EQ.2) THEN
    PeriodicSFCaseMatrix(1:3,1) = 1
    PeriodicSFCaseMatrix(7:9,1) = -1
    DO I = 1,3
      PeriodicSFCaseMatrix(I*3-2,2) = 1
      PeriodicSFCaseMatrix(I*3,2) = -1
    END DO
  END IF
  IF (dim_sf.EQ.3) THEN
    PeriodicSFCaseMatrix(1:9,1) = 1
    PeriodicSFCaseMatrix(19:27,1) = -1
    DO I = 1,3
      PeriodicSFCaseMatrix(I*9-8:I*9-6,2) = 1
      PeriodicSFCaseMatrix(I*9-2:I*9,2) = -1
      DO J = 1,3
        PeriodicSFCaseMatrix((J*3-2)+(I-1)*9,3) = 1
        PeriodicSFCaseMatrix((J*3)+(I-1)*9,3) = -1
      END DO
    END DO
  END IF

  ! Define which of the periodic vectors are used for 2D shape function and display info
  IF(dim_sf.EQ.2)THEN
    IF(GEO%nPeriodicVectors.EQ.1)THEN
      dim_periodic_vec1 = 1
      dim_periodic_vec2 = 0
    ELSEIF(GEO%nPeriodicVectors.EQ.2)THEN
      dim_periodic_vec1 = 1
      dim_periodic_vec2 = 2
    ELSEIF(GEO%nPeriodicVectors.EQ.3)THEN
      dim_periodic_vec1 = dim_sf_dir1
      dim_periodic_vec2 = dim_sf_dir2
    END IF ! GEO%nPeriodicVectors.EQ.1
    CALL PrintOption('Dimension of 1st periodic vector for 2D shape function','INFO',IntOpt=dim_periodic_vec1)
    SWRITE(UNIT_StdOut,*) "1st PeriodicVector =", GEO%PeriodicVectors(1:3,dim_periodic_vec1)
    CALL PrintOption('Dimension of 2nd periodic vector for 2D shape function','INFO',IntOpt=dim_periodic_vec2)
    SWRITE(UNIT_StdOut,*) "2nd PeriodicVector =", GEO%PeriodicVectors(1:3,dim_periodic_vec2)
  END IF ! dim_sf.EQ.2

END IF

END SUBROUTINE InitPeriodicSFCaseMatrix


SUBROUTINE Deposition(doParticle_In, stage_opt)
!============================================================================================================================
! This subroutine performs the deposition of the particle charge and current density to the grid
! following list of distribution methods are implemented
! - shape function       (only one type implemented)
! useVMPF added, therefore, this routine contains automatically the use of variable mpfs
!============================================================================================================================
! USE MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Analyze_Vars ,ONLY: DoVerifyCharge,PartAnalyzeStep
USE MOD_Particle_Vars
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Method        ,ONLY: DepositionMethod
USE MOD_PIC_Analyze           ,ONLY: VerifyDepositedCharge
USE MOD_TimeDisc_Vars         ,ONLY: iter
#if USE_MPI
USE MOD_MPI_Shared            ,ONLY: BARRIER_AND_SYNC
#endif  /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT variable declaration
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
INTEGER,INTENT(IN),OPTIONAL      :: stage_opt ! TODO: definition of this variable
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration
INTEGER         :: stage
!-----------------------------------------------------------------------------------------------------------------------------------
!============================================================================================================================
! Return, if no deposition is required
IF(.NOT.DoDeposition) RETURN
IF (PRESENT(stage_opt)) THEN
  stage = stage_opt
ELSE
  stage = 0
END IF

IF((stage.EQ.0).OR.(stage.EQ.1)) PartSource = 0.0

IF(PRESENT(doParticle_In)) THEN
  CALL DepositionMethod(doParticle_In, stage_opt=stage)
ELSE
  CALL DepositionMethod(stage_opt=stage)
END IF

IF((stage.EQ.0).OR.(stage.EQ.4)) THEN
  IF(MOD(iter,PartAnalyzeStep).EQ.0) THEN
    IF(DoVerifyCharge) CALL VerifyDepositedCharge()
  END IF
END IF

END SUBROUTINE Deposition


PPURE LOGICAL FUNCTION SFMeasureDistance(v1,v2)
!============================================================================================================================
! Check if the two position vectors coincide in the 1D or 2D projection. If yes, then return .FALSE., else return .TRUE.
! If two points coincide in the direction in which the shape function is not deposited, they are ignored (coincide means that the
! real values are equal up to relative precision of 1e-5)
!============================================================================================================================
USE MOD_PICDepo_Vars ,ONLY: dim_sf,dim_sf_dir,dim_sf_dir1,dim_sf_dir2
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: v1(1:3) !< Input vector 1
REAL, INTENT(IN) :: v2(1:3) !< Input vector 2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SFMeasureDistance = .TRUE. ! Default, also used for dim_sf=3 (3D case)

! Depending on the dimensionality
SELECT CASE (dim_sf)
CASE (1)
  SFMeasureDistance = MERGE(.FALSE. , .TRUE. , ALMOSTEQUALRELATIVE(v1(dim_sf_dir) , v2(dim_sf_dir) , 1e-6))
CASE (2)
  SFMeasureDistance = MERGE(.FALSE. , .TRUE. , ALMOSTEQUALRELATIVE(v1(dim_sf_dir1) , v2(dim_sf_dir1) , 1e-6) .AND. &
                                               ALMOSTEQUALRELATIVE(v1(dim_sf_dir2) , v2(dim_sf_dir2) , 1e-6)       )
END SELECT

END FUNCTION SFMeasureDistance


#if USE_MPI
!===================================================================================================================================
!> Exchange the node source container between MPI processes (either during load balance or hdf5 output) and nullify the local charge
!> container NodeSourceExtTmp. Updates the node charge container NodeSourceExt at MPI interfaces.
!===================================================================================================================================
SUBROUTINE ExchangeNodeSourceExtTmp()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: NodeMappingRecv,NodeMappingSend,NodeSourceExtTmp
USE MOD_PICDepo_Vars       ,ONLY: nDepoNodesTotal,nNodeSendExchangeProcs,NodeSendDepoRankToGlobalRank,DepoNodetoGlobalNode
USE MOD_PICDepo_Vars       ,ONLY: nNodeRecvExchangeProcs
USE MOD_PICDepo_Vars       ,ONLY: NodeRecvDepoRankToGlobalRank
#endif  /*USE_MPI*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                        :: iProc
!INTEGER                        :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER                        :: RecvRequest(1:nNodeRecvExchangeProcs),SendRequest(1:nNodeSendExchangeProcs)
!INTEGER                        :: MessageSize
#endif /*USE_MPI*/
INTEGER                        :: globalNode, iNode
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)                :: CounterStart,CounterEnd
REAL(KIND=8)                   :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! 1) Receive charge density
DO iProc = 1, nNodeRecvExchangeProcs
  ! Open receive buffer
  CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeSourceExt(:) &
      , NodeMappingRecv(iProc)%nRecvUniqueNodes            &
      , MPI_DOUBLE_PRECISION                           &
      , NodeRecvDepoRankToGlobalRank(iProc)                &
      , 666                                            &
      , MPI_COMM_WORLD                                 &
      , RecvRequest(iProc)                             &
      , IERROR)
END DO


DO iProc = 1, nNodeSendExchangeProcs
  ! Send message (non-blocking)
  DO iNode = 1, NodeMappingSend(iProc)%nSendUniqueNodes
    NodeMappingSend(iProc)%SendNodeSourceExt(iNode) = NodeSourceExtTmp(NodeMappingSend(iProc)%SendNodeUniqueGlobalID(iNode))
  END DO
  CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeSourceExt(:) &
      , NodeMappingSend(iProc)%nSendUniqueNodes        &
      , MPI_DOUBLE_PRECISION                       &
      , NodeSendDepoRankToGlobalRank(iProc)            &
      , 666                                        &
      , MPI_COMM_WORLD                             &
      , SendRequest(iProc)                         &
      , IERROR)
END DO
! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc = 1, nNodeSendExchangeProcs
  CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO
DO iProc = 1, nNodeRecvExchangeProcs
  CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6) = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
#endif /*defined(MEASURE_MPI_WAIT)*/

! 3) Extract messages
DO iProc = 1, nNodeRecvExchangeProcs
  DO iNode = 1, NodeMappingRecv(iProc)%nRecvUniqueNodes
    ASSOCIATE( NS => NodeSourceExtTmp(NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(iNode)))
      NS = NS + NodeMappingRecv(iProc)%RecvNodeSourceExt(iNode)
    END ASSOCIATE
  END DO
END DO

! Add NodeSourceExtTmp values of the last boundary interaction
DO iNode = 1, nDepoNodesTotal
  globalNode = DepoNodetoGlobalNode(iNode)
  NodeSourceExt(globalNode) = NodeSourceExt(globalNode) + NodeSourceExtTmp(globalNode)
END DO
! Reset local surface charge
NodeSourceExtTmp = 0.
END SUBROUTINE ExchangeNodeSourceExtTmp
#endif /*USE_MPI*/


SUBROUTINE FinalizeDeposition()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic deposition
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,PeriodicSFCaseMatrix
USE MOD_PICDepo_Vars
#if USE_MPI
USE MOD_MPI_Shared_vars    ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
#if USE_LOADBALANCE
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: PartSourceLB,NodeSourceExtEquiLB
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Dielectric_Vars    ,ONLY: DoDielectricSurfaceCharge
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Mesh_Vars          ,ONLY: offsetElem
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
INTEGER,PARAMETER :: N_variables=1
INTEGER           :: iElem
INTEGER           :: NodeID(1:8)
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
SDEALLOCATE(PartSourceConst)
SDEALLOCATE(PartSourceOld)
SDEALLOCATE(PartSourceTmp)
SDEALLOCATE(GaussBorder)
SDEALLOCATE(Vdm_EquiN_GaussN)
SDEALLOCATE(Knots)
SDEALLOCATE(GaussBGMIndex)
SDEALLOCATE(GaussBGMFactor)
SDEALLOCATE(GEO%PeriodicBGMVectors)
SDEALLOCATE(BGMSource)
SDEALLOCATE(GPWeight)
SDEALLOCATE(ElemRadius2_sf)
SDEALLOCATE(Vdm_NDepo_GaussN)
SDEALLOCATE(DDMassInv)
SDEALLOCATE(XiNDepo)
SDEALLOCATE(swGPNDepo)
SDEALLOCATE(wBaryNDepo)
SDEALLOCATE(NDepochooseK)
SDEALLOCATE(tempcharge)
SDEALLOCATE(CellVolWeightFac)
SDEALLOCATE(CellVolWeight_Volumes)
SDEALLOCATE(ChargeSFDone)
!SDEALLOCATE(PartSourceLoc)
SDEALLOCATE(PartSourceGlob)
SDEALLOCATE(PeriodicSFCaseMatrix)

#if USE_MPI
SDEALLOCATE(FlagShapeElem)
SDEALLOCATE(nDepoDOFPerProc)
SDEALLOCATE(nDepoOffsetProc)
SDEALLOCATE(RecvRequestCN)
SDEALLOCATE(SendRequestCN)
SDEALLOCATE(SendElemShapeID)
SDEALLOCATE(CNRankToSendRank)

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

IF(DoDeposition)THEN
  ! Deposition-dependent arrays
  SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean')
    CALL UNLOCK_AND_FREE(NodeVolume_Shared_Win)
  CASE('shape_function_adaptive')
    CALL UNLOCK_AND_FREE(SFElemr2_Shared_Win)
  END SELECT

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

  ADEALLOCATE(NodeVolume_Shared)
END IF ! DoDeposition

! Then, free the pointers or arrays
#endif /*USE_MPI*/

! Deposition-dependent pointers/arrays
SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean')
  CASE('shape_function_adaptive')
    ADEALLOCATE(SFElemr2_Shared)
END SELECT

#if USE_LOADBALANCE
IF ((PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) .AND. DoDeposition) THEN
  SDEALLOCATE(PartSourceLB)
  ALLOCATE(PartSourceLB(1:4,0:PP_N,0:PP_N,0:PP_N,nElems))
  PartSourceLB = PartSource
  IF(DoDielectricSurfaceCharge)THEN
    CALL ExchangeNodeSourceExtTmp()
    !SDEALLOCATE(NodeSourceExtEquiLB)
    !ALLOCATE(NodeSourceExtEquiLB(1:4,0:PP_N,0:PP_N,0:PP_N,nElems))
    ALLOCATE(NodeSourceExtEquiLB(1:N_variables,0:1,0:1,0:1,nElems))
    ! Loop over all elements and store absolute charge values in equidistantly distributed nodes of PP_N=1
    DO iElem=1,PP_nElems
      ! Copy values to equidistant distribution
      NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem+offsetElem)))
      NodeSourceExtEquiLB(1,0,0,0,iElem) = NodeSourceExt(NodeID(1))
      NodeSourceExtEquiLB(1,1,0,0,iElem) = NodeSourceExt(NodeID(2))
      NodeSourceExtEquiLB(1,1,1,0,iElem) = NodeSourceExt(NodeID(3))
      NodeSourceExtEquiLB(1,0,1,0,iElem) = NodeSourceExt(NodeID(4))
      NodeSourceExtEquiLB(1,0,0,1,iElem) = NodeSourceExt(NodeID(5))
      NodeSourceExtEquiLB(1,1,0,1,iElem) = NodeSourceExt(NodeID(6))
      NodeSourceExtEquiLB(1,1,1,1,iElem) = NodeSourceExt(NodeID(7))
      NodeSourceExtEquiLB(1,0,1,1,iElem) = NodeSourceExt(NodeID(8))
    END DO!iElem
    ADEALLOCATE(ElemNodeID_Shared)
  END IF ! DoDielectricSurfaceCharge
END IF
#endif /*USE_LOADBALANCE*/

SDEALLOCATE(PartSource)
SDEALLOCATE(DepoNodetoGlobalNode)
SDEALLOCATE(NodeSource)
SDEALLOCATE(NodeSourceExt)

#if USE_MPI
SDEALLOCATE(NodeSendDepoRankToGlobalRank)
SDEALLOCATE(NodeRecvDepoRankToGlobalRank)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(SendRequest)
SDEALLOCATE(NodeSourceExtTmp)
SDEALLOCATE(NodeMappingSend)
SDEALLOCATE(NodeMappingRecv)
#endif /*USE_MPI*/

END SUBROUTINE FinalizeDeposition

END MODULE MOD_PICDepo
