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
USE MOD_MPI_vars               ,ONLY: offsetElemMPI
USE MOD_MPI_Shared_Vars        ,ONLY: ComputeNodeRootRank
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,NodeToElemInfo,NodeToElemMapping
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank,MPI_COMM_LEADERS_SHARED
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED,myLeaderGroupRank,nLeaderGroupProcs
USE MOD_MPI_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_Restart_Vars           ,ONLY: DoRestart
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
INTEGER                   :: ALLOCSTAT, iElem, i, j, k, kk, ll, mm
REAL                      :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N), DetJac(1,0:1,0:1,0:1)
REAL, ALLOCATABLE         :: Vdm_tmp(:,:)
CHARACTER(255)            :: TimeAverageFile
#if USE_MPI
INTEGER                   :: UniqueNodeID
INTEGER                   :: jElem, NonUniqueNodeID,iNode
INTEGER                   :: SendNodeCount, GlobalElemNode, GlobalElemRank, iProc
INTEGER                   :: TestElemID
LOGICAL,ALLOCATABLE       :: NodeDepoMapping(:,:)
INTEGER                   :: firstNode,lastNode
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
#if USE_MPI
  CALL Allocate_Shared((/4,nUniqueGlobalNodes/),NodeSource_Shared_Win,NodeSource_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeSource_Shared_Win,IERROR)
  NodeSource => NodeSource_Shared
  ALLOCATE(NodeSourceLoc(1:4,1:nUniqueGlobalNodes))

  IF(DoDielectricSurfaceCharge)THEN
    firstNode = INT(REAL( myComputeNodeRank   *nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))+1
    lastNode  = INT(REAL((myComputeNodeRank+1)*nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))

   ! Global, synchronized surface charge contribution (is added to NodeSource AFTER MPI synchronization)
    CALL Allocate_Shared((/nUniqueGlobalNodes/),NodeSourceExt_Shared_Win,NodeSourceExt_Shared)
    CALL MPI_WIN_LOCK_ALL(0,NodeSourceExt_Shared_Win,IERROR)
    NodeSourceExt => NodeSourceExt_Shared
    !ALLOCATE(NodeSourceExtLoc(1:1,1:nUniqueGlobalNodes))
    IF(.NOT.DoRestart)THEN
      DO iNode=firstNode, lastNode
        NodeSourceExt(iNode) = 0.
      END DO
      CALL BARRIER_AND_SYNC(NodeSourceExt_Shared_Win,MPI_COMM_SHARED)
    END IF ! .NOT.DoRestart

   ! Local, non-synchronized surface charge contribution (is added to NodeSource BEFORE MPI synchronization)
    CALL Allocate_Shared((/nUniqueGlobalNodes/),NodeSourceExtTmp_Shared_Win,NodeSourceExtTmp_Shared)
    CALL MPI_WIN_LOCK_ALL(0,NodeSourceExtTmp_Shared_Win,IERROR)
    NodeSourceExtTmp => NodeSourceExtTmp_Shared
    ALLOCATE(NodeSourceExtTmpLoc(1:nUniqueGlobalNodes))
    NodeSourceExtTmpLoc = 0.

    ! this array does not have to be initialized with zero
    ! DO iNode=firstNode, lastNode
    !   NodeSourceExtTmp(iNode) = 0.
    ! END DO
    !CALL BARRIER_AND_SYNC(NodeSourceExtTmp_Shared_Win,MPI_COMM_SHARED)
  END IF ! DoDielectricSurfaceCharge



  IF ((myComputeNodeRank.EQ.0).AND.(nLeaderGroupProcs.GT.1)) THEN
    ALLOCATE(NodeMapping(0:nLeaderGroupProcs-1))
    ALLOCATE(NodeDepoMapping(0:nLeaderGroupProcs-1, 1:nUniqueGlobalNodes))
    NodeDepoMapping = .FALSE.

    DO iElem = 1, nComputeNodeTotalElems
      ! Loop all local nodes
      DO iNode = 1, 8
        NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
        UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)

        ! Loop 1D array [offset + 1 : offset + NbrOfElems]
        ! (all CN elements that are connected to the local nodes)
        DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
          TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
          ! find the compute node
          GlobalElemNode = INT(GlobalElemRank/nComputeNodeProcessors)
          ! check if element for this side is on the current compute-node. Alternative version to the check above
          IF (GlobalElemNode.NE.myLeaderGroupRank) NodeDepoMapping(GlobalElemNode, UniqueNodeID)  = .TRUE.
        END DO
      END DO
    END DO

    DO iProc = 0, nLeaderGroupProcs - 1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      NodeMapping(iProc)%nRecvUniqueNodes = 0
      NodeMapping(iProc)%nSendUniqueNodes = 0
      CALL MPI_IRECV( NodeMapping(iProc)%nRecvUniqueNodes                       &
                  , 1                                                           &
                  , MPI_INTEGER                                                 &
                  , iProc                                                       &
                  , 666                                                         &
                  , MPI_COMM_LEADERS_SHARED                                     &
                  , RecvRequestCN(iProc)                                          &
                  , IERROR)
      DO iNode = 1, nUniqueGlobalNodes
        IF (NodeDepoMapping(iProc,iNode)) NodeMapping(iProc)%nSendUniqueNodes = NodeMapping(iProc)%nSendUniqueNodes + 1
      END DO
      CALL MPI_ISEND( NodeMapping(iProc)%nSendUniqueNodes                         &
                    , 1                                                           &
                    , MPI_INTEGER                                                 &
                    , iProc                                                       &
                    , 666                                                         &
                    , MPI_COMM_LEADERS_SHARED                                     &
                    , SendRequestCN(iProc)                                          &
                    , IERROR)
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      CALL MPI_WAIT(SendRequestCN(iProc),MPISTATUS,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      CALL MPI_WAIT(RecvRequestCN(iProc),MPISTATUS,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        ALLOCATE(NodeMapping(iProc)%RecvNodeUniqueGlobalID(1:NodeMapping(iProc)%nRecvUniqueNodes))
        ALLOCATE(NodeMapping(iProc)%RecvNodeSourceCharge(1:NodeMapping(iProc)%nRecvUniqueNodes))
        ALLOCATE(NodeMapping(iProc)%RecvNodeSourceCurrent(1:3,1:NodeMapping(iProc)%nRecvUniqueNodes))
        IF(DoDielectricSurfaceCharge) ALLOCATE(NodeMapping(iProc)%RecvNodeSourceExt(1:NodeMapping(iProc)%nRecvUniqueNodes))
        CALL MPI_IRECV( NodeMapping(iProc)%RecvNodeUniqueGlobalID                   &
                      , NodeMapping(iProc)%nRecvUniqueNodes                         &
                      , MPI_INTEGER                                                 &
                      , iProc                                                       &
                      , 666                                                         &
                      , MPI_COMM_LEADERS_SHARED                                     &
                      , RecvRequestCN(iProc)                                          &
                      , IERROR)
      END IF
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        ALLOCATE(NodeMapping(iProc)%SendNodeUniqueGlobalID(1:NodeMapping(iProc)%nSendUniqueNodes))
        NodeMapping(iProc)%SendNodeUniqueGlobalID=-1
        ALLOCATE(NodeMapping(iProc)%SendNodeSourceCharge(1:NodeMapping(iProc)%nSendUniqueNodes))
        NodeMapping(iProc)%SendNodeSourceCharge=0.
        ALLOCATE(NodeMapping(iProc)%SendNodeSourceCurrent(1:3,1:NodeMapping(iProc)%nSendUniqueNodes))
        NodeMapping(iProc)%SendNodeSourceCurrent=0.
        IF(DoDielectricSurfaceCharge) ALLOCATE(NodeMapping(iProc)%SendNodeSourceExt(1:NodeMapping(iProc)%nSendUniqueNodes))
        SendNodeCount = 0
        DO iNode = 1, nUniqueGlobalNodes
          IF (NodeDepoMapping(iProc,iNode)) THEN
            SendNodeCount = SendNodeCount + 1
            NodeMapping(iProc)%SendNodeUniqueGlobalID(SendNodeCount) = iNode
          END IF
        END DO
        CALL MPI_ISEND( NodeMapping(iProc)%SendNodeUniqueGlobalID                   &
                      , NodeMapping(iProc)%nSendUniqueNodes                         &
                      , MPI_INTEGER                                                 &
                      , iProc                                                       &
                      , 666                                                         &
                      , MPI_COMM_LEADERS_SHARED                                     &
                      , SendRequestCN(iProc)                                          &
                      , IERROR)
      END IF
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        CALL MPI_WAIT(SendRequestCN(iProc),MPISTATUS,IERROR)
        IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        CALL MPI_WAIT(RecvRequestCN(iProc),MPISTATUS,IERROR)
        IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
    END DO
  END IF
#else
  ALLOCATE(NodeSource(1:4,1:nUniqueGlobalNodes))
  IF(DoDielectricSurfaceCharge)THEN
    ALLOCATE(NodeSourceExt(1:nUniqueGlobalNodes))
    ALLOCATE(NodeSourceExtTmp(1:nUniqueGlobalNodes))
    NodeSourceExt    = 0.
    NodeSourceExtTmp = 0.
  END IF ! DoDielectricSurfaceCharge
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


SUBROUTINE FinalizeDeposition()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic deposition
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Dielectric_Vars    ,ONLY: DoDielectricSurfaceCharge
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,PeriodicSFCaseMatrix
USE MOD_PICDepo_Vars
#if USE_MPI
USE MOD_MPI_Shared_vars    ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
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
SDEALLOCATE(PartSource)

#if USE_MPI
SDEALLOCATE(NodeSourceLoc)
SDEALLOCATE(NodeMapping)
SDEALLOCATE(nDepoDOFPerProc)
SDEALLOCATE(nDepoOffsetProc)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(SendRequest)
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
    CALL UNLOCK_AND_FREE(NodeSource_Shared_Win)
    CALL UNLOCK_AND_FREE(NodeVolume_Shared_Win)

    ! Surface charging arrays
    IF(DoDielectricSurfaceCharge)THEN
      CALL UNLOCK_AND_FREE(NodeSourceExt_Shared_Win)
      CALL UNLOCK_AND_FREE(NodeSourceExtTmp_Shared_Win)
    END IF
  CASE('shape_function_adaptive')
    CALL UNLOCK_AND_FREE(SFElemr2_Shared_Win)
  END SELECT

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

  ADEALLOCATE(NodeSource_Shared)
  ADEALLOCATE(NodeVolume_Shared)
  ADEALLOCATE(NodeSourceExt_Shared)
  ADEALLOCATE(NodeSourceExtTmp_Shared)
  SDEALLOCATE(NodeSourceExtTmpLoc)
END IF ! DoDeposition

! Then, free the pointers or arrays
#endif /*USE_MPI*/

! Deposition-dependent pointers/arrays
SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean')
    ADEALLOCATE(NodeSource)
    ! Surface charging pointers/arrays
    IF(DoDielectricSurfaceCharge)THEN
      ADEALLOCATE(NodeSourceExt)
    END IF ! DoDielectricSurfaceCharge
#if USE_MPI
    ADEALLOCATE(NodeSource_Shared)
#endif /*USE_MPI*/
  CASE('shape_function_adaptive')
    ADEALLOCATE(SFElemr2_Shared)
END SELECT

END SUBROUTINE FinalizeDeposition

END MODULE MOD_PICDepo
