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

INTERFACE InitializeDeposition
  MODULE PROCEDURE InitializeDeposition
END INTERFACE

INTERFACE FinalizeDeposition
  MODULE PROCEDURE FinalizeDeposition
END INTERFACE

PUBLIC:: Deposition, InitializeDeposition, FinalizeDeposition
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeDeposition
!===================================================================================================================================
! Initialize the deposition variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Tools          ,ONLY: CalcCellLocNodeVolumes,ReadTimeAverage,beta
USE MOD_Particle_Vars
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Mesh_Vars              ,ONLY: nElems,Elem_xGP,sJ,nGlobalElems,Vdm_EQ_N
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,MeshVolume
USE MOD_Interpolation_Vars     ,ONLY: xGP,wBary,NodeType,NodeTypeVISU
USE MOD_Basis                  ,ONLY: BarycentricWeights,InitializeVandermonde
USE MOD_Basis                  ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY
USE MOD_PICInterpolation_Vars  ,ONLY: InterpolationType
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueGlobalNodes
#if USE_MPI
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems, nComputeNodeProcessors, myComputeNodeRank, MPI_COMM_LEADERS_SHARED
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED, myLeaderGroupRank, nLeaderGroupProcs
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemNodeID_Shared, NodeInfo_Shared, ElemInfo_Shared, NodeToElemInfo, NodeToElemMapping
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
#endif
USE MOD_ReadInTools            ,ONLY: PrintOption
#if USE_MPI
USE MOD_PICDepo_MPI            ,ONLY: MPIBackgroundMeshInit
#endif /*USE_MPI*/
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_PICDepo_Method         ,ONLY: InitDepositionMethod
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE          :: wBary_tmp(:),Vdm_GaussN_EquiN(:,:)
REAL,ALLOCATABLE          :: xGP_tmp(:),wGP_tmp(:)
INTEGER                   :: ALLOCSTAT, iElem, i, j, k, iBC, kk, ll, mm
REAL                      :: VolumeShapeFunction,r_sf_average
REAL                      :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N), DetJac(1,0:1,0:1,0:1)
REAL, ALLOCATABLE         :: Vdm_tmp(:,:)
CHARACTER(32)             :: hilf, hilf2
CHARACTER(255)            :: TimeAverageFile
REAL                      :: dimFactorSF
INTEGER                   :: nTotalDOF
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND)   :: MPISharedSize
INTEGER                   :: SendNodeCount, GlobalElemNode, GlobalElemRank, iNode, iProc, jElem, NonUniqueNodeID
INTEGER                   :: UniqueNodeID, TestElemID
LOGICAL,ALLOCATABLE       :: NodeDepoMapping(:,:)
INTEGER                   :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1),firstNode,lastNode
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

! Initialize Deposition
!CALL InitDepositionMethod()

!--- Allocate arrays for charge density collection and initialize
#if USE_MPI
MPISharedSize = INT(4*(PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/4*(PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems/),PartSource_Shared_Win,PartSource_Shared)
CALL MPI_WIN_LOCK_ALL(0,PartSource_Shared_Win,IERROR)
PartSource(1:4,0:PP_N,0:PP_N,0:PP_N,1:nComputeNodeTotalElems) => PartSource_Shared(1:4*(PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems)
ALLOCATE(PartSourceProc(   1:4,0:PP_N,0:PP_N,0:PP_N,1:nSendShapeElems))
#else
ALLOCATE(PartSource(1:4,0:PP_N,0:PP_N,0:PP_N,nElems))
#endif
#if USE_MPI
IF(myComputeNodeRank.EQ.0) THEN
#endif
  PartSource=0.
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(PartSource_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif
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
    CALL abort(&
__STAMP__&
,'ERROR in pic_depo.f90: Cannot allocate PartSourceOld!')
  END IF
  PartSourceOld=0.
  OutputSource = .TRUE.
ELSE
  OutputSource = GETLOGICAL('PIC-OutputSource','F')
END IF

!--- check if chargedensity is computed from TimeAverageFile
TimeAverageFile = GETSTR('PIC-TimeAverageFile','none')
IF (TRIM(TimeAverageFile).NE.'none') THEN
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


! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical', 'shape_function_spherical'
IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
  r_sf                  = GETREAL('PIC-shapefunction-radius')
  alpha_sf              = GETINT('PIC-shapefunction-alpha')
  DoSFEqui              = GETLOGICAL('PIC-shapefunction-equi')
  DoSFLocalDepoAtBounds = GETLOGICAL('PIC-shapefunction-local-depo-BC')
  DoSFChargeCons        = GETLOGICAL('PIC-shapefunction-charge-conservation')
  r2_sf = r_sf * r_sf  ! Radius squared
  r2_sf_inv = 1./r2_sf ! Inverse of radius squared


  IF(DoSFLocalDepoAtBounds)THEN ! init cell vol weight
    IF(.NOT.TRIM(DepositionType).EQ.'shape_function_2d') CALL abort(&
        __STAMP__&
        ,' PIC-shapefunction-local-depo-BC=T is currently only implemented for shape_function_2d!')
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
  IF ((TRIM(InterpolationType).NE.'cell_volweight')) THEN
    ALLOCATE(CellVolWeightFac(0:PP_N))
    CellVolWeightFac(0:PP_N) = xGP(0:PP_N)
    CellVolWeightFac(0:PP_N) = (CellVolWeightFac(0:PP_N)+1.0)/2.0
  END IF

  ! Initialize sub-cell volumes around nodes 
  CALL CalcCellLocNodeVolumes()
#if USE_MPI
  MPISharedSize = INT(4*nUniqueGlobalNodes,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/4,nUniqueGlobalNodes/),NodeSource_Shared_Win,NodeSource_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeSource_Shared_Win,IERROR)
  NodeSource => NodeSource_Shared
  ALLOCATE(NodeSourceLoc(1:4,1:nUniqueGlobalNodes))

  IF(DoDielectricSurfaceCharge)THEN

    firstNode = INT(REAL( myComputeNodeRank   *nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))+1
    lastNode  = INT(REAL((myComputeNodeRank+1)*nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))

   ! Global, synchronized surface charge contribution (is added to NodeSource AFTER MPI synchronization)
    MPISharedSize = INT(nUniqueGlobalNodes,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/1,nUniqueGlobalNodes/),NodeSourceExt_Shared_Win,NodeSourceExt_Shared)
    CALL MPI_WIN_LOCK_ALL(0,NodeSourceExt_Shared_Win,IERROR)
    NodeSourceExt => NodeSourceExt_Shared
    !ALLOCATE(NodeSourceExtLoc(1:1,1:nUniqueGlobalNodes))
    DO iNode=firstNode, lastNode
      NodeSourceExt(1,iNode) = 0.
    END DO
    CALL MPI_WIN_SYNC(NodeSourceExt_Shared_Win,IERROR)
    CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

   ! Local, non-synchronized surface charge contribution (is added to NodeSource BEFORE MPI synchronization)
    MPISharedSize = INT(nUniqueGlobalNodes,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/1,nUniqueGlobalNodes/),NodeSourceExtTmp_Shared_Win,NodeSourceExtTmp_Shared)
    CALL MPI_WIN_LOCK_ALL(0,NodeSourceExtTmp_Shared_Win,IERROR)
    NodeSourceExtTmp => NodeSourceExtTmp_Shared
    ALLOCATE(NodeSourceExtTmpLoc(1:1,1:nUniqueGlobalNodes))
    NodeSourceExtTmpLoc = 0.

    DO iNode=firstNode, lastNode
      NodeSourceExtTmp(1,iNode) = 0.
    END DO
    CALL MPI_WIN_SYNC(NodeSourceExtTmp_Shared_Win,IERROR)
    CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
    
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
                  , RecvRequest(iProc)                                          &
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
                    , SendRequest(iProc)                                          &
                    , IERROR)
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        ALLOCATE(NodeMapping(iProc)%RecvNodeUniqueGlobalID(1:NodeMapping(iProc)%nRecvUniqueNodes), &
              NodeMapping(iProc)%RecvNodeSource(1:4,1:NodeMapping(iProc)%nRecvUniqueNodes))
        CALL MPI_IRECV( NodeMapping(iProc)%RecvNodeUniqueGlobalID                   &
                      , NodeMapping(iProc)%nRecvUniqueNodes                         &
                      , MPI_INTEGER                                                 &
                      , iProc                                                       &
                      , 666                                                         &
                      , MPI_COMM_LEADERS_SHARED                                       &
                      , RecvRequest(iProc)                                          &
                      , IERROR)
      END IF
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        ALLOCATE(NodeMapping(iProc)%SendNodeUniqueGlobalID(1:NodeMapping(iProc)%nSendUniqueNodes), &
              NodeMapping(iProc)%SendNodeSource(1:4,1:NodeMapping(iProc)%nSendUniqueNodes))
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
                      , MPI_COMM_LEADERS_SHARED                                       &
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
  END IF
#else
  ALLOCATE(NodeSource(1:4,1:nUniqueGlobalNodes))
  IF(DoDielectricSurfaceCharge)THEN
    ALLOCATE(NodeSourceExt(1:1,1:nUniqueGlobalNodes))
    ALLOCATE(NodeSourceExtTmp(1:1,1:nUniqueGlobalNodes))
    NodeSourceExt    = 0.
    NodeSourceExtTmp = 0.
  END IF ! DoDielectricSurfaceCharge
#endif /*USE_MPI*/

  IF(DoDielectricSurfaceCharge)THEN
    ! Allocate and determine Vandermonde mapping from equidistant (visu) to NodeType node set                                          
    ALLOCATE(Vdm_EQ_N(0:PP_N,0:1))                                                                                                     
    CALL GetVandermonde(1, NodeTypeVISU, PP_N, NodeType, Vdm_EQ_N, modal=.FALSE.) 
  END IF ! DoDielectricSurfaceCharge


!  ! Additional source for cell_volweight_mean (external or surface charge)
!  IF(DoDielectricSurfaceCharge)THEN
!    ALLOCATE(NodeSourceExt(1:nNodes))
!    NodeSourceExt = 0.0
!    ALLOCATE(NodeSourceExtTmp(1:nNodes))
!    NodeSourceExtTmp = 0.0
!  END IF ! DoDielectricSurfaceCharge
CASE('shape_function')
  !ALLOCATE(PartToFIBGM(1:6,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) CALL abort(&
  !    __STAMP__&
  !    ' Cannot allocate PartToFIBGM!')
  !ALLOCATE(ExtPartToFIBGM(1:6,1:PDM%ParticleVecLength),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) THEN
  !  CALL abort(__STAMP__&
  !    ' Cannot allocate ExtPartToFIBGM!')
  BetaFac = beta(1.5, REAL(alpha_sf) + 1.)
  w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) &
                        * (REAL(alpha_sf) + 1.)/(PI*(r_sf**3))

  !-- ResampleAnalyzeSurfCollis
  SFResampleAnalyzeSurfCollis = GETLOGICAL('PIC-SFResampleAnalyzeSurfCollis','.FALSE.')
  IF (SFResampleAnalyzeSurfCollis) THEN
    LastAnalyzeSurfCollis%PartNumberSamp = 0
    LastAnalyzeSurfCollis%PartNumberDepo = 0
    LastAnalyzeSurfCollis%ReducePartNumber = GETLOGICAL('PIC-SFResampleReducePartNumber','.FALSE.')
    LastAnalyzeSurfCollis%PartNumThreshold = GETINT('PIC-PartNumThreshold','0')
    IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN
      WRITE(UNIT=hilf,FMT='(I0)') LastAnalyzeSurfCollis%PartNumThreshold
      LastAnalyzeSurfCollis%PartNumberReduced = GETINT('PIC-SFResamplePartNumberReduced',TRIM(hilf)) !def. PartNumThreshold
    END IF
    WRITE(UNIT=hilf,FMT='(E16.8)') -HUGE(1.0)
    LastAnalyzeSurfCollis%Bounds(1,1)   = MAX(GETREAL('PIC-SFResample-xmin',TRIM(hilf)),GEO%xmin-r_sf)
    LastAnalyzeSurfCollis%Bounds(1,2)   = MAX(GETREAL('PIC-SFResample-ymin',TRIM(hilf)),GEO%ymin-r_sf)
    LastAnalyzeSurfCollis%Bounds(1,3)   = MAX(GETREAL('PIC-SFResample-zmin',TRIM(hilf)),GEO%zmin-r_sf)
    WRITE(UNIT=hilf,FMT='(E16.8)') HUGE(1.0)
    LastAnalyzeSurfCollis%Bounds(2,1)   = MIN(GETREAL('PIC-SFResample-xmax',TRIM(hilf)),GEO%xmax+r_sf)
    LastAnalyzeSurfCollis%Bounds(2,2)   = MIN(GETREAL('PIC-SFResample-ymax',TRIM(hilf)),GEO%ymax+r_sf)
    LastAnalyzeSurfCollis%Bounds(2,3)   = MIN(GETREAL('PIC-SFResample-zmax',TRIM(hilf)),GEO%zmax+r_sf)
    LastAnalyzeSurfCollis%UseFixBounds = GETLOGICAL('PIC-SFResample-UseFixBounds','.TRUE.')
    LastAnalyzeSurfCollis%NormVecOfWall = GETREALARRAY('PIC-NormVecOfWall',3,'1. , 0. , 0.')  !directed outwards
    IF (DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%NormVecOfWall).GT.0.) THEN
      LastAnalyzeSurfCollis%NormVecOfWall = LastAnalyzeSurfCollis%NormVecOfWall &
        / SQRT( DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%NormVecOfWall) )
    END IF
    LastAnalyzeSurfCollis%Restart = GETLOGICAL('PIC-SFResampleRestart','.FALSE.')
    IF (LastAnalyzeSurfCollis%Restart) THEN
      LastAnalyzeSurfCollis%DSMCSurfCollisRestartFile = GETSTR('PIC-SFResampleRestartFile','dummy')
    END IF
    !-- BCs
    LastAnalyzeSurfCollis%NumberOfBCs = GETINT('PIC-SFResampleNumberOfBCs','1')
    SDEALLOCATE(LastAnalyzeSurfCollis%BCs)
    ALLOCATE(LastAnalyzeSurfCollis%BCs(1:LastAnalyzeSurfCollis%NumberOfBCs),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate LastAnalyzeSurfCollis%BCs!')
    END IF
    IF (LastAnalyzeSurfCollis%NumberOfBCs.EQ.1) THEN !already allocated
      LastAnalyzeSurfCollis%BCs = GETINTARRAY('PIC-SFResampleSurfCollisBC',1,'0') ! 0 means all...
    ELSE
      hilf2=''
      DO iBC=1,LastAnalyzeSurfCollis%NumberOfBCs !build default string: 0,0,0,...
        WRITE(UNIT=hilf,FMT='(I0)') 0
        hilf2=TRIM(hilf2)//TRIM(hilf)
        IF (iBC.NE.LastAnalyzeSurfCollis%NumberOfBCs) hilf2=TRIM(hilf2)//','
      END DO
      LastAnalyzeSurfCollis%BCs = GETINTARRAY('PIC-SFResampleSurfCollisBC',LastAnalyzeSurfCollis%NumberOfBCs,hilf2)
    END IF
    !-- spec for dt-calc
    LastAnalyzeSurfCollis%NbrOfSpeciesForDtCalc = GETINT('PIC-SFResampleNbrOfSpeciesForDtCalc','1')
    SDEALLOCATE(LastAnalyzeSurfCollis%SpeciesForDtCalc)
    ALLOCATE(LastAnalyzeSurfCollis%SpeciesForDtCalc(1:LastAnalyzeSurfCollis%NbrOfSpeciesForDtCalc),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate LastAnalyzeSurfCollis%SpeciesForDtCalc!')
    END IF
    IF (LastAnalyzeSurfCollis%NbrOfSpeciesForDtCalc.EQ.1) THEN !already allocated
      LastAnalyzeSurfCollis%SpeciesForDtCalc = GETINTARRAY('PIC-SFResampleSpeciesForDtCalc',1,'0') ! 0 means all...
    ELSE
      hilf2=''
      DO iBC=1,LastAnalyzeSurfCollis%NbrOfSpeciesForDtCalc !build default string: 0,0,0,...
        WRITE(UNIT=hilf,FMT='(I0)') 0
        hilf2=TRIM(hilf2)//TRIM(hilf)
        IF (iBC.NE.LastAnalyzeSurfCollis%NbrOfSpeciesForDtCalc) hilf2=TRIM(hilf2)//','
      END DO
      LastAnalyzeSurfCollis%SpeciesForDtCalc &
        = GETINTARRAY('PIC-SFResampleSpeciesForDtCalc',LastAnalyzeSurfCollis%NbrOfSpeciesForDtCalc,hilf2)
    END IF
  END IF

  VolumeShapeFunction=4./3.*PI*r_sf**3
  nTotalDOF=nGlobalElems*(PP_N+1)**3
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.MeshVolume) &
      CALL abort(&
      __STAMP__&
      ,'ShapeFunctionVolume > MeshVolume')
  END IF

  CALL PrintOption('Average DOFs in Shape-Function','CALCUL.',RealOpt=REAL(nTotalDOF)*VolumeShapeFunction/MeshVolume)

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1.
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

CASE('shape_function_1d','shape_function_2d')
  ! Get deposition direction for 1D or perpendicular direction for 2D
  sf1d_dir = GETINT ('PIC-shapefunction1d-direction')
  ! Distribute the charge over the volume (3D) or line (1D)/area (2D): default is TRUE
  sfDepo3D = GETLOGICAL('PIC-shapefunction-3D-deposition')
  hilf2='volume'

  SELECT CASE(TRIM(DepositionType))
  CASE('shape_function_1d')
    ! Set perpendicular directions
    IF(sf1d_dir.EQ.1)THEN ! Shape function deposits charge in x-direction
      dimFactorSF = (GEO%ymaxglob-GEO%yminglob)*(GEO%zmaxglob-GEO%zminglob)
    ELSE IF (sf1d_dir.EQ.2)THEN ! Shape function deposits charge in y-direction
      dimFactorSF = (GEO%xmaxglob-GEO%xminglob)*(GEO%zmaxglob-GEO%zminglob)
    ELSE IF (sf1d_dir.EQ.3)THEN ! Shape function deposits charge in z-direction
      dimFactorSF = (GEO%xmaxglob-GEO%xminglob)*(GEO%ymaxglob-GEO%yminglob)
    ELSE
      w_sf=2*MeshVolume
    END IF
    IF(sfDepo3D)THEN ! Distribute the charge over the volume (3D)
      ! Set prefix factor
      w_sf = GAMMA(REAL(alpha_sf)+1.5)/(SQRT(PI)*r_sf*GAMMA(REAL(alpha_sf+1))*dimFactorSF)
    ELSE ! Distribute the charge over the line (1D)
      ! Set prefix factor
      w_sf = GAMMA(REAL(alpha_sf)+1.5)/(SQRT(PI)*r_sf*GAMMA(REAL(alpha_sf+1)))
      ! Set shape function length (1D volume)
      !VolumeShapeFunction=2*r_sf
      hilf2='line'
    END IF
    ! Set shape function length (3D volume)
    VolumeShapeFunction=2*r_sf*dimFactorSF
    ! Calculate number of 1D DOF (assume second and third direction with 1 cell layer and area given by dimFactorSF)
    nTotalDOF=nGlobalElems*(PP_N+1)
    hilf='1D'
  CASE('shape_function_2d')
    ! Set perpendicular direction
    IF(sf1d_dir.EQ.1)THEN ! Shape function deposits charge in y-z-direction (const. in x)
      dimFactorSF = (GEO%xmaxglob-GEO%xminglob)
    ELSE IF (sf1d_dir.EQ.2)THEN ! Shape function deposits charge in x-z-direction (const. in y)
      dimFactorSF = (GEO%ymaxglob-GEO%yminglob)
    ELSE IF (sf1d_dir.EQ.3)THEN! Shape function deposits charge in x-y-direction (const. in z)
      dimFactorSF = (GEO%zmaxglob-GEO%zminglob)
    ELSE
      w_sf=2*MeshVolume
    END IF
    IF(sfDepo3D)THEN ! Distribute the charge over the volume (3D)
      ! Set prefix factor
      w_sf = (REAL(alpha_sf)+1.0)/(PI*r2_sf*dimFactorSF)
    ELSE ! Distribute the charge over the area (2D)
      ! Set prefix factor
      w_sf = (REAL(alpha_sf)+1.0)/(PI*r2_sf)
      ! Set shape function length (2D volume)
      !VolumeShapeFunction=PI*(r_sf**2)
      hilf2='area'
    END IF
    ! Set shape function length (3D volume)
    VolumeShapeFunction=PI*(r_sf**2)*dimFactorSF
    ! Calculate number of 2D DOF (assume third direction with 1 cell layer and width dimFactorSF)
    nTotalDOF=nGlobalElems*(PP_N+1)**2
    hilf='2D'
  END SELECT

  ASSOCIATE(nTotalDOFin3D             => nGlobalElems*(PP_N+1)**3   ,&
            VolumeShapeFunctionSphere => 4./3.*PI*r_sf**3           )
    SWRITE(UNIT_stdOut,'(A)') ' The complete charge is '//TRIM(hilf2)//' distributed (deposition function is '//TRIM(hilf)//')'
    IF(.NOT.sfDepo3D)THEN
      SWRITE(UNIT_stdOut,'(A)') ' Note that the integral of the charge density over the volume is larger than the complete charge!'
    END IF
    IF(TRIM(DepositionType).EQ.'shape_function_1d')THEN
      CALL PrintOption('Shape function volume (3D box)', 'CALCUL.', RealOpt=VolumeShapeFunction)
    ELSE
      CALL PrintOption('Shape function volume (3D cylinder)', 'CALCUL.', RealOpt=VolumeShapeFunction)
    END IF
    !CALL PrintOption('Shape function volume ('//TRIM(hilf)//')'                , 'CALCUL.' , RealOpt=VolumeShapeFunction)
    IF(MPIRoot)THEN
      IF(VolumeShapeFunction.GT.MeshVolume)THEN
        CALL PrintOption('Mesh volume ('//TRIM(hilf)//')'                , 'CALCUL.' , RealOpt=MeshVolume)
        WRITE(UNIT_stdOut,'(A)') ' Maybe wrong perpendicular direction (PIC-shapefunction1d-direction)?'
        CALL abort(&
        __STAMP__&
        ,'ShapeFunctionVolume > MeshVolume ('//TRIM(hilf)//' shape function)')
      END IF
    END IF
    CALL PrintOption('Average DOFs in Shape-Function '//TRIM(hilf2)//' ('//TRIM(hilf)//')'       , 'CALCUL.' , RealOpt=&
         REAL(nTotalDOF)*VolumeShapeFunction/MeshVolume)
    CALL PrintOption('Average DOFs in Shape-Function (corresponding 3D sphere)' , 'CALCUL.' , RealOpt=&
         REAL(nTotalDOFin3D)*VolumeShapeFunctionSphere/MeshVolume)
  END ASSOCIATE

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1.
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

CASE('shape_function_cylindrical','shape_function_spherical')
  !IF(.NOT.DoRefMapping) CALL abort(&
  !  __STAMP__&
  !  ,' Shape function has to be used with ref element tracking.')
  !ALLOCATE(PartToFIBGM(1:6,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) CALL abort(&
  !    __STAMP__&
  !    ' Cannot allocate PartToFIBGM!')
  !ALLOCATE(ExtPartToFIBGM(1:6,1:PDM%ParticleVecLength),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) THEN
  !  CALL abort(__STAMP__&
  !    ' Cannot allocate ExtPartToFIBGM!')
  IF(TRIM(DepositionType).EQ.'shape_function_cylindrical')THEN
    SfRadiusInt=2
  ELSE
    SfRadiusInt=3
  END IF
  r_sf0      = GETREAL   ('PIC-shapefunction-radius0','1.')
  r_sf_scale = GETREAL   ('PIC-shapefunction-scale','0.')
  BetaFac    = beta(1.5, REAL(alpha_sf) + 1.)
  w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) &
                        * (REAL(alpha_sf) + 1.)!/(PI*(r_sf**3))

  r_sf_average=0.5*(r_sf+r_sf0)
  VolumeShapeFunction=4./3.*PI*r_sf_average*r_sf_average
  nTotalDOF=nGlobalElems*(PP_N+1)**3
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.MeshVolume) &
      CALL abort(&
      __STAMP__&
      ,'ShapeFunctionVolume > MeshVolume')
  END IF

  CALL PrintOption('Average DOFs in Shape-Function','CALCUL.',RealOpt=REAL(nTotalDOF)*VolumeShapeFunction/MeshVolume)

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1.
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Unknown DepositionType in pic_depo.f90')
END SELECT

IF (PartSourceConstExists) THEN
  ALLOCATE(PartSourceConst(1:4,0:PP_N,0:PP_N,0:PP_N,nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
__STAMP__&
,'ERROR in pic_depo.f90: Cannot allocate PartSourceConst!')
  END IF
  PartSourceConst=0.
END IF

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition


SUBROUTINE Deposition(doParticle_In)
!============================================================================================================================
! This subroutine performs the deposition of the particle charge and current density to the grid
! following list of distribution methods are implemented
! - shape function       (only one type implemented)
! useVMPF added, therefore, this routine contains automatically the use of variable mpfs
!============================================================================================================================
! USE MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Analyze_Vars       ,ONLY: DoVerifyCharge,PartAnalyzeStep
USE MOD_Particle_Vars
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Method              ,ONLY: DepositionMethod
USE MOD_PIC_Analyze                 ,ONLY: VerifyDepositedCharge
USE MOD_TimeDisc_Vars               ,ONLY: iter
#if USE_MPI
USE MOD_MPI_Shared_Vars             ,ONLY: myComputeNodeRank,MPI_COMM_SHARED
#endif  /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT variable declaration
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
!============================================================================================================================
! Return, if no deposition is required
IF(.NOT.DoDeposition) RETURN

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif  /*USE_MPI*/
  PartSource = 0.0
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(PartSource_Shared_Win, IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
#endif  /*USE_MPI*/

IF(PRESENT(doParticle_In)) THEN
  CALL DepositionMethod(doParticle_In)
ELSE
  CALL DepositionMethod()
END IF

IF(MOD(iter,PartAnalyzeStep).EQ.0) THEN
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF
RETURN
END SUBROUTINE Deposition


SUBROUTINE FinalizeDeposition()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic deposition
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars ,ONLY: Geo
USE MOD_Dielectric_Vars    ,ONLY: DoDielectricSurfaceCharge
#if USE_MPI
USE MOD_MPI_Shared_vars    ,ONLY: MPI_COMM_SHARED
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
SDEALLOCATE(GaussBorder)
SDEALLOCATE(ElemDepo_xGP)
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

#if USE_MPI
SDEALLOCATE(PartSourceProc)

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

IF(DoDeposition)THEN
  CALL MPI_WIN_UNLOCK_ALL(PartSource_Shared_Win, iError)
  CALL MPI_WIN_FREE(      PartSource_Shared_Win, iError)

  ! Deposition-dependent arrays
  SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean')
    CALL MPI_WIN_UNLOCK_ALL(NodeSource_Shared_Win, iError)
    CALL MPI_WIN_FREE(      NodeSource_Shared_Win, iError)
    ADEALLOCATE(NodeSource_Shared)
    ! Surface charging arrays
    IF(DoDielectricSurfaceCharge)THEN
      CALL MPI_WIN_UNLOCK_ALL(NodeSourceExt_Shared_Win, iError)
      CALL MPI_WIN_FREE(      NodeSourceExt_Shared_Win, iError)
      ADEALLOCATE(NodeSourceExt_Shared)
    END IF ! DoDielectricSurfaceCharge
  END SELECT

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
END IF ! DoDeposition

! Then, free the pointers or arrays
ADEALLOCATE(PartSource_Shared)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
ADEALLOCATE(PartSource)

! Deposition-dependent pointers/arrays
SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean')
    ADEALLOCATE(NodeSource)
  ! Surface charging pointers/arrays
  IF(DoDielectricSurfaceCharge)THEN
      ADEALLOCATE(NodeSourceExt)
  END IF ! DoDielectricSurfaceCharge
END SELECT


END SUBROUTINE FinalizeDeposition

END MODULE MOD_PICDepo
