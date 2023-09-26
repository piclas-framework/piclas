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

MODULE MOD_Particle_Mesh_Build
!===================================================================================================================================
! Contains all build routines required for particle mesh
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

PUBLIC:: BuildElementRadiusTria,BuildElemTypeAndBasisTria,BuildEpsOneCell,BuildBCElemDistance
PUBLIC:: BuildNodeNeighbourhood,BuildElementOriginShared,BuildElementBasisAndRadius
PUBLIC:: BuildSideOriginAndRadius,BuildLinearSideBaseVectors, BuildMesh2DInfo
PUBLIC:: BuildSideSlabAndBoundingBox
!===================================================================================================================================

CONTAINS

SUBROUTINE BuildMesh2DInfo()
!===================================================================================================================================
!> Routine determines a symmetry side and calculates the 2D (area faces in symmetry plane) and axisymmetric volumes (cells are
!> revolved around the symmetry axis). The symmetry side will be used later on to determine in which two directions the quadtree
!> shall refine the mesh, skipping the z-dimension to avoid an unnecessary refinement.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared, SideInfo_Shared
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
#if USE_MPI
USE MOD_MPI_Shared             
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED, nComputeNodeTotalElems
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemSideNodeID2D_Shared_Win, SideNormalEdge2D_Shared_Win
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank, nComputeNodeProcessors
#else
USE MOD_Mesh_Vars               ,ONLY: nElems
#endif /*USE_MPI*/
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemSideNodeID2D_Shared, SideNormalEdge2D_Shared
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, iLocSide, iNode, iELem
REAL                            :: VecCell(2), FaceMidPoint(2), NormVec(2), EdgeVec(2), nVal
INTEGER                         :: firstElem,lastElem, GlobalElemID, tmpNode
LOGICAL                         :: DefineSide
!===================================================================================================================================
#if USE_MPI
CALL Allocate_Shared((/2,6,nComputeNodeTotalElems/),ElemSideNodeID2D_Shared_Win,ElemSideNodeID2D_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemSideNodeID2D_Shared_Win,IERROR)
CALL Allocate_Shared((/4,6,nComputeNodeTotalElems/),SideNormalEdge2D_Shared_Win,SideNormalEdge2D_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideNormalEdge2D_Shared_Win,IERROR)

firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(ElemSideNodeID2D_Shared(1:2,1:6,1:nElems))
ALLOCATE(SideNormalEdge2D_Shared(1:4,1:6,1:nElems))
firstElem = 1
lastElem  = nElems
#endif

DO iElem = firstElem, lastElem
  GlobalElemID = GetGlobalElemID(iElem)
  DO iLocSide = 1, 6   
    DefineSide = .TRUE.
    SideID=GetGlobalNonUniqueSideID(GlobalElemID,iLocSide)
    IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
      IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%SymmetryBC) THEN
        ElemSideNodeID2D_Shared(:,iLocSide, iElem) = -1
        SideNormalEdge2D_Shared(:,iLocSide, iElem) = 0.
        DefineSide = .FALSE.
      END IF
    END IF
    
    IF (DefineSide) THEN
      tmpNode = 0 
      DO iNode = 1, 4
        IF(NodeCoords_Shared(3,ElemSideNodeID_Shared(iNode,iLocSide,iElem)+1).GT.(GEO%zmaxglob+GEO%zminglob)/2.) THEN
          tmpNode = tmpNode + 1
          ElemSideNodeID2D_Shared(tmpNode,iLocSide, iElem) = ElemSideNodeID_Shared(iNode,iLocSide,iElem)+1
        END IF
      END DO      
      EdgeVec(1:2) = NodeCoords_Shared(1:2,ElemSideNodeID2D_Shared(2,ilocSide,iElem))-NodeCoords_Shared(1:2,ElemSideNodeID2D_Shared(1,ilocSide,iElem))
      NormVec(1) = -EdgeVec(2)
      NormVec(2) = EdgeVec(1) 
      FaceMidPoint(1:2) = (NodeCoords_Shared(1:2,ElemSideNodeID_Shared(1,iLocSide,iElem)+1) &
                        + NodeCoords_Shared(1:2,ElemSideNodeID_Shared(2,iLocSide,iElem)+1)&
                        + NodeCoords_Shared(1:2,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
                        + NodeCoords_Shared(1:2,ElemSideNodeID_Shared(4,iLocSide,iElem)+1)) / 4.
      VecCell(1:2) = FaceMidPoint(1:2) - ElemBaryNGeo(1:2,iElem) ! vector from elem bary to side
      IF (DOT_PRODUCT(VecCell,NormVec).GT.0.0) THEN
        SideNormalEdge2D_Shared(1:2,iLocSide, iElem) = NormVec(1:2)
      ELSE
        SideNormalEdge2D_Shared(1,iLocSide, iElem) = EdgeVec(2)
        SideNormalEdge2D_Shared(2,iLocSide, iElem) = -EdgeVec(1)
      END IF
      nVal = SQRT(SideNormalEdge2D_Shared(1,iLocSide,iElem)**2. + SideNormalEdge2D_Shared(2,iLocSide,iElem)**2.)
      SideNormalEdge2D_Shared(1,iLocSide,iElem)= SideNormalEdge2D_Shared(1,iLocSide,iElem)/ nVal
      SideNormalEdge2D_Shared(2,iLocSide,iElem)= SideNormalEdge2D_Shared(2,iLocSide,iElem) / nVal
      nVal = SQRT(EdgeVec(1)**2. + EdgeVec(2)**2.)
      SideNormalEdge2D_Shared(3,iLocSide,iElem) = EdgeVec(1) / nVal
      SideNormalEdge2D_Shared(4,iLocSide,iElem)  = EdgeVec(2) / nVal
    END IF
  END DO
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemSideNodeID2D_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideNormalEdge2D_Shared_Win,MPI_COMM_SHARED)
#endif

END SUBROUTINE BuildMesh2DInfo

SUBROUTINE BuildElementRadiusTria()
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo,ElemRadius2NGeo, ElemRadiusNGeo
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo_Shared,ElemRadius2NGeo_Shared,ElemRadiusNGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared_Win,ElemRadius2NGeo_Shared_Win, ElemRadiusNGeo_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: nELems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,iNode
REAL                           :: xPos(3),Radius
INTEGER                        :: firstElem, lastElem
!================================================================================================================================

#if USE_MPI
  CALL Allocate_Shared((/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)
  CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemRadius2NGeo_Shared_Win,ElemRadius2NGEO_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadius2NGeo_Shared_Win,IERROR)
  ElemRadius2NGeo    => ElemRadius2NGeo_Shared
  ElemBaryNGeo       => ElemBaryNGeo_Shared
  IF(StringBeginsWith(DepositionType,'shape_function'))THEN
    CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemRadiusNGeo_Shared_Win,ElemRadiusNGEO_Shared)
    CALL MPI_WIN_LOCK_ALL(0,ElemRadiusNGeo_Shared_Win,IERROR)
    ElemRadiusNGeo    => ElemRadiusNGeo_Shared
  END IF

#else
ALLOCATE(ElemBaryNGeo(1:3,nElems) &
        ,ElemRadius2NGeo( nElems))
IF(StringBeginsWith(DepositionType,'shape_function')) ALLOCATE(ElemRadiusNGeo(nElems))
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemRadius2NGeo = 0.
  IF(StringBeginsWith(DepositionType,'shape_function')) ElemRadiusNGeo = 0.
#if USE_MPI
END IF
IF(StringBeginsWith(DepositionType,'shape_function')) CALL BARRIER_AND_SYNC(ElemRadiusNGeo_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemRadius2NGeo_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI*/

#if USE_MPI
firstElem=INT(REAL(myComputeNodeRank    )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif /*USE_MPI*/

DO iElem=firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  Radius=0.
  xPos  =0.

  DO iNode=1,8
    xPos = xPos + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode)
  END DO
    ElemBaryNGeo(:,iElem) = xPos/8.
  DO iNode=1,8
    xPos   = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode) - ElemBaryNGeo(:,iElem)
    Radius = MAX(Radius,VECNORM(xPos))
  END DO
  ElemRadius2NGeo(iElem) = Radius*Radius
  IF(StringBeginsWith(DepositionType,'shape_function')) ElemRadiusNGeo(iElem) = Radius
END DO ! iElem

#if USE_MPI
IF(StringBeginsWith(DepositionType,'shape_function')) CALL BARRIER_AND_SYNC(ElemRadiusNGeo_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemRadius2NGeo_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemBaryNGeo_Shared_Win   ,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementRadiusTria


SUBROUTINE BuildElemTypeAndBasisTria()
!===================================================================================================================================
!> Dummy routine to fill the ElemCurved array with TriaTracking
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis                   ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Mesh_Vars               ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars               ,ONLY: nElems
USE MOD_Mesh_Vars               ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Photon_TrackingVars     ,ONLY: UsePhotonTriaTracking
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,iDir
INTEGER                        :: i,j,k
REAL                           :: xPos(3)
REAL                           :: Xi(3,6),Lag(1:3,0:NGeo)
INTEGER                        :: firstElem, lastElem
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved ...'

! elements
#if USE_MPI
IF(UsePhotonTriaTracking)THEN
  CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
END IF ! UsePhotonTriaTracking
CALL Allocate_Shared((/3,6,nComputeNodeTotalElems/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
CALL Allocate_Shared((/6,nComputeNodeTotalElems/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
IF(UsePhotonTriaTracking)THEN
  ElemCurved         => ElemCurved_Shared
END IF ! UsePhotonTriaTracking
XiEtaZetaBasis     => XiEtaZetaBasis_Shared
slenXiEtaZetaBasis => slenXiEtaZetaBasis_Shared

ASSOCIATE(XCL_NGeo     => XCL_NGeo_Shared)

#else
IF(UsePhotonTriaTracking)THEN
  ALLOCATE(ElemCurved(            1:nElems))
END IF ! UsePhotonTriaTracking
ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:nElems))
ALLOCATE(slenXiEtaZetaBasis(1:6,1:nElems))
#endif /*USE_MPI*/

IF(UsePhotonTriaTracking)THEN
  ! only CN root nullifies
#if USE_MPI
  IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
    ElemCurved   = .FALSE.
#if USE_MPI
  END IF
#endif /*USE_MPI*/
END IF ! UsePhotonTriaTracking

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos = 0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos = xPos+XCL_NGeo(:,i,j,k,ElemID)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    XiEtaZetaBasis(1:3,iDir,iElem) = xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(:,iDir,iElem)   = XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem) = 1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
  END DO ! iDir = 1, 6
END DO

#if USE_MPI
END ASSOCIATE
IF(UsePhotonTriaTracking)THEN
  CALL BARRIER_AND_SYNC(ElemCurved_Shared_Win      ,MPI_COMM_SHARED)
END IF ! UsePhotonTriaTracking
CALL BARRIER_AND_SYNC(XiEtaZetaBasis_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(slenXiEtaZetaBasis_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE BuildElemTypeAndBasisTria


SUBROUTINE BuildEpsOneCell()
!===================================================================================================================================
! Build epsOneCell for each element
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeCL,NodeType
USE MOD_Mesh_Vars              ,ONLY: sJ
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: RefMappingEps
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoRef
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: dXCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ_Shared,ElemEpsOneCell_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ_Shared_Win,ElemEpsOneCell_Shared_Win
#else
USE MOD_Particle_Mesh_Vars     ,ONLY: nComputeNodeElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,firstElem,lastElem
REAL                           :: scaleJ,maxScaleJ
#if USE_MPI
INTEGER                        :: ElemID
INTEGER                        :: i,j,k
! Vandermonde matrices
REAL                           :: Vdm_CLNGeo_NGeoRef(0:NgeoRef,0:NGeo)
REAL                           :: Vdm_NGeoRef_N(     0:PP_N   ,0:NGeoRef)
! Jacobian on CL N and NGeoRef
REAL                           :: detJac_Ref(1  ,0:NGeoRef,0:NGeoRef,0:NGeoRef)
REAL                           :: DetJac_N(  1  ,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on CL N
REAL                           :: dX_NGeoRef(3,3,0:NGeoRef,0:NGeoRef,0:NGeoRef)

INTEGER                        :: ElemLocID
#endif /*USE_MPI*/
REAL                           :: StartT,EndT
!===================================================================================================================================

IF(MPIRoot)THEN
  LBWRITE(UNIT_StdOut,'(132("-"))')
  LBWRITE(UNIT_StdOut,'(A)',ADVANCE='NO') ' Building EpsOneCell for all elements ...'
  GETTIME(StartT)
END IF ! MPIRoot

! build sJ for all elements not on local proc
#if USE_MPI
CALL Allocate_Shared((/(PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems/),ElemsJ_Shared_Win,ElemsJ_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemsJ_Shared_Win,IERROR)
ElemsJ(0:PP_N,0:PP_N,0:PP_N,1:nComputeNodeTotalElems) => ElemsJ_Shared

IF (myComputeNodeRank.EQ.0) THEN
  ElemsJ_Shared = 0.
END IF

CALL BARRIER_AND_SYNC(ElemsJ_Shared_Win,MPI_COMM_SHARED)

firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

#if USE_MPI
! Calculate sJ for elements not inside current proc, otherwise copy local values
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNGeo_NGeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NGeoRef_N     , modal=.TRUE.)

DO iElem = firstElem,lastElem
  ElemID    = GetGlobalElemID(iElem)
  ElemLocID = ElemID-offsetElem
  ! element on local proc, sJ already calculated in metrics.f90
  IF ((ElemLocID.GT.0) .AND. (ElemLocID.LE.nElems)) THEN
    ElemsJ(:,:,:,iElem) = sJ(:,:,:,ElemLocID)

  ! element not on local proc, calculate sJ frm dXCL_NGeo_Shared
  ELSE
    detJac_Ref = 0.
    ! Compute Jacobian on NGeo and then interpolate:
    ! required to guarantee conservation when restarting with N<NGeo
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,1,:,:,:,ElemID),dX_NGeoRef(:,1,:,:,:))
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,2,:,:,:,ElemID),dX_NGeoRef(:,2,:,:,:))
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,3,:,:,:,ElemID),dX_NGeoRef(:,3,:,:,:))
    DO k=0,NGeoRef; DO j=0,NGeoRef; DO i=0,NGeoRef
      detJac_Ref(1,i,j,k)=detJac_Ref(1,i,j,k) &
        + dX_NGeoRef(1,1,i,j,k)*(dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(3,3,i,j,k) - dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(2,3,i,j,k))  &
        + dX_NGeoRef(2,1,i,j,k)*(dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(1,3,i,j,k) - dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(3,3,i,j,k))  &
        + dX_NGeoRef(3,1,i,j,k)*(dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(2,3,i,j,k) - dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(1,3,i,j,k))
    END DO; END DO; END DO !i,j,k=0,NgeoRef

    ! interpolate detJac_ref to the solution points
    CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:),DetJac_N)

    ! assign to global Variable sJ
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ElemsJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
    END DO; END DO; END DO !i,j,k=0,PP_N
  END IF
END DO

CALL BARRIER_AND_SYNC(ElemsJ_Shared_Win,MPI_COMM_SHARED)
#else
ElemsJ => sJ
#endif /* USE_MPI*/

! Exit routine here if TriaTracking is active
IF (TrackingMethod.EQ.TRIATRACKING)THEN
  IF(MPIRoot)THEN
    GETTIME(EndT)
    CALL DisplayMessageAndTime(EndT-StartT, 'DONE!')
  END IF ! MPIRoot
  RETURN
END IF

! allocate epsOneCell
#if USE_MPI
CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemEpsOneCell_Shared_Win,ElemEpsOneCell_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemEpsOneCell_Shared_Win,IERROR)
ElemEpsOneCell => ElemEpsOneCell_Shared
#else
ALLOCATE(ElemEpsOneCell(1:nComputeNodeElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemEpsOneCell = -1.
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(ElemEpsOneCell_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI*/

maxScaleJ = 0.
DO iElem = firstElem,lastElem
  scaleJ = MAXVAL(ElemsJ(:,:,:,iElem))/MINVAL(ElemsJ(:,:,:,iElem))
  ElemEpsOneCell(iElem) = 1.0 + SQRT(3.0*scaleJ*RefMappingEps)
  maxScaleJ  =MAX(scaleJ,maxScaleJ)
END DO ! iElem = firstElem,lastElem

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemEpsOneCell_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI*/

!IF(CalcMeshInfo)THEN
!  CALL AddToElemData(ElementOut,'epsOneCell',RealArray=epsOneCell(1:nElems))
!END IF
IF(MPIRoot)THEN
  GETTIME(EndT)
  CALL DisplayMessageAndTime(EndT-StartT, 'DONE!')
END IF ! MPIRoot

END SUBROUTINE BuildEpsOneCell


SUBROUTINE BuildBCElemDistance()
!===================================================================================================================================
! get the distance of each BC face
!> 1) identify BC elements to be handled by Tracing
!> 2) build mapping global elem ID to BC elem ID
!> 3) calculate distance from element origin to each side
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
#if USE_MPI
USE MOD_Globals_Vars           ,ONLY: c
#endif /*USE_MPI*/
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides,SideBCMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo,ElemRadiusNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueBCSides
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Utils                  ,ONLY: InsertionSort
#if USE_MPI
USE MOD_TimeDisc_Vars          ,ONLY: ManualTimeStep
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides_Shared,SideBCMetrics_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides_Shared_Win,SideBCMetrics_Shared_Win
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps,halo_eps_velo,SafetyFactor
#if ! (USE_HDG)
USE MOD_CalcTimeStep           ,ONLY: CalcTimeStep
#endif
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars          ,ONLY: nRKStages,RK_c
#endif
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: nComputeNodeElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: p
INTEGER                        :: ElemID,SideID
INTEGER                        :: iElem,firstElem,lastElem
INTEGER                        :: iSide,firstSide,lastSide
INTEGER                        :: nComputeNodeBCSides
INTEGER                        :: nBCSidesElem,nBCSidesProc,offsetBCSidesProc,offsetBCSides
INTEGER                        :: iBCSide,BCElemID,BCSideID
INTEGER                        :: CNElemID,BCCNElemID
REAL                           :: origin(1:3),vec(1:3)
REAL                           :: BC_halo_eps
LOGICAL                        :: fullMesh
REAL,ALLOCATABLE               :: tmpSideBCMetrics(:,:)
REAL,ALLOCATABLE               :: tmpSideBCDistance(:)
INTEGER,ALLOCATABLE            :: intSideBCMetrics(:)
#if USE_MPI
REAL                           :: BC_halo_eps_velo,BC_halo_diag,deltaT
INTEGER                        :: sendbuf,recvbuf
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
INTEGER                        :: iStage
#endif
#endif /*USE_MPI*/
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_StdOut,'(A)') ' Identifying BC sides and calculating side metrics ...'

! elements
#if USE_MPI
CALL Allocate_Shared((/2,nComputeNodeTotalElems/),ElemToBCSides_Shared_Win,ElemToBCSides_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBCSides_Shared_Win,IERROR)
ElemToBCSides => ElemToBCSides_Shared
#else
ALLOCATE(ElemToBCSides(1:2,1:nComputeNodeElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemToBCSides = -1
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(ElemToBCSides_Shared_Win,MPI_COMM_SHARED)

firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

! if running on one node, halo_eps is meaningless. Get a representative BC_halo_eps for BC side identification
fullMesh = .FALSE.
IF (halo_eps.EQ.0) THEN
  ! reconstruct halo_eps_velo
  IF (halo_eps_velo.EQ.0) THEN
    BC_halo_eps_velo = c
  ELSE
    BC_halo_eps_velo = halo_eps_velo
  END IF

  ! reconstruct deltaT
  deltaT = 0.
  IF (ManualTimeStep.GT.0.) THEN
    deltaT    = ManualTimeStep
#if ! (USE_HDG)
  ELSE
    deltaT    = CalcTimeStep()
#endif
  END IF

  ! calculate halo_eps
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  BC_halo_eps = RK_c(2)
  DO iStage=2,nRKStages-1
    BC_halo_eps = MAX(BC_halo_eps,RK_c(iStage+1)-RK_c(iStage))
  END DO
  BC_halo_eps = MAX(BC_halo_eps,1.-RK_c(nRKStages))
  BC_halo_eps = BC_halo_eps*BC_halo_eps_velo*deltaT*SafetyFactor
#else
  BC_halo_eps = BC_halo_eps_velo*deltaT*SafetyFactor
#endif

  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  BC_halo_diag = VECNORM(vec)

  ! compare halo_eps against global diagonal and reduce if necessary
  IF (.NOT.ALMOSTZERO(BC_halo_eps).AND.(BC_halo_diag.GE.BC_halo_eps)) THEN
    LBWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to ',BC_halo_eps
  ELSEIF (.NOT.ALMOSTZERO(BC_halo_eps).AND.(BC_halo_diag.LT.BC_halo_eps)) THEN
    fullMesh = .TRUE.
    BC_halo_eps = BC_halo_diag
    LBWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to global diag with ',BC_halo_eps
  ! halo_eps still at zero. Set it to global diagonal
  ELSE
    fullMesh = .TRUE.
    BC_halo_eps = BC_halo_diag
    LBWRITE(UNIT_stdOUt,'(A,F11.3)') ' | No halo_eps given and could not be reconstructed. Using global diag with ',BC_halo_eps
  END IF
ELSE
  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  BC_halo_diag = VECNORM(vec)

  IF (BC_halo_diag.LE.halo_eps) fullMesh = .TRUE.
  BC_halo_eps = halo_eps
END IF

#else
! get distance of diagonal of mesh
vec(1)   = GEO%xmaxglob-GEO%xminglob
vec(2)   = GEO%ymaxglob-GEO%yminglob
vec(3)   = GEO%zmaxglob-GEO%zminglob
BC_halo_eps = DOT_PRODUCT(vec,vec)
fullMesh = .TRUE.

firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

nBCSidesProc      = 0
offsetBCSides     = 0

! for fullMesh, each element requires ALL BC faces
IF (fullMesh) THEN
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

    DO iBCSide = 1,nUniqueBCSides
      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.ElemID) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO ! iBCSide

    ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
    IF (nBCSidesElem.GT.0) THEN
      ElemToBCSides(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
      ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides
    END IF

    offsetBCSides = nBCSidesProc
  END DO ! iElem

! .NOT. fullMesh
ELSE
  ! sum up all BC sides in range of BC_halo_eps
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)
      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

    ! loop over all sides. Check distance from every local side to total sides.
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)
      BCCNElemID = GetCNElemID(BCElemID)

      ! Ignore elements not on the compute node
      IF (BCCNElemID.EQ.-1) CYCLE

      ! Ignore the same element
      IF (BCElemID.EQ.ElemID) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,iElem) - ElemBaryNGeo(:,BCCNElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(iElem) + ElemRadiusNGeo(BCCNElemID))) CYCLE

      ! loop over all local sides of the element
      IF (VECNORM(ElemBaryNGeo(:,iElem) - BCSideMetrics(1:3,iBCSide)) &
        .LE. (BC_halo_eps + ElemRadiusNGeo(iElem) + BCSideMetrics(4,iBCSide))) THEN
           nBCSidesElem = nBCSidesElem + 1
           nBCSidesProc = nBCSidesProc + 1
      END IF
    END DO ! iBCSide

    ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
    IF (nBCSidesElem.GT.0) THEN
      ElemToBCSides(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
      ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides
    END IF

    offsetBCSides = nBCSidesProc
  END DO ! iElem
END IF ! fullMesh

! Find CN global number of BC sides and write Elem to BC Side mapping into shared array
#if USE_MPI
sendbuf = nBCSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetBCSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetBCSidesProc + nBCSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeBCSides = sendbuf

ElemToBCSides(ELEM_FIRST_BCSIDE,firstElem:lastElem) = ElemToBCSides(ELEM_FIRST_BCSIDE,firstElem:lastElem) + offsetBCSidesProc
#else
offsetBCSidesProc   = 0
nComputeNodeBCSides = nBCSidesProc
#endif /*USE_MPI*/

! Allocate shared array for BC sides
#if USE_MPI
CALL Allocate_Shared((/7,nComputeNodeBCSides/),SideBCMetrics_Shared_Win,SideBCMetrics_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideBCMetrics_Shared_Win,IERROR)
SideBCMetrics => SideBCMetrics_Shared
#else
ALLOCATE(SideBCMetrics(1:7,1:nComputeNodeBCSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  SideBCMetrics = -1.
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(ElemToBCSides_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideBCMetrics_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI*/

nBCSidesProc      = 0

! We did not know the number of BC sides before. Therefore, we need to do the check again and build the final mapping
! for fullMesh, each element requires ALL BC faces
IF (fullMesh) THEN
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(ElemID)
    END DO

    DO iBCSide = 1,nUniqueBCSides
      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.ElemID) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(ElemID)
    END DO ! iBCSide
  END DO ! iElem

! .NOT. fullMesh
ELSE
  ! sum up all BC sides in range of BC_halo_eps
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(ElemID)
    END DO

    ! loop over all sides. Check distance from every local side to total sides. Once a side has been flagged,
    ! it must not be counted again
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)
      BCCNElemID = GetCNElemID(BCElemID)

      ! Ignore elements not on the compute node
      IF (BCCNElemID.EQ.-1) CYCLE

      ! Ignore the same element
      IF (BCElemID.EQ.ElemID) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,iElem) - ElemBaryNGeo(:,BCCNElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(iElem) + ElemRadiusNGeo(BCCNElemID))) CYCLE

      ! loop over all local sides of the element
      IF (VECNORM(ElemBaryNGeo(:,iElem) - BCSideMetrics(1:3,iBCSide)) &
        .LE. (BC_halo_eps + ElemRadiusNGeo(iElem) + BCSideMetrics(4,iBCSide))) THEN
          nBCSidesProc = nBCSidesProc + 1
          SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
          SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(ElemID)
      END IF
    END DO ! iBCSide
  END DO ! iElem
END IF ! fullMesh

#if USE_MPI
CALL BARRIER_AND_SYNC(SideBCMetrics_Shared_Win,MPI_COMM_SHARED)

firstSide = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeBCSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeBCSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nComputeNodeBCSides
#endif /*USE_MPI*/

! calculate origin, radius and distance to sides
DO iSide = firstSide,lastSide
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,iSide))
  BCSideID = SideID2BCSide(SideID)
  CNElemID = GetCNElemID(INT(SideBCMetrics(BCSIDE_ELEMID,iSide)))

  !> get origin and radius from BC Side
  SideBCMetrics(5:7          ,iSide) = BCSideMetrics(1:3,BCSideID)
  SideBCMetrics(BCSIDE_RADIUS,iSide) = BCSideMetrics(4  ,BCSideID)

  !> build side distance
  origin(1:3) = ElemBaryNGeo(1:3,CNElemID)
  vec(1:3)    = origin(1:3) - BCSideMetrics(1:3,BCSideID)
  SideBCMetrics(BCSIDE_DISTANCE,iSide) = SQRT(DOT_PRODUCT(vec,vec))-ElemRadiusNGeo(CNElemID)-BCSideMetrics(4,BCSideID)
END DO ! iSide

#if USE_MPI
CALL BARRIER_AND_SYNC(SideBCMetrics_Shared_Win,MPI_COMM_SHARED)
#endif

! finally, sort by distance to help speed up BC tracking
!> allocate dummy array to hold variables
ALLOCATE(tmpSideBCMetrics(1:7,1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))
ALLOCATE(tmpSideBCDistance(   1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))
ALLOCATE(intSideBCMetrics(    1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))

DO iElem = firstElem,lastElem
  ! skip elements with no BC sides
  IF (ElemToBCSides(ELEM_NBR_BCSIDES,iElem).LE.0) CYCLE

  ! save values in temporary array
  firstSide    = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) + 1
  lastSide     = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) + ElemToBCSides(ELEM_NBR_BCSIDES,iElem)
  nBCSidesElem = ElemToBCSides(ELEM_NBR_BCSIDES,iElem)

  tmpSideBCMetrics(:,1:nBCSidesElem) = SideBCMetrics(:,firstSide:lastSide)
  tmpSideBCDistance( 1:nBCSidesElem) = SideBCMetrics(BCSIDE_DISTANCE,firstSide:lastSide)
  intSideBCMetrics(  1:nBCSidesElem) = (/((p),p=1,nBCSidesElem)/)

  ! sort SideID according to distance
  CALL InsertionSort(tmpSideBCDistance(1:nBCSidesElem),intSideBCMetrics(1:nBCSidesElem),nBCSidesElem)

  ! write back dummy array with variables
  DO iSide = 1,nBCSidesElem
    SideID = intSideBCMetrics(iSide)
    SideBCMetrics(:,firstSide+iSide-1) = tmpSideBCMetrics(:,SideID)
  END DO
END DO

DEALLOCATE(tmpSideBCMetrics)
DEALLOCATE(tmpSideBCDistance)
DEALLOCATE(intSideBCMetrics)

#if USE_MPI
CALL BARRIER_AND_SYNC(SideBCMetrics_Shared_Win,MPI_COMM_SHARED)
#endif

END SUBROUTINE BuildBCElemDistance


SUBROUTINE BuildNodeNeighbourhood()
!----------------------------------------------------------------------------------------------------------------------------------!
! Build shared arrays for mapping
! 1. UniqueNodeID -> all connected CN element IDs:    NodeToElemMapping(1,:) = OffsetNodeToElemMapping(:)
!                                                     NodeToElemMapping(2,:) = NbrOfElemsOnUniqueNode(:)
!
!                                                     NodeToElemInfo(nNodeToElemMapping)
!
! 2. CN elment ID -> all connected CN element IDs:    ElemToElemMapping(1,:) = Offset
!                                                     ElemToElemMapping(2,:) = Number of elements
!
!                                                     ElemToElemInfo(nElemToElemMapping)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals            ,ONLY: abort,UNIT_stdOUt,DisplayMessageAndTime!,myRank
USE MOD_Particle_Mesh_Vars ,ONLY: nUniqueGlobalNodes
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping,NodeToElemInfo,ElemToElemMapping,ElemToElemInfo
#if USE_MPI
USE MOD_Globals            ,ONLY: MPIRoot
USE MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping_Shared,NodeToElemInfo_Shared,ElemToElemMapping_Shared,ElemToElemInfo_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping_Shared_Win,NodeToElemInfo_Shared_Win
USE MOD_Particle_Mesh_Vars ,ONLY: ElemToElemMapping_Shared_Win,ElemToElemInfo_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE            :: NbrOfElemsOnUniqueNode(:),OffsetNodeToElemMapping(:)
INTEGER                        :: UniqueNodeID,NonUniqueNodeID,iElem,iNode,TestElemID,jElem,kElem
INTEGER                        :: OffsetCounter,OffsetElemToElemMapping,OffsetElemToElemCounter
INTEGER                        :: nNodeToElemMapping,iUniqueNode,firstElem,lastElem,nElemToElemMapping,CountElems
INTEGER,ALLOCATABLE            :: CheckedElemIDs(:)
#if USE_MPI
INTEGER                        :: sendbuf,recvbuf,iError
#endif /*USE_MPI*/
REAL                           :: StartT,EndT
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(A)',ADVANCE='NO') ' Building node neighbourhood ...'
GETTIME(StartT)

! 1.1 Get number of CN elements attached to each UNIQUE node and store in NbrOfElemsOnUniqueNode(UniqueNodeID)
! 1.2 Store the total number of counted elements in nNodeToElemMapping = SUM(NbrOfElemsOnUniqueNode)
! 1.3 Store the element offset for each UNIQUE node in OffsetNodeToElemMapping(iUniqueNode) by summing up the number of elements
!     from the first the previous (iUniqueNode-1) node

! Only the node leader fills the arrays
#if USE_MPI
IF(myComputeNodeRank.EQ.0)THEN
#endif /*USE_MPI*/
  ALLOCATE(NbrOfElemsOnUniqueNode(nUniqueGlobalNodes))
  NbrOfElemsOnUniqueNode=0
  ! Loop all CN elements (iElem is CNElemID)
#if USE_MPI
  DO iElem = 1,nComputeNodeTotalElems
#else
  DO iElem = 1,nElems
#endif /*USE_MPI*/
    ! Loop all local nodes
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
      UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
      NbrOfElemsOnUniqueNode(UniqueNodeID) = NbrOfElemsOnUniqueNode(UniqueNodeID) + 1
    END DO ! iNode = 1, 8
  END DO ! iElem = 1, nComputeNodeTotalElems
  nNodeToElemMapping = SUM(NbrOfElemsOnUniqueNode)

  ALLOCATE(OffsetNodeToElemMapping(nUniqueGlobalNodes))
  OffsetNodeToElemMapping = 0
  DO iUniqueNode = 2, nUniqueGlobalNodes
    OffsetNodeToElemMapping(iUniqueNode) = SUM(NbrOfElemsOnUniqueNode(1:iUniqueNode-1))
  END DO ! iUniqueNode = 1, nUniqueGlobalNodes
#if USE_MPI
END IF ! myComputeNodeRank.EQ.0
CALL MPI_BCAST(nNodeToElemMapping,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/


! 2. Allocate shared arrays for mapping
!    UniqueNodeID -> all CN element IDs to which it is connected      : NodeToElemMapping = [OffsetNodeToElemMapping(:),
!                                                                                            NbrOfElemsOnUniqueNode(:)]
!    NodeToElemMapping (offset and number of elements) -> CN element IDs : NodeToElemInfo = [CN elem IDs]
#if USE_MPI
! NodeToElemMapping
CALL Allocate_Shared((/2,nUniqueGlobalNodes/),NodeToElemMapping_Shared_Win,NodeToElemMapping_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeToElemMapping_Shared_Win,IERROR)
NodeToElemMapping => NodeToElemMapping_Shared
! NodeToElemInfo
CALL Allocate_Shared((/nNodeToElemMapping/),NodeToElemInfo_Shared_Win,NodeToElemInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeToElemInfo_Shared_Win,IERROR)
NodeToElemInfo => NodeToElemInfo_Shared
#else
ALLOCATE(NodeToElemMapping(2,nUniqueGlobalNodes))
ALLOCATE(NodeToElemInfo(nNodeToElemMapping))
#endif /*USE_MPI*/


! 3.1 Fill NodeToElemMapping = [OffsetNodeToElemMapping(:), NbrOfElemsOnUniqueNode(:)]
! 3.2 Store the CN element IDs in NodeToElemInfo(NodeToElemMapping(1,UniqueNodeID)+1 :
!                                                   NodeToElemMapping(1,UniqueNodeID)+NodeToElemMapping(2,UniqueNodeID))
! Now all CN elements attached to a UniqueNodeID can be accessed

! Only the node leader fills the arrays
#if USE_MPI
IF(myComputeNodeRank.EQ.0)THEN
#endif /*USE_MPI*/
  NodeToElemMapping = 0
  NodeToElemInfo = 0

  NodeToElemMapping(1,:) = OffsetNodeToElemMapping(:)
  DEALLOCATE(OffsetNodeToElemMapping)
  NodeToElemMapping(2,:) = NbrOfElemsOnUniqueNode(:)

  NbrOfElemsOnUniqueNode=0
  ! Loop all CN elements (iElem is CNElemID)
#if USE_MPI
  DO iElem = 1, nComputeNodeTotalElems
#else
  DO iElem = 1, nElems
#endif
    ! Loop all local nodes
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
      UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
      NbrOfElemsOnUniqueNode(UniqueNodeID) = NbrOfElemsOnUniqueNode(UniqueNodeID) + 1
      NodeToElemInfo(NodeToElemMapping(1,UniqueNodeID)+NbrOfElemsOnUniqueNode(UniqueNodeID)) = iElem
    END DO ! iNode = 1, 8
  END DO ! iElem = 1, nComputeNodeTotalElems
  DEALLOCATE(NbrOfElemsOnUniqueNode)
#if USE_MPI
END IF ! myComputeNodeRank.EQ.0

CALL BARRIER_AND_SYNC(NodeToElemInfo_Shared_Win   ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeToElemMapping_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

! 4. Allocate shared array for mapping
!    CN element ID -> all CN element IDs to which it is connected : ElemToElemMapping = [offset, Nbr of CN elements]
#if USE_MPI
! ElemToElemMapping
CALL Allocate_Shared((/2,nComputeNodeTotalElems/),ElemToElemMapping_Shared_Win,ElemToElemMapping_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToElemMapping_Shared_Win,IERROR)
ElemToElemMapping => ElemToElemMapping_Shared
IF (myComputeNodeRank.EQ.0) ElemToElemMapping = 0
CALL BARRIER_AND_SYNC(ElemToElemMapping_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(ElemToElemMapping(2,nElems))
ElemToElemMapping = 0
#endif /*USE_MPI*/


! 5. Fill ElemToElemMapping = [offset, Nbr of CN elements]
!    Note that the number of elements stored in ElemToElemMapping(2,iElem) must be shifted after communication with other procs
#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

OffsetCounter = 0
ALLOCATE(CheckedElemIDs(500))

! Loop all CN elements (iElem is CNElemID)
DO iElem = firstElem,lastElem
  CountElems = 0
  CheckedElemIDs = 0
  ! Loop all local nodes
  DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)

    ! Loop 1D array [offset + 1 : offset + NbrOfElems]
    ! (all CN elements that are connected to the local nodes)
    Elemloop: DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
      TestElemID = NodeToElemInfo(jElem)

      IF(iElem.EQ.TestElemID) CYCLE Elemloop

      ! Check previously stored element IDs and cycle if an already stored ID is encountered
      DO kElem = 1, CountElems
        IF(CheckedElemIDs(kElem).EQ.TestElemID) CYCLE Elemloop
      END DO ! kElem = 1, CountElems

      CountElems = CountElems + 1

      IF(CountElems.GT.500) CALL abort(&
      __STAMP__&
      ,'CountElems > 500. Inrease the number and try again!')

      CheckedElemIDs(CountElems) = TestElemID
      ! Note that the number of elements stored in ElemToElemMapping(2,iElem) must be shifted after communication with other procs
      ElemToElemMapping(2,iElem) = CountElems
    END DO Elemloop
  END DO ! iNode = 1, 8
  ElemToElemMapping(1,iElem) = OffsetCounter
  OffsetCounter = OffsetCounter + CountElems
END DO ! iElem = firstElem, lastElem

! 6. Find CN global number of connected CN elements (=nElemToElemMapping) and write into shared array
#if USE_MPI
sendbuf = OffsetCounter
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
OffsetElemToElemMapping   = recvbuf
! last proc knows CN total number of connected CN elements
sendbuf = OffsetElemToElemMapping + OffsetCounter
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nElemToElemMapping = sendbuf

!ElemToElemMapping(1,firstElem:lastElem) = ElemToElemMapping(1,firstElem:lastElem)
ElemToElemMapping(1,firstElem:lastElem) = ElemToElemMapping(1,firstElem:lastElem) + OffsetElemToElemMapping
#else
OffsetElemToElemMapping = 0
nElemToElemMapping = OffsetCounter
!ElemToElemMapping(:,firstElem:lastElem) = ElemToElemMapping(:,firstElem:lastElem)
#endif /*USE_MPI*/


! 7. Allocate shared array for mapping
!    CN element ID -> all CN element IDs to which it is connected : ElemToElemInfo = [CN elem IDs]
#if USE_MPI
! ElemToElemInfo
CALL Allocate_Shared((/nElemToElemMapping/),ElemToElemInfo_Shared_Win,ElemToElemInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToElemInfo_Shared_Win,IERROR)
ElemToElemInfo => ElemToElemInfo_Shared
#else
ALLOCATE(ElemToElemInfo(nElemToElemMapping))
#endif /*USE_MPI*/


! 8. Fill ElemToElemInfo = [CN elem IDs] by finding all nodes connected to a CN element
!    (and all subsequent CN elements that are connected to those nodes)
!    Store the CN element IDs in ElemToElemInfo(ElemToElemMapping(1,iElem)+1 :
!                                               ElemToElemMapping(1,iElem)+ElemToElemMapping(2,iElem))
OffsetElemToElemCounter = OffsetElemToElemMapping
! Loop all CN elements (iElem is CNElemID)
DO iElem = firstElem, lastElem
  CountElems = 0
  CheckedElemIDs = 0
  ! Loop all local nodes
  DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)

    ! Loop all CN elements that are connected to the local nodes
    Elemloop2: DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
      TestElemID = NodeToElemInfo(jElem)

      IF(iElem.EQ.TestElemID) CYCLE Elemloop2

      DO kElem = 1, CountElems
        IF(CheckedElemIDs(kElem).EQ.TestElemID) CYCLE Elemloop2
      END DO ! kElem = 1, CountElems

      CountElems = CountElems + 1
      OffsetElemToElemCounter = OffsetElemToElemCounter + 1

      IF(CountElems.GT.500) CALL abort(__STAMP__,'CountElems > 500. Inrease the number and try again!')

      CheckedElemIDs(CountElems) = TestElemID
      ElemToElemInfo(OffsetElemToElemCounter) = TestElemID

    END DO ElemLoop2
  END DO ! iNode = 1, 8

END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemToElemInfo_Shared_Win   ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemToElemMapping_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!')

END SUBROUTINE BuildNodeNeighbourhood


SUBROUTINE BuildElementOriginShared()
!================================================================================================================================
! compute the element origin at xi=(0,0,0)^T and set it as ElemBaryNGeo
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars          ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars ,ONLY: ElemBaryNGeo
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars, ONLY: ElemBaryNGeo_Shared,ElemBaryNGeo_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Mesh_Vars ,ONLY: nComputeNodeElems
USE MOD_Mesh_Vars          ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,i,j,k
REAL                           :: XPos(3),buf
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
!================================================================================================================================
#if USE_MPI
CALL Allocate_Shared((/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)
ElemBaryNGeo => ElemBaryNGeo_Shared

ASSOCIATE(XCL_NGeo => XCL_NGeo_Shared)

! Set ranges
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(ElemBaryNGeo(1:3,nComputeNodeElems))
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemBaryNGeo = 0.
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(ElemBaryNGeo_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI*/

! evaluate the polynomial at origin: Xi=(/0.0,0.0,0.0/)
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  xPos=0.

  DO k=0,NGeo
    DO j=0,NGeo
      buf=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        xPos=xPos+XCL_NGeo(1:3,i,j,k,ElemID)*Lag(1,i)*buf
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo

  ElemBaryNGeo(:,iElem)=xPos
END DO ! iElem

#if USE_MPI
END ASSOCIATE
CALL BARRIER_AND_SYNC(ElemBaryNGeo_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementOriginShared


SUBROUTINE BuildElementBasisAndRadius()
!================================================================================================================================
! Build the element local basis system, where the origin is located at xi=(0,0,0)^T and each local coordinate system is pointing
! to an element side
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis                  ,ONLY: DeCasteljauInterpolation
USE MOD_Basis                  ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars              ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadiusNGEO_Shared,ElemRadiusNGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo_Shared,ElemRadius2NGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: nELems
USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,SideID
INTEGER                        :: i,j,k,ilocSide
INTEGER                        :: iDir
REAL                           :: Xi(3,6),xPos(3),Radius
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
!================================================================================================================================
#if USE_MPI
  CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemRadiusNGeo_Shared_Win,ElemRadiusNGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadiusNGeo_Shared_Win,IERROR)
  CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemRadius2NGeo_Shared_Win,ElemRadius2NGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadius2NGeo_Shared_Win,IERROR)
  CALL Allocate_Shared((/3,6,nComputeNodeTotalElems/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
  CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
  CALL Allocate_Shared((/6,nComputeNodeTotalElems/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
  CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
  ElemRadiusNGeo     => ElemRadiusNGeo_Shared
  ElemRadius2NGeo    => ElemRadius2NGeo_Shared
  XiEtaZetaBasis     => XiEtaZetaBasis_Shared
  slenXiEtaZetaBasis => slenXiEtaZetaBasis_Shared

ASSOCIATE(XCL_NGeo     => XCL_NGeo_Shared)

#else
  ALLOCATE(ElemRadiusNGeo(          nElems) &
          ,ElemRadius2NGeo(         nElems) &
          ,XiEtaZetaBasis(1:3,1:6,1:nElems) &
          ,slenXiEtaZetaBasis(1:6,1:nElems))
#endif /*USE_MPI*/

#if USE_MPI
IF(myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ElemRadiusNGeo =0.
  ElemRadius2NGeo=0.
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(ElemRadiusNGeo_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemRadius2NGeo_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

#if USE_MPI
firstElem=INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif /*USE_MPI*/

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

! iElem is CN elem
DO iElem=firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos=0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos=xPos+XCL_NGeo(:,i,j,k,ElemID)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    XiEtaZetaBasis(1:3,iDir,iElem)=xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(:,iDir,iElem)=XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
  END DO ! iDir = 1, 6

  Radius=0.
  DO ilocSide=1,6
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)
    IF(SideID.EQ.-1) CYCLE
    DO j=0,NGeo
      DO i=0,NGeo
        xPos=BezierControlPoints3D(:,i,j,SideID)-ElemBaryNGeo(:,iElem)
        Radius=MAX(Radius,SQRT(DOT_PRODUCT(xPos,xPos)))
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO ! ilocSide
  ElemRadiusNGeo (iElem)=Radius
  ElemRadius2NGeo(iElem)=Radius*Radius
END DO ! iElem

#if USE_MPI
END ASSOCIATE
CALL BARRIER_AND_SYNC(ElemRadiusNGeo_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemRadius2NGeo_Shared_Win   ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(XiEtaZetaBasis_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(slenXiEtaZetaBasis_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementBasisAndRadius


SUBROUTINE BuildSideOriginAndRadius()
!===================================================================================================================================
! Globally identifies all BC sides and build side origin and radius
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                  ,ONLY: DeCasteljauInterpolation
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,nUniqueBCSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID_Shared,SideID2BCSide_Shared,BCSideMetrics_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID_Shared_Win,SideID2BCSide_Shared_Win,BCSideMetrics_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER                 :: xi(1:2) = 0
INTEGER                        :: p,q
INTEGER                        :: iSide,firstSide,lastSide,BCSideID
INTEGER                        :: nUniqueBCSidesProc,offsetUniqueBCSidesProc
REAL                           :: origin(1:3),radius,radiusMax,vec(1:3)
#if USE_MPI
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   )*REAL(nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! Count number of BC sides in range
nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
  ! ignore elements not on the compute node
  IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
END DO

! Find global number of BC sides and write side <=> BCSide mapping into shared array
#if USE_MPI
sendbuf = nUniqueBCSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetUniqueBCSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetUniqueBCSidesProc + nUniqueBCSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nUniqueBCSides = sendbuf

CALL Allocate_Shared((/nUniqueBCSides/),BCSide2SideID_Shared_Win,BCSide2SideID_Shared)
CALL MPI_WIN_LOCK_ALL(0,BCSide2SideID_Shared_Win,IERROR)
CALL Allocate_Shared((/nNonUniqueGlobalSides/),SideID2BCSide_Shared_Win,SideID2BCSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideID2BCSide_Shared_Win,IERROR)
BCSide2SideID => BCSide2SideID_Shared
SideID2BCSide => SideID2BCSide_Shared

! Also allocate array to hold BC Side metrics
CALL Allocate_Shared((/4,nUniqueBCSides/),BCSideMetrics_Shared_Win,BCSideMetrics_Shared)
CALL MPI_WIN_LOCK_ALL(0,BCSideMetrics_Shared_Win,IERROR)
BCSideMetrics => BCSideMetrics_Shared
#else
offsetUniqueBCSidesProc = 0
nUniqueBCSides = nUniqueBCSidesProc

ALLOCATE(BCSide2SideID(    1:nUniqueBCSides)        &
        ,SideID2BCSide(    1:nNonUniqueGlobalSides))
! Also allocate array to hold BC Side metrics
ALLOCATE(BCSideMetrics(1:4,1:nUniqueBCSides))
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  BCSide2SideID = -1
  SideID2BCSide = -1
  BCSideMetrics = -1
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(BCSide2SideID_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideID2BCSide_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BCSideMetrics_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
  ! ignore elements not on the compute node
  IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
  BCSideID           = offsetUniqueBCSidesProc + nUniqueBCSidesProc
  BCSide2SideID(BCSideID) = iSide
  SideID2BCSide(iSide)    = BCSideID

  ! calculate origin, radius for all BC sides
  !> build side origin
  ! TODO: BezierControlPoints are allocated with global side ID, so this SHOULD work. Breaks if we reduce the halo region
  CALL DeCasteljauInterpolation(NGeo,xi,iSide,origin)
  BCSideMetrics(1:3,BCSideID) = origin(1:3)

  !> build side radius
  radiusMax = 0.
  DO q = 0,NGeo
    DO p = 0,NGeo
      vec(1:3) = BezierControlPoints3D(:,p,q,iSide) - origin
      radius   = DOTPRODUCT(Vec)
      radiusMax= MAX(radiusMax,radius)
    END DO
  END DO
  BCSideMetrics(4,BCSideID) = SQRT(RadiusMax)
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(BCSide2SideID_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideID2BCSide_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BCSideMetrics_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE BuildSideOriginAndRadius


SUBROUTINE BuildLinearSideBaseVectors()
!===================================================================================================================================
! computes the face base vector for linear (planar or bilinear) face intersection calculation
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,BaseVectorsScale
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectorsScale_Shared,BaseVectorsScale_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectors0_Shared,BaseVectors1_Shared,BaseVectors2_Shared,BaseVectors3_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectors0_Shared_Win,BaseVectors1_Shared_Win,BaseVectors2_Shared_Win,BaseVectors3_Shared_Win
#endif /* USE_MPI */
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,firstSide,lastSide
REAL                           :: crossVec(3)
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)', ADVANCE='NO') ' GET LINEAR SIDE BASEVECTORS...'
#if USE_MPI
CALL Allocate_Shared((/3,nNonUniqueGlobalSides/),BaseVectors0_Shared_Win,BaseVectors0_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors0_Shared_Win,IERROR)
CALL Allocate_Shared((/3,nNonUniqueGlobalSides/),BaseVectors1_Shared_Win,BaseVectors1_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors1_Shared_Win,IERROR)
CALL Allocate_Shared((/3,nNonUniqueGlobalSides/),BaseVectors2_Shared_Win,BaseVectors2_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors2_Shared_Win,IERROR)
CALL Allocate_Shared((/3,nNonUniqueGlobalSides/),BaseVectors3_Shared_Win,BaseVectors3_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors3_Shared_Win,IERROR)
CALL Allocate_Shared((/nNonUniqueGlobalSides/),BaseVectorsScale_Shared_Win,BaseVectorsScale_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectorsScale_Shared_Win,IERROR)
BaseVectors0 => BaseVectors0_Shared
BaseVectors1 => BaseVectors1_Shared
BaseVectors2 => BaseVectors2_Shared
BaseVectors3 => BaseVectors3_Shared
BaseVectorsScale => BaseVectorsScale_Shared

firstSide = INT(REAL (myComputeNodeRank   )*REAL(nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE( BaseVectors0(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors1(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors2(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors3(1:3,1:nNonUniqueGlobalSides),&
          BaseVectorsScale(1:nNonUniqueGlobalSides))

firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,lastSide
    BaseVectors0(:,iSide) = (+BezierControlPoints3DElevated(:,0,0           ,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             +BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3DElevated(:,0,0           ,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             -BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3DElevated(:,0,0           ,iSide)-BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             +BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3DElevated(:,0,0           ,iSide)-BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             -BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
ELSE
  DO iSide=firstSide,lastSide
    BaseVectors0(:,iSide) = (+BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
END IF

#if USE_MPI
CALL BARRIER_AND_SYNC(BaseVectors0_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BaseVectors1_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BaseVectors2_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BaseVectors3_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BaseVectorsScale_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI */

LBWRITE(UNIT_stdOut,'(A)')' DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE BuildLinearSideBaseVectors


!===================================================================================================================================
!> Builds SideSlabNormals_Shared, SideSlabIntervals_Shared and  BoundingBoxIsEmpty_Shared from Bezier control points
!===================================================================================================================================
SUBROUTINE BuildSideSlabAndBoundingBox()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated,SideSlabNormals,SideSlabIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty,BezierElevation
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID,GetGlobalSideID
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /* USE_MPI */
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: firstSide,lastSide,iSide,SideID
#if !USE_MPI
INTEGER          :: ALLOCSTAT
#endif
#ifdef CODE_ANALYZE
! TODO
! REAL             :: dx,dy,dz
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

#if USE_MPI
CALL Allocate_Shared((/3,3,nComputeNodeTotalSides/),SideSlabNormals_Shared_Win,SideSlabNormals_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideSlabNormals_Shared_Win,IERROR)
CALL Allocate_Shared((/6,nComputeNodeTotalSides/),SideSlabIntervals_Shared_Win,SideSlabIntervals_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideSlabIntervals_Shared_Win,IERROR)
CALL Allocate_Shared((/nComputeNodeTotalSides/),BoundingBoxIsEmpty_Shared_Win,BoundingBoxIsEmpty_Shared)
CALL MPI_WIN_LOCK_ALL(0,BoundingBoxIsEmpty_Shared_Win,IERROR)
firstSide = INT(REAL (myComputeNodeRank   )*REAL(nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
SideSlabNormals    => SideSlabNormals_Shared
SideSlabIntervals  => SideSlabIntervals_Shared
BoundingBoxIsEmpty => BoundingBoxIsEmpty_Shared
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(SideSlabNormals(1:3,1:3,1:nNonUniqueGlobalSides) &
        ,SideSlabIntervals(  1:6,1:nNonUniqueGlobalSides) &
        ,BoundingBoxIsEmpty(     1:nNonUniqueGlobalSides) &
        ,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate SideMetrics arrays!')
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /* USE_MPI */
! TODO: bounding box volumes must be calculated for all unique sides.
!#ifdef CODE_ANALYZE
!    ALLOCATE(SideBoundingBoxVolume(nSides))
!#endif /*CODE_ANALYZE*/

IF (BezierElevation.GT.0) THEN
  DO iSide = firstSide,LastSide
    ! ignore sides that are not on the compute node
    ! IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

    SideID = GetGlobalSideID(iSide)

    ! Ignore small mortar sides attached to big mortar sides
    IF (SideInfo_Shared(SIDE_LOCALID,SideID).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,SideID).GT.6) CYCLE

    ! BezierControlPoints are always on nonUniqueGlobalSide
    CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,SideID) &
                                       ,SideSlabNormals(   1:3,1:3,iSide)                                       &
                                       ,SideSlabInterVals( 1:6    ,iSide)                                       &
                                       ,BoundingBoxIsEmpty(iSide))
  END DO
ELSE
  DO iSide=firstSide,LastSide
    ! ignore sides that are not on the compute node
    ! IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

    SideID = GetGlobalSideID(iSide)

    ! Ignore small mortar sides attached to big mortar sides
    IF (SideInfo_Shared(SIDE_LOCALID,SideID).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,SideID).GT.6) CYCLE

    ! BezierControlPoints are always on nonUniqueGlobalSide
    CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)                         &
                                       ,SideSlabNormals(   1:3,1:3,iSide)                                       &
                                       ,SideSlabInterVals( 1:6    ,iSide)                                       &
                                       ,BoundingBoxIsEmpty(iSide))
  END DO
END IF
#if USE_MPI
CALL BARRIER_AND_SYNC(SideSlabNormals_Shared_Win   ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideSlabIntervals_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(BoundingBoxIsEmpty_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI */
!#ifdef CODE_ANALYZE
! TODO: bounding box volumes must be calculated for all unique sides.
!               offsetSideID = ElemInfo_Shared(SideIf
!               DO iSide=offsetMPISides_YOUR,LastSide
!                 dx=ABS(SideSlabIntervals(2)-SideSlabIntervals(1))
!                 dy=ABS(SideSlabIntervals(4)-SideSlabIntervals(3))
!                 dz=ABS(SideSlabIntervals(6)-SideSlabIntervals(5))
!                 SideID = SideInfo
!                 SideBoundingBoxVolume(SideID)=dx*dy*dz
!               END DO
!#endif /*CODE_ANALYZE*/


END SUBROUTINE BuildSideSlabAndBoundingBox


END MODULE MOD_Particle_Mesh_Build
