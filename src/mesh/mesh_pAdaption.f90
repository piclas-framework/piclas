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

MODULE MOD_Mesh_pAdaption
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitpAdaption
#if !(PP_TimeDiscMethod==700)
PUBLIC::Set_N_DG_Mapping
#endif /*!(PP_TimeDiscMethod==700)*/

INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_DBG2 = -2 ! Debugging
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_DBG  = -1 ! Debugging
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_ZERO = 0  ! deactivate
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_RDN  = 1  ! random
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_NPB  = 2  ! Elements with non-periodic boundary sides get NMax
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_HH   = 3  ! Elements in the lower half domain in x-direction are set to NMin and the upper half are set to NMax. Origin must be at x=0.

INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_LVL_MINTWO  = -2 ! Directly connected elements are set to NMax and 2nd layer elements are set to NMin+1
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_LVL_MINONE  = -1 ! Directly connected elements are set to NMin+1
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_LVL_DEFAULT = 1 ! Directly connected elements are set to NMax
INTEGER,PARAMETER,PUBLIC :: PRM_P_ADAPTION_LVL_TWO     = 2 ! Directly connected elements and elements connected with these are set to NMax
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get local N for each element and side
!===================================================================================================================================
SUBROUTINE InitpAdaption()
! MODULES
USE MOD_Globals
USE MOD_PreProc
!USE MOD_DG_Vars            ,ONLY: DG_Elems_master,DG_Elems_slave
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars            ,ONLY: N_DG,pAdaptionType,pAdaptionBCLevel,NDGAllocationIsDone
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_IO_HDF5            ,ONLY: AddToElemData,ElementOut,ElementOutNloc
USE MOD_Mesh_Vars          ,ONLY: nElems,SideToElem,nBCSides,Boundarytype,BC,readFEMconnectivity,NodeCoords,nGlobalDOFs
USE MOD_Mesh_Vars          ,ONLY: NMaxGlobal,NMinGlobal
USE MOD_ReadInTools        ,ONLY: GETINTFROMSTR
USE MOD_Interpolation_Vars ,ONLY: NMax,NMin
USE MOD_Mesh_Vars          ,ONLY: readFEMconnectivity, offsetElem, nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,BCSideID,BCType
REAL    :: RandVal,x
LOGICAL :: SetBCElemsToNMax
INTEGER(KIND=8) :: nLocalDOFs
!===================================================================================================================================
! Set defaults
SetBCElemsToNMax = .FALSE. ! Initialize

#if !(PP_TimeDiscMethod==700)
pAdaptionBCLevel = -1
NDGAllocationIsDone = .FALSE.

! Read p-adaption specific input data
pAdaptionType = GETINTFROMSTR('pAdaptionType')
#if USE_HDG && !(USE_PETSC)
! Sanity check: p-adaption without petsc is not allowed
IF (pAdaptionType.NE.PRM_P_ADAPTION_ZERO) CALL abort(__STAMP__,'p-adaption is only implemented with petsc. Set LIBS_USE_PETSC=ON')
#endif /*USE_HDG && !(USE_PETSC)*/

! Allocate arrays and initialize local polynomial degree
! This happens here because nElems is determined here and N_DG is required below for the mesh initialisation
ALLOCATE(N_DG(1:nElems))
N_DG = -1
!CALL AddToElemData(ElementOut,'Nloc',IntArray=N_DG_Mapping(2,1+offSetElem:nElems+offSetElem)) ! TODO: Why does this not work?
! Add array containing the local polynomial degree to the hdf5 output
CALL AddToElemData(ElementOut,'Nloc',IntArray=N_DG)
! Add the same container to a separate object for writing it to special .h5 files that contain other fields using Nloc
CALL AddToElemData(ElementOutNloc,'Nloc',IntArray=N_DG)

SELECT CASE(pAdaptionType)
CASE(PRM_P_ADAPTION_DBG2) ! Debugging
  ! Nloc = 2,3,4,5
  DO iElem=1,nElems
    N_DG(iElem) = 2
    IF (iElem+offsetElem.eq.1) THEN
      N_DG(iElem) = 4
    ELSEIF (iElem+offsetElem.eq.2) THEN
      N_DG(iElem) = 5
    ELSEIF (iElem+offsetElem.eq.4) THEN
      N_DG(iElem) = 3
    END IF ! iElem+offset
  END DO
CASE(PRM_P_ADAPTION_DBG) ! Debugging
  ! Nloc = 2 for all elements, except for iGlogalElemID=1,6,10
  DO iElem=1,nElems
    N_DG(iElem) = 2
    IF (iElem+offsetElem.eq.1) THEN
      N_DG(iElem) = 1
    ELSEIF (iElem+offsetElem.eq.6) THEN
      N_DG(iElem) = 1
    ELSEIF (iElem+offsetElem.eq.10) THEN
      N_DG(iElem) = 1
    END IF ! iElem+offset
  END DO
CASE(PRM_P_ADAPTION_ZERO)
  N_DG = PP_N ! By default, the initial degree is set to PP_N
CASE(PRM_P_ADAPTION_RDN) ! Random between NMin and NMax
  DO iElem=1,nElems
    CALL RANDOM_NUMBER(RandVal)
    N_DG(iElem) = NMin + INT(RandVal*(NMax-NMin+1))
  END DO
CASE(PRM_P_ADAPTION_NPB) ! Non-periodic BCs are set to NMax
  ! Get depth of increased polynomial degree
  pAdaptionBCLevel = GETINTFROMSTR('pAdaptionBCLevel')
  N_DG = Nmin ! By default, the initial degree is set to Nmin
  SetBCElemsToNMax = .TRUE.
CASE(PRM_P_ADAPTION_HH) ! Elements in the lower half domain in x-direction are set to NMin and the upper half are set to NMax
  DO iElem=1,nElems
    x = (MAXVAL(NodeCoords(1,:,:,:,iElem))+MINVAL(NodeCoords(1,:,:,:,iElem)))/2.0
    IF(x.GT.0.0)THEN
      N_DG(iElem) = NMax
    ELSE
      N_DG(iElem) = NMin
    END IF ! x.GT.0.0
  END DO
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown pAdaptionType!' ,IntInfo=pAdaptionType)
END SELECT

! Check if BC elements are to be set to a higher polynomial degree
IF(ABS(pAdaptionBCLevel).GT.1.OR.&
  (ABS(pAdaptionBCLevel).GT.0.AND.readFEMconnectivity))THEN
  CALL SetpAdaptionBCLevel()
ELSE
  ! Check if all BC elements are set to Nmax
  IF(SetBCElemsToNMax)THEN
    DO BCSideID=1,nBCSides
      BCType=Boundarytype(BC(BCSideID),BC_TYPE)
      IF(BCType.EQ.1) CYCLE ! Skip periodic sides
      iElem       = SideToElem(S2E_ELEM_ID,BCSideID)
      N_DG(iElem) = NMax
    END DO ! BCSideID=1,nBCSides
  END IF ! SetBCElemsToNMax
END IF

! Sanity check and nGlobalDOFs calculation
nLocalDOFs=0
DO iElem=1,nElems
  IF((N_DG(iElem).LT.NMin).OR.(N_DG(iElem).GT.NMax))THEN
    IPWRITE(*,*) "iElem       = ", iElem
    IPWRITE(*,*) "N_DG(iElem) = ", N_DG(iElem)
    IPWRITE(*,*) "NMin        = ", NMin
    IPWRITE(*,*) "NMax        = ", NMax
  END IF
  IF(N_DG(iElem).LT.NMin) CALL abort(__STAMP__,'N_DG(iElem)<NMin')
  IF(N_DG(iElem).GT.NMax) CALL abort(__STAMP__,'N_DG(iElem)>NMax')
  nLocalDOFs = nLocalDOFs + INT((N_DG(iElem)+1)**3,8)
END DO

NMinGlobal = MINVAL(N_DG)
NMaxGlobal = MAXVAL(N_DG)
#if USE_MPI
! Get global min/max polynomial degree that are actually present (not the theoretical limits)
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE , NMinGlobal , 1 , MPI_INTEGER , MPI_MIN , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(MPI_IN_PLACE , NMaxGlobal , 1 , MPI_INTEGER , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
ELSE
  CALL MPI_REDUCE(NMinGlobal   , 0          , 1 , MPI_INTEGER , MPI_MIN , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(NMaxGlobal   , 0          , 1 , MPI_INTEGER , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ! in this case the receive value is not relevant.
END IF
! Calculate the sum across all processes. Only the MPI root process needs this information
IF(MPIRoot)THEN
  CALL MPI_REDUCE(nLocalDOFs , nGlobalDOFs , 1 , MPI_INTEGER8 , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  ! Sanity check
  IF(nGlobalDOFs.LE.0) CALL abort(__STAMP__,'nGlobalDOFs.LE.0')
ELSE
  CALL MPI_REDUCE(nLocalDOFs , 0           , 1 , MPI_INTEGER8 , MPI_SUM , 0 , MPI_COMM_PICLAS , IError)
END IF ! MPIRoot
#else
nGlobalDOFs = nLocalDOFs
#endif /*USE_MPI*/

! Initialize element containers
CALL Build_N_DG_Mapping()
#else
NMinGlobal = PP_N
NMaxGlobal = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/

END SUBROUTINE InitpAdaption


#if !(PP_TimeDiscMethod==700)
!===================================================================================================================================
!> Set Nloc of BC elements to a higher polynomial degree than in the volume
!>
!> 1. Loop over the processor-local elements
!> 2. Loop over the corner vertices of the element
!> 3. Use VertexConnectInfo to get the neighbour element index and vertex for all possible connection (also periodic)
!> 4. Check the sides of connected to the neighbour node and find out if the side is a BC side
!> 5. For further layers, loop over the elements again and check for already marked elements instead of the sides
!===================================================================================================================================
SUBROUTINE SetpAdaptionBCLevel()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: readFEMconnectivity, offsetElem, nElems
USE MOD_Mesh_Vars          ,ONLY: ElemInfo,VertexInfo,VertexConnectInfo
USE MOD_Mesh_Vars          ,ONLY: BoundaryType
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_DG_Vars            ,ONLY: N_DG,pAdaptionBCLevel,N_DG_Mapping
USE MOD_Interpolation_Vars ,ONLY: NMax,NMin
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping_Shared_Win
USE MOD_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,BCType,NonUniqueGlobalSideID,iGlobalElemID,BCIndex,ElemType,OffsetCounter
INTEGER :: iVertexConnect,GlobalNbElemID,GlobalNbLocVertexID,LocSideList(3),iLocSideList,iLocSide
INTEGER :: FirstElemInd,LastElemInd
INTEGER :: FirstVertexInd,LastVertexInd,FirstVertexConnectInd,LastVertexConnectInd
!===================================================================================================================================
! Do not re-allocate during load balance here as it is communicated between the processors
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

! Sanity check: This routine requires FEM connectivity
IF(.NOT.readFEMconnectivity) CALL abort(__STAMP__,'Error in p-adaption init: readFEMconnectivity=T is required!')

! Element index
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems

! Loop over the process-local global elements indices
GlobalElemIDLoop: DO iGlobalElemID = FirstElemInd, LastElemInd
  iElem = iGlobalElemID - offsetElem
  ElemType = ElemInfo(ELEM_TYPE,iGlobalElemID)
  ! Sanity check: currently only hexahedral elements are implemented
  SELECT CASE(ElemType)
  CASE(108,118,208)
    ! Hexahedral elements
  CASE DEFAULT
    CALL abort(__STAMP__,'Element type not implemented: ElemType =',IntInfoOpt=ElemType)
  END SELECT
  ! Get local VertexInfo of current element
  FirstVertexInd = ElemInfo(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1
  LastVertexInd  = ElemInfo(ELEM_LASTVERTEXIND,iGlobalElemID)
  ! Get local vertex connectivity
  FirstVertexConnectInd = VertexInfo(VERTEX_FIRSTCONNECTIND,FirstVertexInd)+1
  LastVertexConnectInd  = VertexInfo(VERTEX_LASTCONNECTIND,LastVertexInd)
  VertexConnectLoop: DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
    ! Check if current element has already been flagged
    IF(N_DG(iElem).EQ.NMax) EXIT VertexConnectLoop
    ! Get neighbour infos
    GlobalNbElemID      = ABS(VertexConnectInfo(VERTEXCONNECT_NBELEMID   ,iVertexConnect))
    GlobalNbLocVertexID = VertexConnectInfo(VERTEXCONNECT_NBLOCNODEID,iVertexConnect)
    ! Set sides depending on the element type: Only implemented for Hexahedral elements
    CALL GetLocSideList(ElemType,GlobalNbLocVertexID,LocSideList)
    LocSideListLoop: DO iLocSideList = 1, 3
      ! Check if current element has already been flagged
      IF(N_DG(iElem).EQ.NMax) EXIT VertexConnectLoop
      iLocSide = LocSideList(iLocSideList)
      NonUniqueGlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalNbElemID) + iLocSide
      BCIndex = SideInfo_Shared(SIDE_BCID,NonUniqueGlobalSideID)
      IF(BCIndex.LE.0) CYCLE LocSideListLoop ! Skip inner sides
      BCType = BoundaryType(BCIndex,BC_TYPE)
      IF(BCType.LE.1) CYCLE LocSideListLoop ! Skip periodic sides
      IF(BCType.EQ.10) CYCLE LocSideListLoop ! Skip Neumann sides
      IF(pAdaptionBCLevel.EQ.-1)THEN
        N_DG(iElem) = NMin+1
      ELSE
        N_DG(iElem) = NMax
      END IF ! pAdaptionBCLevel.EQ.-1
    END DO LocSideListLoop ! iLocSideList = 1, 3
  END DO VertexConnectLoop ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
END DO GlobalElemIDLoop ! iGlobalElemID = FirstElemInd, LastElemInd


! For further layers, loop over the elements again and check for already marked elements instead of the sides
IF(ABS(pAdaptionBCLevel).GT.1)THEN
  ! Allocate the shared memory container and associate pointer: N_DG_Mapping
  CALL Allocate_N_DG_Mapping()

  ! Set Nloc in N_DG_Mapping
  CALL Set_N_DG_Mapping(OffsetCounter)

#if USE_MPI
  ! Synchronize shared array before utilization as GlobalNbElemID can be outside of the processors local elements
  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

  ! Loop over the process-local global elements indices
  iGlobalElemID_loop: DO iGlobalElemID = FirstElemInd, LastElemInd
    iElem = iGlobalElemID - offsetElem
    IF(N_DG(iElem).GT.NMin) CYCLE iGlobalElemID_loop
    ! Get local VertexInfo of current element
    FirstVertexInd = ElemInfo(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1
    LastVertexInd  = ElemInfo(ELEM_LASTVERTEXIND,iGlobalElemID)
    ! Get local vertex connectivity
    FirstVertexConnectInd = VertexInfo(VERTEX_FIRSTCONNECTIND,FirstVertexInd)+1
    LastVertexConnectInd  = VertexInfo(VERTEX_LASTCONNECTIND,LastVertexInd)
    iVertexConnect_loop: DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
      ! Check if current element has already been flagged
      IF(N_DG(iElem).GT.NMin) EXIT iVertexConnect_loop
      ! Get neighbour infos
      GlobalNbElemID = ABS(VertexConnectInfo(VERTEXCONNECT_NBELEMID   ,iVertexConnect))
      IF(N_DG_Mapping(2,GlobalNbElemID).EQ.NMax)THEN
        IF(pAdaptionBCLevel.EQ.-2)THEN
          N_DG(iElem) = NMin+1
        ELSE
          N_DG(iElem) = NMax
        END IF ! pAdaptionBCLevel.EQ.-2
      END IF ! N_DG_Mapping(2,GlobalNbElemID).EQ.NMax
    END DO iVertexConnect_loop ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
  END DO iGlobalElemID_loop ! iGlobalElemID = FirstElemInd, LastElemInd

#if USE_MPI
  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/
END IF ! pAdaptionBCLevel.GT.1
END SUBROUTINE SetpAdaptionBCLevel


!===================================================================================================================================
!> Allocate the shared memory container N_DG_Mapping_Shared
!===================================================================================================================================
SUBROUTINE Allocate_N_DG_Mapping()
! MODULES
USE MOD_Globals
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping,NDGAllocationIsDone
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_Mesh_Vars       ,ONLY: nGlobalElems
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping_Shared,N_DG_Mapping_Shared_Win
USE MOD_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED, myComputeNodeRank
#else
USE MOD_Mesh_Vars       ,ONLY: nElems
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
CALL Allocate_Shared((/3,nGlobalElems/),N_DG_Mapping_Shared_Win,N_DG_Mapping_Shared)
CALL MPI_WIN_LOCK_ALL(0,N_DG_Mapping_Shared_Win,IERROR)
N_DG_Mapping => N_DG_Mapping_Shared
IF (myComputeNodeRank.EQ.0) N_DG_Mapping = 0
CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(N_DG_Mapping(3,nElems))
N_DG_Mapping = 0
#endif /*USE_MPI*/
NDGAllocationIsDone = .TRUE.
END SUBROUTINE Allocate_N_DG_Mapping


!===================================================================================================================================
!> Loop over all local elements and set N_DG_Mapping(1:2,iElem+offSetElem).
!> N_DG_Mapping(1,:) can only be utilized after the barrier & sync, and the communication between the processes in Build_N_DG_Mapping
!> N_DG_Mapping(2,:) requires at least a barrier & sync, or each process can only access its own local elements
!===================================================================================================================================
SUBROUTINE Set_N_DG_Mapping(OffsetCounter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY: nElems,offSetElem
USE MOD_DG_Vars   ,ONLY: N_DG_Mapping,N_DG
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: OffsetCounter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: locDofs,iElem,iGlobalElem,Nloc
!===================================================================================================================================
OffsetCounter = 0
! Loop over local elements using the global ID
DO iGlobalElem = 1+offSetElem,nElems+offSetElem
  iElem = iGlobalElem-offSetElem
  Nloc = N_DG(iElem)
  locDofs = (Nloc+1)**3
  N_DG_Mapping(2,iGlobalElem) = Nloc
  N_DG_Mapping(1,iGlobalElem) = OffsetCounter
  OffsetCounter = OffsetCounter + locDofs
END DO ! iGlobalElem
END SUBROUTINE Set_N_DG_Mapping


!===================================================================================================================================
!> Create shared memory array N_DG_Mapping containing the global element information
!>   N_DG_Mapping(1,nElems+offSetElem): DOF offset
!>   N_DG_Mapping(2,nElems+offSetElem): element polynomial degree Nloc
!===================================================================================================================================
SUBROUTINE Build_N_DG_Mapping()
! MODULES
USE MOD_Globals
USE MOD_DG_Vars            ,ONLY: nDofsMapping
USE MOD_DG_Vars            ,ONLY: NDGAllocationIsDone
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem,nGlobalElems
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping_Shared_Win,displsDofs,recvcountDofs,N_DG_Mapping
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_LEADERS_SHARED, MPI_COMM_SHARED, myComputeNodeRank, myleadergrouprank
USE MOD_MPI_Shared_Vars    ,ONLY: nLeaderGroupProcs,nComputeNodeProcessors
USE MOD_Particle_Mesh_Vars ,ONLY: offsetComputeNodeElem
USE MOD_MPI_Shared
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: OffsetCounter,OffsetN_DG_Mapping
#if USE_MPI
INTEGER             :: iProc
INTEGER             :: sendbuf,recvbuf
#endif
!===================================================================================================================================
! Do not re-allocate during load balance here as it is communicated between the processors
#if USE_LOADBALANCE
IF(PerformLoadBalance)THEN
#endif /*USE_LOADBALANCE*/

  ! N_DG_Mapping is already set
  !OffsetCounter = N_DG_Mapping(1,nElems+offSetElem) + (N_DG_Mapping(2,nElems+offSetElem)+1)**3

#if USE_LOADBALANCE
ELSE
#endif /*USE_LOADBALANCE*/
  IF(.NOT.NDGAllocationIsDone)THEN
    ! Allocate the shared memory container and associate pointer: N_DG_Mapping
    CALL Allocate_N_DG_Mapping()
  END IF ! .NOT.NDGAllocationIsDone

  ! Set Nloc in N_DG_Mapping
  CALL Set_N_DG_Mapping(OffsetCounter)

#if USE_MPI
  ! Set the correct offsets in N_DG_Mapping(1,:) for parallel simulation
  sendbuf = OffsetCounter
  recvbuf = 0
  CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
  OffsetN_DG_Mapping   = recvbuf
  ! The last process knows CN total number of connected CN elements
  sendbuf = OffsetN_DG_Mapping + OffsetCounter
  CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nProcessors-1,MPI_COMM_PICLAS,iError)
  nDofsMapping = sendbuf

  N_DG_Mapping(1,1+offSetElem:nElems+offSetElem) = N_DG_Mapping(1,1+offSetElem:nElems+offSetElem) + OffsetN_DG_Mapping
  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)

  ! Communication between nodes
  IF (nComputeNodeProcessors.NE.nProcessors.AND.myComputeNodeRank.EQ.0) THEN
    ! Arrays for the compute node to hold the elem offsets
    ALLOCATE(displsDofs(   0:nLeaderGroupProcs-1), recvcountDofs(0:nLeaderGroupProcs-1))
    displsDofs(myLeaderGroupRank) = offsetComputeNodeElem
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsDofs,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountDofs(iProc-1) = displsDofs(iProc)-displsDofs(iProc-1)
    END DO
    recvcountDofs(nLeaderGroupProcs-1) = nGlobalElems - displsDofs(nLeaderGroupProcs-1)

    CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
        , 0                             &
        , MPI_DATATYPE_NULL             &
        , N_DG_Mapping               &
        , 3*recvcountDofs   &
        , 3*displsDofs      &
        , MPI_INTEGER          &
        , MPI_COMM_LEADERS_SHARED       &
        , IERROR)

    displsDofs(myLeaderGroupRank) = N_DG_Mapping(1,1+offSetElem)
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsDofs,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountDofs(iProc-1) = displsDofs(iProc)-displsDofs(iProc-1)
    END DO
    recvcountDofs(nLeaderGroupProcs-1) = nDofsMapping - displsDofs(nLeaderGroupProcs-1)
    END IF

  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win ,MPI_COMM_SHARED)

#else
  OffsetN_DG_Mapping = 0
  nDofsMapping = OffsetCounter
#endif /*USE_MPI*/

#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/
END SUBROUTINE Build_N_DG_Mapping


!===================================================================================================================================
!> Returns a list of sides (depending on the CGNS ordering) for a given local corner node index
!===================================================================================================================================
SUBROUTINE GetLocSideList(ElemType,iLocNode,LocSideList)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)  :: ElemType
INTEGER, INTENT(IN)  :: iLocNode
INTEGER, INTENT(OUT) :: LocSideList(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Set sides depending on the element type: Only implemented for Hexahedral elements
SELECT CASE(ElemType)
CASE(108,118,208)
  ! Hexahedral elements
  SELECT CASE(iLocNode)
  CASE(1)
    LocSideList=(/1,2,5/)
  CASE(2)
    LocSideList=(/1,2,3/)
  CASE(3)
    LocSideList=(/1,3,4/)
  CASE(4)
    LocSideList=(/1,4,5/)
  CASE(5)
    LocSideList=(/2,5,6/)
  CASE(6)
    LocSideList=(/2,3,6/)
  CASE(7)
    LocSideList=(/3,4,6/)
  CASE(8)
    LocSideList=(/4,5,6/)
  CASE DEFAULT
    CALL abort(__STAMP__,'Wrong iLocNode',IntInfoOpt=iLocNode)
  END SELECT
CASE DEFAULT
  CALL abort(__STAMP__,'Element type not implemented: ElemType =',IntInfoOpt=ElemType)
END SELECT
END SUBROUTINE GetLocSideList
#endif /*!(PP_TimeDiscMethod==700)*/

END MODULE MOD_Mesh_pAdaption