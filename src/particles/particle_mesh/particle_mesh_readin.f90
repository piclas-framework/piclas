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

!===================================================================================================================================
! Builds the mesh for particle tracking, separate from the DG mesh
!===================================================================================================================================
MODULE MOD_Particle_Mesh_Readin
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables
!-----------------------------------------------------------------------------------------------------------------------------------
! Interfaces
INTERFACE ReadMeshBasics
  MODULE PROCEDURE ReadMeshBasics
END INTERFACE

INTERFACE ReadMeshSideNeighbors
  MODULE PROCEDURE ReadMeshSideNeighbors
END INTERFACE

INTERFACE StartCommunicateMeshReadin
  MODULE PROCEDURE StartCommunicateMeshReadin
END INTERFACE

INTERFACE FinishCommunicateMeshReadin
  MODULE PROCEDURE FinishCommunicateMeshReadin
END INTERFACE

PUBLIC :: ReadMeshBasics
PUBLIC :: ReadMeshSideNeighbors
PUBLIC :: StartCommunicateMeshReadin
PUBLIC :: FinishCommunicateMeshReadin
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadMeshBasics()
!===================================================================================================================================
! Read basic global counters from mesh file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input                ,ONLY: File_ID,ReadAttribute
USE MOD_Particle_Mesh_Vars        ,ONLY: nNonUniqueGlobalSides,nNonUniqueGlobalNodes
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

!CALL ReadAttribute(File_ID,'nUniqueSides',1,IntScalar=nGlobalUniqueSidesFromMesh)
CALL ReadAttribute(File_ID,'nSides'      ,1,IntScalar=nNonUniqueGlobalSides)
CALL ReadAttribute(File_ID,'nNodes'      ,1,IntScalar=nNonUniqueGlobalNodes)

END SUBROUTINE ReadMeshBasics


SUBROUTINE ReadMeshSideNeighbors(ElemID,SideID)
!===================================================================================================================================
! Fills temporary array to add side neighbors to SideInfo(_Shared)
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
USE MOD_Particle_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
INTEGER,INTENT(IN)             :: SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (ElemID.EQ.0) THEN
  ! no connection
  SideInfo_Shared_tmp(SideID) = 0
ELSE
#if USE_MPI
  IF (.NOT.ElementOnNode(ElemID)) THEN
    ! neighbour element is outside of compute-node
    SideInfo_Shared_tmp(SideID) = 2
  ELSE
#endif /*USE_MPI*/
    SideInfo_Shared_tmp(SideID) = 1
#if USE_MPI
  END IF
#endif /*USE_MPI*/
END IF

END SUBROUTINE ReadMeshSideNeighbors


SUBROUTINE StartCommunicateMeshReadin()
!===================================================================================================================================
! Communicates the readin mesh between MPI leaders
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: StartT
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: nSideIDs,offsetSideID
#if USE_MPI
INTEGER                        :: iProc
INTEGER                        :: offsetNodeID!,nNodeIDs
#endif /*USE_MPI*/
!===================================================================================================================================

! Start timer: finished in FinishCommunicateMeshReadin()
GETTIME(StartT)

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
! nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
#else
FirstElemInd = 1
LastElemInd  = nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  IF (myComputeNodeRank.EQ.0) THEN
    LBWRITE(UNIT_stdOut,'(A)',ADVANCE="NO") ' Updating mesh on shared memory...'

    ! Arrays for the compute node to hold the elem offsets
    ALLOCATE(displsElem(   0:nLeaderGroupProcs-1), recvcountElem(0:nLeaderGroupProcs-1))
    displsElem(myLeaderGroupRank) = offsetComputeNodeElem
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
    END DO
    recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
  END IF

  ! Broadcast compute node side offset on node
  offsetComputeNodeSide=offsetSideID
  CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

  IF (myComputeNodeRank.EQ.0) THEN
    ! Arrays for the compute node to hold the side offsets
    ALLOCATE(displsSide(   0:nLeaderGroupProcs-1), recvcountSide(0:nLeaderGroupProcs-1))
    displsSide(myLeaderGroupRank) = offsetComputeNodeSide
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
    END DO
    recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)

    ! Gather mesh information in a non-blocking way
    ALLOCATE(MPI_COMM_LEADERS_REQUEST(1))
    ! ElemInfo_Shared only needs ELEM_RANK updated, performed in loadbalance_tools.f90
    ! CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE       *recvcountElem  &
    !     ,ELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
    CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)   *recvcountSide  &
        ,(SIDEINFOSIZE+1)*displsSide ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)

  END IF

  ! Broadcast compute node node offset on node
  offsetComputeNodeNode=offsetNodeID
  CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

  RETURN
END IF
#endif /*USE_LOADBALANCE*/

LBWRITE(UNIT_stdOut,'(A)',ADVANCE="NO") ' Communicating mesh on shared memory...'

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the elem offsets
  ALLOCATE(displsElem(   0:nLeaderGroupProcs-1), recvcountElem(0:nLeaderGroupProcs-1))
  displsElem(myLeaderGroupRank) = offsetComputeNodeElem
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
  END DO
  recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
END IF

! Broadcast compute node side offset on node
offsetComputeNodeSide=offsetSideID
CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the side offsets
  ALLOCATE(displsSide(   0:nLeaderGroupProcs-1), recvcountSide(0:nLeaderGroupProcs-1))
  displsSide(myLeaderGroupRank) = offsetComputeNodeSide
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
  END DO
  recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)
END IF

! Broadcast compute node node offset on node
offsetComputeNodeNode=offsetNodeID
CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the node offsets
  ALLOCATE(displsNode(   0:nLeaderGroupProcs-1), recvcountNode(0:nLeaderGroupProcs-1))
  displsNode(myLeaderGroupRank) = offsetComputeNodeNode
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsNode,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountNode(iProc-1) = displsNode(iProc)-displsNode(iProc-1)
  END DO
  recvcountNode(nLeaderGroupProcs-1) = nNonUniqueGlobalNodes - displsNode(nLeaderGroupProcs-1)

  ! Gather mesh information in a non-blocking way
  ALLOCATE(MPI_COMM_LEADERS_REQUEST(1:4))
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE    *recvcountElem  &
      ,ELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)*recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(2),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeInfo_Shared,                 recvcountNode  &
      ,displsNode                  ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(3),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeCoords_Shared,3             *recvcountNode  &
      ,3*displsNode                ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(4),IERROR)
END IF
#endif /*USE_MPI*/

END SUBROUTINE StartCommunicateMeshReadin


SUBROUTINE FinishCommunicateMeshReadin()
!===================================================================================================================================
! Communicates the readin mesh between MPI leaders
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: CommMeshReadinWallTime,StartT
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: FirstElemInd,LastElemInd
INTEGER :: iElem,NbElemID
INTEGER :: nSideIDs,offsetSideID
INTEGER :: iSide,sideCount
INTEGER :: iLocSide,jLocSide,nlocSides,nlocSidesNb,NbSideID
REAL    :: EndT
!===================================================================================================================================

#if USE_MPI
! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#else
FirstElemInd = 1
LastElemInd  = nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  ! Finish non-blocking mesh communication
  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_WAITALL(1,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,IERROR)
    DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
  END IF

  ! final sync of all mesh shared arrays
  CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  ! Write compute-node local SIDE_NBELEMTYPE
  IF (myComputeNodeRank.EQ.0) THEN
    SideInfo_Shared(SIDE_NBELEMTYPE,:) = 0
  END IF
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  SideInfo_Shared(SIDE_NBELEMTYPE,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp
  DEALLOCATE(SideInfo_Shared_tmp)
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  RETURN
END IF
#endif  /*USE_LOADBALANCE*/

#if USE_MPI
! Finish non-blocking mesh communication
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_WAITALL(4,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,IERROR)
  DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
END IF

! Ensure communication for determination of SIDE_LOCALID
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

DO iElem = FirstElemInd,LastElemInd
  iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  SideInfo_Shared(SIDE_ELEMID,iSide+1:ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)) = iElem
  sideCount = 0
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  DO iLocSide = 1,nlocSides
    iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
    ! Big mortar side
    ! 104:bilinear,204:curved
    IF (SideInfo_Shared(SIDE_TYPE,iSide).LE.100) THEN
      sideCount = sideCount + 1
      SideInfo_Shared(SIDE_LOCALID,iSide) = sideCount
    ELSE
      ! Mortar case
      SideInfo_Shared(SIDE_LOCALID,iSide) = -1
    END IF
    ! Check all sides on the small element side to find the small mortar side pointing back
    NbElemID    = SideInfo_Shared(SIDE_NBELEMID,iSide)
    IF(NbElemID.EQ.0) THEN
      SideInfo_Shared(SIDE_NBSIDEID,iSide) = 0
    ELSE IF (NbElemID.LE.-1) THEN
      SideInfo_Shared(SIDE_NBSIDEID,iSide) = -1
    ELSE
      nlocSidesNb = ElemInfo_Shared(ELEM_LASTSIDEIND,NbElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID)
      DO jLocSide = 1,nlocSidesNb
        NbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID) + jLocSide
        IF (ABS(SideInfo_Shared(SIDE_ID,iSide)).EQ.ABS(SideInfo_Shared(SIDE_ID,NbSideID))) THEN
          SideInfo_Shared(SIDE_NBSIDEID,iSide) = NbSideID
          EXIT
        END IF
      END DO
    END IF
  END DO
END DO

#if USE_MPI
! Perform second communication step to distribute updated SIDE_LOCALID
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)*recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
END IF

! Write compute-node local SIDE_NBELEMTYPE
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  SideInfo_Shared(SIDE_NBELEMTYPE,:) = 0
END IF
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
SideInfo_Shared(SIDE_NBELEMTYPE,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp

! final sync of all mesh shared arrays
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeCoords_Shared_Win,MPI_COMM_SHARED)
#endif  /*USE_MPI*/

nUniqueGlobalNodes = MAXVAL(NodeInfo_Shared)
SDEALLOCATE(SideInfo_Shared_tmp)

GETTIME(EndT)
CommMeshReadinWallTime=EndT-StartT
CALL DisplayMessageAndTime(CommMeshReadinWallTime, 'DONE!')

END SUBROUTINE FinishCommunicateMeshReadin


END MODULE MOD_Particle_Mesh_Readin
