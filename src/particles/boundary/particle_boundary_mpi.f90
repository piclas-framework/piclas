#include "boltzplatz.h"

MODULE MOD_part_boundary_mpi
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! PUBLIC
#ifdef MPI
!  INTEGER,ALLOCATABLE           :: NodeIndex(:)
  PUBLIC                        :: Initialize
!  PUBLIC                        :: processParallelBoundary
!===================================================================================================================================

CONTAINS

SUBROUTINE Initialize()
!===================================================================================================================================
! This Routine should be called from DomainUpdate so that for each new
! domain-decomposition the boundaries are updated
! What this Routine does is it collects information from its MPI neighbors about
! reflecting or periodic boundaries in the vicinity of the MPI boundaries. The
! geometrical information about these boundary sides is then saved in lists to
! be used once a particle crosses an MPI boundary near the reflecting or
! periodic boundary.
! The boundary sides are sorted according to the FastInitBackgoundMesh cells so
! that only those boundary sides will have to be checked (whether the particle
! has crossed them) that are associated to the FastInitBackgoundMesh cell
! crossed by the particle.
!===================================================================================================================================
! MODULES
  USE MOD_Globals!,                ONLY : Logging,UNIT_logout
  USE MOD_Mesh_Vars,              ONLY : nNodes, nSides, nBCSides, nInnerSides, nElems
  USE MOD_MPI_Vars,               ONLY : nNbProcs, offsetMPISides_MINE, offsetMPISides_YOUR
  USE MOD_Particle_Vars,          ONLY : GEO
  USE MOD_part_MPI_Vars,          ONLY : MPIGEO, PMPIVAR, PMPIExchange, FIBGMCellPadding, &
                                         SafetyFactor, halo_eps
  USE MOD_CalcTimeStep,           ONLY:CalcTimeStep
  !USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_Equation_Vars,          ONLY : c
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                :: iProc
  INTEGER                :: MPIProcPIC
  INTEGER                :: ALLOCSTAT
  INTEGER                :: MPIGROUPMAP
  INTEGER                :: nConnectedInnerSides(0:PMPIVAR%nProcs-1,1:2)
  INTEGER                :: nConnectedMPISides(0:PMPIVAR%nProcs-1,1:2)
  INTEGER                :: nBoundarySides(0:PMPIVAR%nProcs-1,1:2)
  INTEGER                :: nGhostElements(0:PMPIVAR%nProcs-1,1:2)
  INTEGER,ALLOCATABLE    :: NodeIndex(:)
!===================================================================================================================================

  !--- Each proc receives all MPI sides that lie in FIBGM cells within an 
  !    eps-distance of FIBGM cells containing any of my own MPI sides.
  !--- Check whether old geometry information from other procs exists.
  !    If yes, delete it
  IF (ASSOCIATED(MPIGEO)) THEN 
    IF (ALLOCATED(MPIGEO%ElemSideNodeID)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%ElemSideNodeID,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%NodeCoords)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%NodeCoords,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%ElemToSide)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%ElemToSide,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%SideToElem)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%SideToElem,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%NativeElemID)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%NativeElemID,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%BC)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%BC,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%ElemMPIID)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%ElemMPIID,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%PeriodicElemSide)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%PeriodicElemSide,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF
    IF (ALLOCATED(MPIGEO%ConcaveElemSide)) THEN
      LOGWRITE(*,*)'  - Deleting existing halo Cells...'
      DEALLOCATE(MPIGEO%ConcaveElemSide,STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate halo Cell array',999,999.)
      LOGWRITE(*,*)'    ...Done'
    END IF

  ELSE 
    ALLOCATE(MPIGEO)
  END IF

  ALLOCATE(NodeIndex(1:nNodes))
  ALLOCATE(MPIGEO%haloMPINbSide(1:(nSides-nBCSides-nInnerSides)))
  MPIGEO%haloMPINbSide = 0

  DO iProc=0,PMPIVAR%nProcs-1
    IF (iProc.EQ.PMPIVAR%iProc) CYCLE
    LOGWRITE(*,*)'  - Identify non-immediate MPI-Neighborhood...'
    !--- Search for sides in the neighborhood of non-immediate MPI neighbors
    !--- AS: identifies which of my node have to be sent to iProc w.r.t. to 
    !        eps vicinity region.
    CALL IdentifyNonImmediateMPINeighborhood(iProc,NodeIndex)
    LOGWRITE(*,*)'    ...Done'

    LOGWRITE(*,*)'  - Exchange Geometry of MPI-Neighborhood...'
    CALL ExchangeMPINeighborhoodGeometry(iProc,NodeIndex)
    LOGWRITE(*,*)'    ...Done'
    NodeIndex(:)=0
  END DO

  DEALLOCATE(NodeIndex)

  !--- Setup PMPIVAR%MPINeighbor array, marking each MPI Proc with which
  !    I have to exchange particle data each time step
!  PMPIVAR%MPINeighbor(:)=.FALSE.
!  DO iProc=0,PMPIVAR%nProcs-1
!    IF (MPIGEO%nElems(iProc).GT.0) PMPIVAR%MPINeighbor(iProc)=.TRUE. ! ATTENTION: EXTRA CHECK REQUIRED FOR SHAPEFUNCTION!
!  END DO

  !--- Additional error checking
!  nConnectedInnerSides(0:PMPIVAR%nProcs-1,1:2) = 0
!  nConnectedMPISides(0:PMPIVAR%nProcs-1,1:2)   = 0
!  nBoundarySides(0:PMPIVAR%nProcs-1,1:2)       = 0
!  nGhostElements(0:PMPIVAR%nProcs-1,1:2)       = 0
!  DO iProc=0,PMPIVAR%nProcs-1
!    IF (iProc.EQ.PMPIVAR%iProc) THEN
!      Elem=>MESH%firstElem
!    ELSE IF (ASSOCIATED(FirstElem(iProc)%ep)) THEN
!      Elem=>FirstElem(iProc)%ep
!    ELSE
!      LOGWRITE(*,'(A,I3,A)')'Proc ',iProc,'    :   InnerSides     MPISides BoundarySides'
!      LOGWRITE(*,'(A,3(1X,I12))')'Connected   :', 0,0,0
!      LOGWRITE(*,'(A,3(1X,I12))')'Disconnected:', 0,0,0
!      CYCLE
!    END IF
!    DO WHILE (ASSOCIATED(Elem))
!      nGhostElements(iProc,1) = nGhostElements(iProc,1)+1
!      Side => Elem%firstSide
!      DO WHILE (ASSOCIATED(Side))
!        IF (ASSOCIATED(Side%MPIConnect)) THEN
!          IF (ASSOCIATED(Side%MPIConnect%Connection)) THEN
!            nConnectedMPISides(iProc,1)=nConnectedMPISides(iProc,1)+1
!          ELSE
!            nConnectedMPISides(iProc,2)=nConnectedMPISides(iProc,2)+1
!          END IF
!        ELSE IF (ASSOCIATED(Side%BC)) THEN
!          IF (Side%BC%PartBCType.EQ.PIC%PartBound%MPINeighborhoodBC) THEN
!            nBoundarySides(iProc,2)=nBoundarySides(iProc,2)+1
!          ELSE
!            nBoundarySides(iProc,1)=nBoundarySides(iProc,1)+1
!          END IF
!        ELSE IF (ASSOCIATED(Side%Connection)) THEN
!          nConnectedInnerSides(iProc,1)=nConnectedInnerSides(iProc,1)+1
!        ELSE
!          nConnectedInnerSides(iProc,2)=nConnectedInnerSides(iProc,2)+1
!        END IF
!        Side => Side%nextElemSide
!      END DO
!      Elem => Elem%nextElem
!    END DO
!    LOGWRITE(*,'(A,I3,A)')'Proc ',iProc,'    :   InnerSides     MPISides BoundarySides'
!    LOGWRITE(*,'(A,3(1X,I12))')'Connected   :', &
!                                             nConnectedInnerSides(iProc,1), nConnectedMPISides(iProc,1), nBoundarySides(iProc,1)
!    LOGWRITE(*,'(A,3(1X,I12))')'Disconnected:', &
!                                             nConnectedInnerSides(iProc,2), nConnectedMPISides(iProc,2), nBoundarySides(iProc,2)
!  END DO
!
!  IF(CALC%Logging) FLUSH(UNIT_LogOut)
!
!  nGhostElements(PMPIVAR%iProc,1) = 0
!  DO iProc=0,PMPIVAR%nProcs-1
!    IF (PMPIVAR%iProc.LT.iProc) THEN
!      CALL MPI_SEND(nGhostElements(iProc,1),1,MPI_INTEGER,iProc,1191,PMPIVAR%COMM,IERROR)
!      CALL MPI_RECV(nGhostElements(iProc,2),1,MPI_INTEGER,iProc,1191,PMPIVAR%COMM,MPISTATUS,IERROR)
!    ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
!      CALL MPI_RECV(nGhostElements(iProc,2),1,MPI_INTEGER,iProc,1191,PMPIVAR%COMM,MPISTATUS,IERROR)
!      CALL MPI_SEND(nGhostElements(iProc,1),1,MPI_INTEGER,iProc,1191,PMPIVAR%COMM,IERROR)
!    END IF
!    IF      ((nGhostElements(iProc,1).GT.0).AND.(nGhostElements(iProc,2).EQ.0)) THEN
!      CALL abort(__STAMP__, &
!                 'MPI Neighbor mismatch (1)!',iProc,999.)
!    ELSE IF ((nGhostElements(iProc,2).GT.0).AND.(nGhostElements(iProc,1).EQ.0)) THEN
!      CALL abort(__STAMP__, &
!                 'MPI Neighbor mismatch (2)!',iProc,999.)
!    ELSE IF ((nGhostElements(iProc,1).EQ.0).AND.(PMPIVAR%MPINeighbor(iProc))) THEN
!      CALL abort(__STAMP__, &
!                 'MPI Neighbor mismatch (3)!',iProc,999.)
!    ELSE IF ((nGhostElements(iProc,1).GT.0).AND.(.NOT.PMPIVAR%MPINeighbor(iProc))) THEN
!      CALL abort(__STAMP__, &
!                 'MPI Neighbor mismatch (4)!',iProc,999.)
!    END IF
!  END DO
!
!  DO iProc=0,PMPIVAR%nProcs-1
!    IF (iProc.EQ.PMPIVAR%iProc) THEN
!!      CALL WriteDebugMesh(MESH%firstElem,iProc)
!    ELSE IF (ASSOCIATED(firstElem(iProc)%ep)) THEN
!!      CALL WriteDebugMesh(firstElem(iProc)%ep,iProc)
!    END IF
!  END DO
!  CALL WriteDebugFIBGM(999)

  SWRITE(*,*)'...DONE'
  SWRITE(*,*)'Allocating send and receive requests and message buffers for communicated particle numbers...'
  CALL Finalize()
  !--- allocate message buffers and request buffers
  ALLOCATE(PMPIExchange%send_message(0:PMPIVAR%nProcs-1),&
           PMPIExchange%send_request(0:PMPIVAR%nProcs-1,1:2),&
           PMPIExchange%nbrOfSendParticles(0:PMPIVAR%nProcs-1,1:2),&
           PMPIExchange%recv_request(0:PMPIVAR%nProcs-1,1:2),STAT=ALLOCSTAT)
           !PMPIExchange%recv_message(0:PMPIVAR%nProcs-1),&
           !PMPIExchange%nbrOfRecvParticles(0:PMPIVAR%nProcs-1),STAT=ALLOCSTAT)
  IF (allocStat .NE. 0) THEN
    WRITE(*,*)'ERROR in MPIParallel Initialize: cannot allocate PMPIExchange variables!'; STOP
  END IF
  !PMPIExchange%nbrOfRecvParticles(:)=0
  PMPIExchange%nbrOfSendParticles(:,:)=0
  SWRITE(*,*)'...DONE'

  !--- In case of Shape Function Deposition, messages are sent to my own proc also.
!  IF (PIC%DepositionType.EQ.'shape_function') PMPIVAR%MPINeighbor(PMPIVAR%iProc)=.TRUE.

  RETURN
END SUBROUTINE Initialize

SUBROUTINE ExchangeMPINeighborhoodGeometry(iProc,NodeIndex)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_part_MPI_Vars, ONLY : MPIGEO, PMPIVAR
  USE MOD_Mesh_Vars    , ONLY : nElems, nNodes, nSides, nBCSides, nInnerSides, ElemToSide, &
                                SideToElem, BC
  USE MOD_Particle_Vars, ONLY : GEO, PartBound
  USE MOD_MPI_Vars     , ONLY : offsetElemMPI, nNbProcs, NbProc, &
                                offsetMPISides_MINE, offsetMPISides_YOUR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)             :: iProc       ! MPI proc with which the local proc is to exchange boundary information
  INTEGER, INTENT(INOUT)          :: NodeIndex(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  TYPE tMPISideMessage
    REAL,ALLOCATABLE          :: NodeCoords(:,:)    
    INTEGER,ALLOCATABLE       :: ElemToSide(:,:,:) 
    INTEGER,ALLOCATABLE       :: ElemSideNodeID(:,:,:)
    INTEGER,ALLOCATABLE       :: SideToElem(:,:)
    INTEGER,ALLOCATABLE       :: BC(:,:)
    INTEGER,ALLOCATABLE       :: NativeElemID(:)
    INTEGER,ALLOCATABLE       :: PeriodicElemSide(:,:)
    LOGICAL,ALLOCATABLE       :: ConcaveElemSide(:,:)
    INTEGER                   :: nSides                 ! number of sides to send
    INTEGER                   :: nElems                 ! number of elems to send
    INTEGER                   :: nNodes                 ! number of nodes to send
  END TYPE
  TYPE(tMPISideMessage)       :: SendMsg
  TYPE(tMPISideMessage)       :: RecvMsg
  INTEGER                     :: iElem,iSide,iLocSide,jSide,SideID,iNode,jNode,jProc,iBC
  INTEGER                     :: nNeighborhoodNodes, iIndex, iNbProc, CurrentNbProc
  INTEGER                     :: BCalphaInd
  INTEGER                     :: connectionInd
  INTEGER                     :: ALLOCSTAT
  REAL                        :: t1(1:3),t2(1:3)
  LOGICAL                     :: MPINeighborhood
  INTEGER, ALLOCATABLE        :: TEMPARRAY0(:), TEMPARRAY(:,:), TEMPARRAY2(:,:,:)
  REAL   , ALLOCATABLE        :: TEMPARRAY_R(:,:)
  LOGICAL, ALLOCATABLE        :: TEMPARRAYL(:,:)
  INTEGER, ALLOCATABLE        :: ElemIndex(:)
  INTEGER, ALLOCATABLE        :: SideIndex(:)
  INTEGER                     :: HALOoffsetNode, HALOoffsetSide, HALOoffsetElem
  INTEGER                     :: countDoku
!===================================================================================================================================

  ALLOCATE(ElemIndex(1:nElems), SideIndex(1:nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate Elem- or SideIndex',999,999.)
  SideIndex = 0
  ElemIndex = 0

  !--- First, count marker node indices (nNeighborhoodNodes are within eps distance of at least one MPI-neighbor's node)
  nNeighborhoodNodes=MAXVAL(NodeIndex)
  SendMsg%nNodes=0
  !--- For each MPI neighbor, identify the number of nodes, sides (trias and quads) and elements to be sent
  SendMsg%nElems=0
  SendMsg%nSides=0
  SendMsg%nNodes=nNeighborhoodNodes
  LOGWRITE(*,*)'nNeighborhoodNodes=',nNeighborhoodNodes
  DO iElem=1,nElems
    MPINeighborhood=.FALSE.
    DO iNode=1,8
      IF ( (NodeIndex(GEO%ElemToNodeID(iNode,iElem)).LE.nNeighborhoodNodes) .AND. &
           (NodeIndex(GEO%ElemToNodeID(iNode,iElem)).GT.0)                ) MPINeighborhood=.TRUE.
    END DO
    IF (MPINeighborhood) THEN
      !--- update element-count
      SendMsg%nElems=SendMsg%nElems+1
      ElemIndex(iElem) = SendMsg%nElems
      DO iLocSide = 1,6
        IF (SideIndex(ElemToSide(E2S_SIDE_ID,iLocSide,iElem)).EQ.0) THEN
          !--- update side-count
          SendMsg%nSides=SendMsg%nSides+1
          SideIndex(ElemToSide(E2S_SIDE_ID,iLocSide,iElem)) = SendMsg%nSides
        END IF
      END DO
      !--- name nodes and update node-count
      DO iNode=1,8
        IF (NodeIndex(GEO%ElemToNodeID(iNode,iElem)).EQ.0) THEN
          SendMsg%nNodes=SendMsg%nNodes+1
          NodeIndex(GEO%ElemToNodeID(iNode,iElem)) = SendMsg%nNodes
        END IF
      END DO
    END IF
  END DO
  !WRITE(*,*) "Nodes:", SendMsg%nNodes,"Sides:",SendMsg%nSides,"elems:",SendMsg%nElems,"iProc",PMPIVAR%iProc
  !--- Communicate number of sides (trias,quads), elems (tets,hexas) and nodes to each MPI proc
  
  IF (PMPIVAR%iProc.LT.iProc) THEN
    CALL MPI_SEND(SendMsg%nElems,1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,IERROR)
    CALL MPI_SEND(SendMsg%nSides,1,MPI_INTEGER,iProc,1102,PMPIVAR%COMM,IERROR)
    CALL MPI_SEND(SendMsg%nNodes,1,MPI_INTEGER,iProc,1103,PMPIVAR%COMM,IERROR)
    CALL MPI_RECV(RecvMsg%nElems,1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)
    CALL MPI_RECV(RecvMsg%nSides,1,MPI_INTEGER,iProc,1102,PMPIVAR%COMM,MPISTATUS,IERROR)
    CALL MPI_RECV(RecvMsg%nNodes,1,MPI_INTEGER,iProc,1103,PMPIVAR%COMM,MPISTATUS,IERROR)
  ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
    CALL MPI_RECV(RecvMsg%nElems,1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)
    CALL MPI_RECV(RecvMsg%nSides,1,MPI_INTEGER,iProc,1102,PMPIVAR%COMM,MPISTATUS,IERROR)
    CALL MPI_RECV(RecvMsg%nNodes,1,MPI_INTEGER,iProc,1103,PMPIVAR%COMM,MPISTATUS,IERROR)
    CALL MPI_SEND(SendMsg%nElems,1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,IERROR)
    CALL MPI_SEND(SendMsg%nSides,1,MPI_INTEGER,iProc,1102,PMPIVAR%COMM,IERROR)
    CALL MPI_SEND(SendMsg%nNodes,1,MPI_INTEGER,iProc,1103,PMPIVAR%COMM,IERROR)
  END IF

  !--- allocate send buffers for nodes, sides and elements for each MPI neighbor
  !--- ElemToSide Mapping ------------------------------------------------------!
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%ElemToSide(1:2,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%ElemToSide',SendMsg%nElems,999.)
    SendMsg%ElemToSide(:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ElemToSide(1:2,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems,999.)
    RecvMsg%ElemToSide(:,:,:)=0
  END IF

  !--- NodeCoords --------------------------------------------------------------!
  IF (SendMsg%nNodes.GT.0) THEN
    ALLOCATE(SendMsg%NodeCoords(1:3,1:SendMsg%nNodes),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%NodeCoords',SendMsg%nNodes,999.)
    SendMsg%NodeCoords(:,:)=0.
  END IF
  IF (RecvMsg%nNodes.GT.0) THEN
    ALLOCATE(RecvMsg%NodeCoords(1:3,1:RecvMsg%nNodes),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%NodeCoords',RecvMsg%nNodes,999.)
    RecvMsg%NodeCoords(:,:)=0.
  END IF

  !--- ElemSideNodeID Mapping ------------------------------------------------------!
  IF (SendMsg%nElems.GT.0) THEN       ! ElemSideNodeID(1:iNode,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%ElemSideNodeID(1:4,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save NodeID for each iLocSide  
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%ElemSideNodeID',SendMsg%nElems,999.)
    SendMsg%ElemSideNodeID(:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ElemSideNodeID(1:4,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%ElemSideNodeID',RecvMsg%nElems,999.)
    RecvMsg%ElemSideNodeID(:,:,:)=0
  END IF

  !--- SideToElem Mapping ------------------------------------------------------!
  IF (SendMsg%nSides.GT.0) THEN       ! SideToElem(1:2,1:nSides) 
    ALLOCATE(SendMsg%SideToElem(1:5,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%SideToElem',SendMsg%nSides,999.)
    SendMsg%SideToElem(:,:)=0
  END IF
  IF (RecvMsg%nSides.GT.0) THEN
    ALLOCATE(RecvMsg%SideToElem(1:5,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%SideToElem',RecvMsg%nSides,999.)
    RecvMsg%SideToElem(:,:)=0
  END IF

  !--- BC Mapping ------------------------------------------------------!
  ! BC(1:4,1:nSides) 
  ! 1:BC, 
  ! 2:NBProc, 
  ! 3:1=Mine/2=Yours 
  ! 4:SIDE_ID-MPI_Offset(NBProc)  
  IF (SendMsg%nSides.GT.0) THEN       
    ALLOCATE(SendMsg%BC(1:4,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%BC',SendMsg%nSides,999.)
    SendMsg%BC(:,:)=0
  END IF
  IF (RecvMsg%nSides.GT.0) THEN
    ALLOCATE(RecvMsg%BC(1:4,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%BC',RecvMsg%nSides,999.)
    RecvMsg%BC(:,:)=0
  END IF

  !--- NativeElemID ------------------------------------------------------!
  IF (SendMsg%nElems.GT.0) THEN 
    ALLOCATE(SendMsg%NativeElemID(1:SendMsg%nElems),STAT=ALLOCSTAT)  
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%NativeElemID',SendMsg%nElems,999.)
    SendMsg%NativeElemID(:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%NativeElemID(1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%NativeElemID',RecvMsg%nElems,999.)
    RecvMsg%NativeElemID(:)=0
  END IF

  !--- PeriodicElemSide Mapping ------------------------------------------------------!
  IF (SendMsg%nElems.GT.0) THEN       ! PeriodicElemSide(1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%PeriodicElemSide(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%PeriodicElemSide',SendMsg%nElems,999.)
    SendMsg%PeriodicElemSide(:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%PeriodicElemSide(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%PeriodicElemSide',RecvMsg%nElems,999.)
    RecvMsg%PeriodicElemSide(:,:)=0
  END IF
  !--- ConcaveElemSide Mapping ------------------------------------------------------!
  IF (SendMsg%nElems.GT.0) THEN       ! PeriodicElemSide(1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%ConcaveElemSide(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate SendMsg%PeriodicElemSide',SendMsg%nElems,999.)
    SendMsg%ConcaveElemSide(:,:)=.FALSE.
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ConcaveElemSide(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate RecvMsg%PeriodicElemSide',RecvMsg%nElems,999.)
    RecvMsg%ConcaveElemSide(:,:)=.FALSE.
  END IF

  !--- fill send buffers with node, side and element data (including connectivity!)
  !--- ElemtoSide --------------------------------------------------------------!
  DO iElem = 1,nElems
    IF (ElemIndex(iElem).NE.0) THEN
      DO iLocSide = 1,6
        SendMsg%ElemToSide(1,iLocSide,ElemIndex(iElem)) = &
                 SideIndex(ElemToSide(E2S_SIDE_ID,iLocSide,iElem))
        SendMsg%ElemToSide(2,iLocSide,ElemIndex(iElem)) = &
                ElemToSide(2,iLocSide,iElem)
      END DO
    END IF
  END DO
  !--- NodeCoords --------------------------------------------------------------!
  DO iNode = 1,nNodes
    IF (NodeIndex(iNode).NE.0) THEN
      SendMsg%NodeCoords(1:3,NodeIndex(iNode))=GEO%NodeCoords(1:3,iNode)
    END IF
  END DO
  !--- ElemSideNodeID Mapping ------------------------------------------------------!
  DO iElem = 1,nElems
    IF (ElemIndex(iElem).NE.0) THEN
      DO iLocSide = 1,6
        DO iNode = 1,4
          SendMsg%ElemSideNodeID(iNode,iLocSide,ElemIndex(iElem)) = &
             NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,iElem))
        END DO
      END DO
    END IF
  END DO 
  !--- SideToElem Mapping ------------------------------------------------------!
  DO iSide = 1,nSides
    IF (SideIndex(iSide).NE.0) THEN
      DO iIndex = 1,2   ! S2E_ELEM_ID, S2E_NB_ELEM_ID
        IF (SideToElem(iIndex,iSide).GT.0) THEN
           SendMsg%SideToElem(iIndex,SideIndex(iSide)) = &
                   ElemIndex(SideToElem(iIndex,iSide))     ! gets 0 if pointing to a non-Halo-cell.
                                                           ! if 0 and no BC is associated this 
                                                           ! is a Halo-bound-cell
        END IF
      END DO ! S2E_LOC_SIDE_ID, S2E_NB_LOC_SIDE_ID, S2E_FLIP
      SendMsg%SideToElem(3:5,SideIndex(iSide)) = &
          SideToElem(3:5,iSide)
    END IF
  END DO
  !--- BC Mapping ------------------------------------------------------!
  DO iSide = 1,nBCSides  ! no need to go through all side since BC(1:nBCSides)
    IF (SideIndex(iSide).NE.0) THEN
      SendMsg%BC(1,SideIndex(iSide)) = BC(iSide)
    END IF
  END DO
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  DO iSide = nBCSides+nInnerSides+1,nSides ! only go through MPI-Sides 
    IF(SideIndex(iSide).NE.0) THEN
      DO iNbProc =1,nNbProcs  ! ??? iNbProc=1 start at nBCSides+nInnerSides ???
        ! first check "mine"-sides:
        SendMsg%BC(1,SideIndex(iSide)) = -1 
        IF ((iSide.GE.offsetMPISides_MINE(iNbProc-1)+1).AND. &
            (iSide.LE.offsetMPISides_MINE(iNbProc))) THEN
            SendMsg%BC(2,SideIndex(iSide)) = NbProc(iNbProc)  ! global NbProcID
            SendMsg%BC(3,SideIndex(iSide)) = 1                ! 1=Mine/2=Yours  
            SendMsg%BC(4,SideIndex(iSide)) = iSide-offsetMPISides_MINE(iNbProc-1) ! MPITag 
        ! then check "your"-sides:
        ELSE IF ((iSide.GE.offsetMPISides_YOUR(iNbProc-1)+1).AND. &
            (iSide.LE.offsetMPISides_YOUR(iNbProc))) THEN
            SendMsg%BC(2,SideIndex(iSide)) = NbProc(iNbProc)  ! global NbProcID
            SendMsg%BC(3,SideIndex(iSide)) = 2                ! 1=Mine/2=Yours  
            SendMsg%BC(4,SideIndex(iSide)) = iSide-offsetMPISides_YOUR(iNbProc-1) ! MPITag 
        END IF
      END DO
    END IF
  END DO
  ! proceed MPI sides
  DO iSide = nBCSides+1,nBCSides+nInnerSides
    IF(SideIndex(iSide).NE.0) THEN
      IF ((SendMsg%SideToElem(1,SideIndex(iSide)).EQ.0).OR. &
          (SendMsg%SideToElem(2,SideIndex(iSide)).EQ.0)) THEN
        SendMsg%BC(1,SideIndex(iSide)) = 424242        ! EPS-Region-Bound in innersides
      END IF
    END IF
  END DO
  !--- NativeElemID ------------------------------------------------------!
  DO iElem = 1,nElems
    IF (ElemIndex(iElem).NE.0) THEN
      SendMsg%NativeElemID(ElemIndex(iElem)) = iElem
    END IF
  END DO 

 !--- PeriodicElemSide ------------------------------------------------------!
  IF(GEO%nPeriodicVectors.GT.0)THEN
    DO iElem = 1,nElems
      IF (ElemIndex(iElem).NE.0) THEN
        DO iLocSide = 1,6
          SendMsg%PeriodicElemSide(iLocSide,ElemIndex(iElem)) = &
               GEO%PeriodicElemSide(iLocSide,iElem)
        END DO
      END IF
    END DO
  END IF
 !--- ConcaveElemSide ------------------------------------------------------!
 DO iElem = 1,nElems
   IF (ElemIndex(iElem).NE.0) THEN
     DO iLocSide = 1,6
       SendMsg%ConcaveElemSide(iLocSide,ElemIndex(iElem)) = &
            GEO%ConcaveElemSide(iLocSide,iElem)
     END DO
   END IF
 END DO


  IF (PMPIVAR%iProc.LT.iProc) THEN
    ! Send:
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nNodes.GT.0) CALL MPI_SEND(SendMsg%NodeCoords,SendMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1105,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemSideNodeID,SendMsg%nElems*4*6,MPI_INTEGER   ,iProc,1106,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1107,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides*4,MPI_INTEGER                 ,iProc,1108,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1109,PMPIVAR%COMM,IERROR)
    IF(GEO%nPeriodicVectors.GT.0)THEN
      IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%PeriodicElemSide,SendMsg%nElems*6,MPI_INTEGER ,iProc,1110,PMPIVAR%COMM,IERROR)
    END IF
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ConcaveElemSide,SendMsg%nElems*6,MPI_LOGICAL ,iProc,1111,PMPIVAR%COMM,IERROR)
    ! Receive:
    IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nNodes.GT.0) &
      CALL MPI_RECV(RecvMsg%NodeCoords,RecvMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1105,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemSideNodeID,RecvMsg%nElems*4*6,MPI_INTEGER   ,iProc,1106,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1107,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides*4,MPI_INTEGER                 ,iProc,1108,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1109,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF(GEO%nPeriodicVectors.GT.0)THEN
      IF (RecvMsg%nElems.GT.0) &
           CALL MPI_RECV(RecvMsg%PeriodicElemSide,RecvMsg%nElems*6,MPI_INTEGER,iProc,1110,PMPIVAR%COMM,MPISTATUS,IERROR)
    END IF
    IF (RecvMsg%nElems.GT.0) &
         CALL MPI_RECV(RecvMsg%ConcaveElemSide,RecvMsg%nElems*6,MPI_LOGICAL,iProc,1111,PMPIVAR%COMM,MPISTATUS,IERROR)
  ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
    ! Receive:
    IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nNodes.GT.0) &
      CALL MPI_RECV(RecvMsg%NodeCoords,RecvMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1105,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemSideNodeID,RecvMsg%nElems*4*6,MPI_INTEGER   ,iProc,1106,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1107,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides*4,MPI_INTEGER                 ,iProc,1108,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1109,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF(GEO%nPeriodicVectors.GT.0)THEN
      IF (RecvMsg%nElems.GT.0) &
           CALL MPI_RECV(RecvMsg%PeriodicElemSide,RecvMsg%nElems*6,MPI_INTEGER,iProc,1110,PMPIVAR%COMM,MPISTATUS,IERROR)
    END IF
    IF (RecvMsg%nElems.GT.0) &
         CALL MPI_RECV(RecvMsg%ConcaveElemSide,RecvMsg%nElems*6,MPI_LOGICAL,iProc,1111,PMPIVAR%COMM,MPISTATUS,IERROR)
    ! Send:
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nNodes.GT.0) CALL MPI_SEND(SendMsg%NodeCoords,SendMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1105,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemSideNodeID,SendMsg%nElems*4*6,MPI_INTEGER   ,iProc,1106,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1107,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides*4,MPI_INTEGER                 ,iProc,1108,PMPIVAR%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1109,PMPIVAR%COMM,IERROR)
    IF(GEO%nPeriodicVectors.GT.0)THEN
      IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%PeriodicElemSide,SendMsg%nElems*6,MPI_INTEGER ,iProc,1110,PMPIVAR%COMM,IERROR)
    END IF
    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ConcaveElemSide,SendMsg%nElems*6,MPI_LOGICAL ,iProc,1111,PMPIVAR%COMM,IERROR)
  END IF

!  IF (SendMsg%nElems.GT.0) LOGWRITE(*,*)'SendElems=',SendMsg%Elems
!  IF (RecvMsg%nElems.GT.0) LOGWRITE(*,*)'RecvElems=',RecvMsg%Elems
!  IF (SendMsg%nSides.GT.0) WRITE(UNIT_LogOut,*)'SendSides(',iProc,')=',&
!     SendMsg%Sides(1:SendMsg%nSides,6:7)
!  IF (RecvMsg%nSides.GT.0) WRITE(UNIT_LogOut,*)'RecvSides(',iProc,')=',&
!     RecvMsg%Sides(1:RecvMsg%nSides,6:7)
!  IF (SendMsg%nNodes.GT.0) LOGWRITE(*,*)'SendNodes=',SendMsg%Nodes
!  IF (RecvMsg%nNodes.GT.0) LOGWRITE(*,*)'RecvNodes=',RecvMsg%Nodes

  IF ((RecvMsg%nElems.EQ.0) .AND. ((RecvMsg%nSides.GT.0).OR.(RecvMsg%nNodes.GT.0))) THEN
      ERRWRITE(*,*)'ERROR: nElems=0 when nSides=',RecvMsg%nSides,' and nNodes=',RecvMsg%nNodes,'!'
      CALL abort(__STAMP__,'nElems=0 while nNodes=',RecvMsg%nNodes,999.)
  END IF
  IF (.NOT.RecvMsg%nNodes.EQ.0) THEN
  IF (ALLOCATED(MPIGEO%ElemSideNodeID)) THEN   !--- Resize Arrays
    HALOoffsetNode=SIZE(MPIGEO%NodeCoords,2)
    HALOoffsetSide=SIZE(MPIGEO%SideToElem,2)
    HALOoffsetElem=SIZE(MPIGEO%ElemSideNodeID,3)
    ALLOCATE(TEMPARRAY2(1:4,1:6,1:HALOoffsetElem),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (1)!',iProc,999.)
    TEMPARRAY2(1:4,1:6,1:HALOoffsetElem)=MPIGEO%ElemSideNodeID(1:4,1:6,1:HALOoffsetElem)
    DEALLOCATE(MPIGEO%ElemSideNodeID)
    ALLOCATE(MPIGEO%ElemSideNodeID(1:4,1:6,1:HALOoffsetElem+RecvMsg%nElems), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (2)!',iProc,999.)
    MPIGEO%ElemSideNodeID(:,:,:)=0
    MPIGEO%ElemSideNodeID(1:4,1:6,1:HALOoffsetElem)=TEMPARRAY2(1:4,1:6,1:HALOoffsetElem)
    DEALLOCATE(TEMPARRAY2)

    ALLOCATE(TEMPARRAY_R(1:3,1:HALOoffsetNode),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (3)!',iProc,999.)
    TEMPARRAY_R(1:3,1:HALOoffsetNode)=MPIGEO%NodeCoords(1:3,1:HALOoffsetNode)
    DEALLOCATE(MPIGEO%NodeCoords)
    ALLOCATE(MPIGEO%NodeCoords(1:3,1:HALOoffsetNode+RecvMsg%nNodes), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (4)!',iProc,999.)
    MPIGEO%NodeCoords(:,:)=0.
    MPIGEO%NodeCoords(1:3,1:HALOoffsetNode)=TEMPARRAY_R(1:3,1:HALOoffsetNode)
    DEALLOCATE(TEMPARRAY_R)

    ALLOCATE(TEMPARRAY2(1:2,1:6,1:HALOoffsetElem),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (5)!',iProc,999.)
    TEMPARRAY2(1:2,1:6,1:HALOoffsetElem)=MPIGEO%ElemToSide(1:2,1:6,1:HALOoffsetElem)
    DEALLOCATE(MPIGEO%ElemToSide)
    ALLOCATE(MPIGEO%ElemToSide(1:2,1:6,1:HALOoffsetElem+RecvMsg%nElems), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (6)!',iProc,999.)
    MPIGEO%ElemToSide(:,:,:)=0
    MPIGEO%ElemToSide(1:2,1:6,1:HALOoffsetElem)=TEMPARRAY2(1:2,1:6,1:HALOoffsetElem)
    DEALLOCATE(TEMPARRAY2)

    ALLOCATE(TEMPARRAY(1:5,1:HALOoffsetSide),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (7)!',iProc,999.)
    TEMPARRAY(1:5,1:HALOoffsetSide)=MPIGEO%SideToElem(1:5,1:HALOoffsetSide)
    DEALLOCATE(MPIGEO%SideToElem)
    ALLOCATE(MPIGEO%SideToElem(1:5,1:HALOoffsetSide+RecvMsg%nSides), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (8)!',iProc,999.)
    MPIGEO%SideToElem(:,:)=-1  ! consistent with SideToElem
    MPIGEO%SideToElem(1:5,1:HALOoffsetSide)=TEMPARRAY(1:5,1:HALOoffsetSide)
    DEALLOCATE(TEMPARRAY)

    ALLOCATE(TEMPARRAY0(1:HALOoffsetElem),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (9)!',iProc,999.)
    TEMPARRAY0(1:HALOoffsetElem)=MPIGEO%ElemMPIID(1:HALOoffsetElem)
    DEALLOCATE(MPIGEO%ElemMPIID)
    ALLOCATE(MPIGEO%ElemMPIID(1:HALOoffsetElem+RecvMsg%nElems), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (10)!',iProc,999.)
    MPIGEO%ElemMPIID(:)=0
    MPIGEO%ElemMPIID(1:HALOoffsetElem)=TEMPARRAY0(1:HALOoffsetElem)
    DEALLOCATE(TEMPARRAY0)

    ALLOCATE(TEMPARRAY(1:5,1:HALOoffsetSide),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (11)!',iProc,999.)
    TEMPARRAY(1:5,1:HALOoffsetSide)=MPIGEO%BC(1:5,1:HALOoffsetSide)
    DEALLOCATE(MPIGEO%BC)
    ALLOCATE(MPIGEO%BC(1:5,1:HALOoffsetSide+RecvMsg%nSides), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (12)!',iProc,999.)
    MPIGEO%BC(1:5,:)=0
    MPIGEO%BC(1:5,1:HALOoffsetSide)=TEMPARRAY(1:5,1:HALOoffsetSide)
    DEALLOCATE(TEMPARRAY)

    ALLOCATE(TEMPARRAY0(1:HALOoffsetElem),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (9)!',iProc,999.)
    TEMPARRAY0(1:HALOoffsetElem)=MPIGEO%NativeElemID(1:HALOoffsetElem)
    DEALLOCATE(MPIGEO%NativeElemID)
    ALLOCATE(MPIGEO%NativeElemID(1:HALOoffsetElem+RecvMsg%nElems), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (10)!',iProc,999.)
    MPIGEO%NativeElemID(:)=0
    MPIGEO%NativeElemID(1:HALOoffsetElem)=TEMPARRAY0(1:HALOoffsetElem)
    DEALLOCATE(TEMPARRAY0)

    IF(GEO%nPeriodicVectors.GT.0)THEN
      ALLOCATE(TEMPARRAY(1:6,1:HALOoffsetElem),STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (9)!',iProc,999.)
      TEMPARRAY(1:6,1:HALOoffsetElem)=MPIGEO%PeriodicElemSide(1:6,1:HALOoffsetElem)
      DEALLOCATE(MPIGEO%PeriodicElemSide)
      ALLOCATE(MPIGEO%PeriodicElemSide(1:6,1:HALOoffsetElem+RecvMsg%nElems), STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (10)!',iProc,999.)
      MPIGEO%PeriodicElemSide(:,:)=0
      MPIGEO%PeriodicElemSide(1:6,1:HALOoffsetElem)=TEMPARRAY(1:6,1:HALOoffsetElem)
      DEALLOCATE(TEMPARRAY)
    END IF

    ALLOCATE(TEMPARRAYL(1:6,1:HALOoffsetElem),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (9)!',iProc,999.)
    TEMPARRAYL(1:6,1:HALOoffsetElem)=MPIGEO%ConcaveElemSide(1:6,1:HALOoffsetElem)
    DEALLOCATE(MPIGEO%ConcaveElemSide)
    ALLOCATE(MPIGEO%ConcaveElemSide(1:6,1:HALOoffsetElem+RecvMsg%nElems), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not allocate (10)!',iProc,999.)
    MPIGEO%ConcaveElemSide(:,:)=.FALSE.
    MPIGEO%ConcaveElemSide(1:6,1:HALOoffsetElem)=TEMPARRAYL(1:6,1:HALOoffsetElem)
    DEALLOCATE(TEMPARRAYL)


  ELSE
    IF(GEO%nPeriodicVectors.GT.0) ALLOCATE(MPIGEO%PeriodicElemSide(1:6,1:RecvMsg%nElems))
    ALLOCATE(MPIGEO%ElemSideNodeID(1:4,1:6,1:RecvMsg%nElems), &
             MPIGEO%NodeCoords(1:3,1:RecvMsg%nNodes),         &
             MPIGEO%ElemToSide(1:2,1:6,1:RecvMsg%nElems),     &
             MPIGEO%SideToElem(1:5,1:RecvMsg%nSides),         &
             MPIGEO%BC(1:5,1:RecvMsg%nSides),                 &
             MPIGEO%ElemMPIID(1:RecvMsg%nElems),              &
             MPIGEO%NativeElemID(1:RecvMsg%nElems),           &
             MPIGEO%ConcaveElemSide(1:6,1:RecvMsg%nElems),    &
             STAT=ALLOCSTAT                                        )
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__,'Could not allocate receive messages for HALO-cells!',iProc,999.)
    END IF
    MPIGEO%ElemSideNodeID(:,:,:) = 0
    MPIGEO%NodeCoords(:,:)       = 0.
    MPIGEO%ElemToSide(:,:,:)     = 0 
    MPIGEO%SideToElem(:,:)       = -1
    MPIGEO%BC(:,:)               = 0
    MPIGEO%ElemMPIID(:)          = 0
    MPIGEO%NativeElemID(:)       = 0
    MPIGEO%ConcaveElemSide(:,:)  = .FALSE.
    IF (GEO%nPeriodicVectors.GT.0) MPIGEO%PeriodicElemSide(:,:) = 0
    HALOoffsetNode               = 0
    HALOoffsetElem               = 0
    HALOoffsetSide               = 0
  END IF !--- end of resize
 
  ! Fill new grid

  !--- PeriodicElemSide  --------------------------------------------------------------!
  IF(GEO%nPeriodicVectors.GT.0)THEN
    MPIGEO%PeriodicElemSide(1:6,HALOoffsetElem+1:HALOoffsetElem+RecvMsg%nElems) = &
         RecvMsg%PeriodicElemSide(1:6,1:RecvMsg%nElems)
  END IF
  !--- ConcaveElemSide  --------------------------------------------------------------!
  MPIGEO%ConcaveElemSide(1:6,HALOoffsetElem+1:HALOoffsetElem+RecvMsg%nElems) = &
       RecvMsg%ConcaveElemSide(1:6,1:RecvMsg%nElems)
  !--- ElemMPIID  --------------------------------------------------------------!
  MPIGEO%ElemMPIID(HALOoffsetElem+1:HALOoffsetElem+RecvMsg%nElems) = iProc
  !--- NativeElemID  --------------------------------------------------------------!
  MPIGEO%NativeElemID(HALOoffsetElem+1:HALOoffsetElem+RecvMsg%nElems) = &
                      RecvMsg%NativeElemID(1:RecvMsg%nElems)
  !--- NodeCoords --------------------------------------------------------------!
  DO iNode=1,RecvMsg%nNodes
    MPIGEO%NodeCoords(1:3,HALOoffsetNode+iNode) = RecvMsg%NodeCoords(1:3,iNode)
  END DO
  !--- ElemtoSide ------------------------------------------------------!
  DO iElem = 1,RecvMsg%nElems
    DO iLocSide = 1,6
        MPIGEO%ElemToSide(1,iLocSide,HALOoffsetElem+iElem) = &   ! SIDE_ID
           HALOoffsetSide + RecvMsg%ElemToSide(1,iLocSide,iElem)
        MPIGEO%ElemToSide(2,iLocSide,HALOoffsetElem+iElem) = &   ! FLIP
                            RecvMsg%ElemToSide(1,iLocSide,iElem)
      END DO
  END DO
  !--- ElemSideNodeID Mapping ------------------------------------------------------!
  DO iElem = 1,RecvMsg%nElems
    DO iLocSide = 1,6
      DO iNode = 1,4
        MPIGEO%ElemSideNodeID(iNode,iLocSide,HALOoffsetElem+iElem) = &
           HALOoffsetNode + RecvMsg%ElemSideNodeID(iNode,iLocSide,iElem) 
      END DO
    END DO
  END DO 
  !--- SideToElem Mapping ------------------------------------------------------!
  DO iSide = 1,RecvMsg%nSides
    DO iIndex = 1,2   ! S2E_ELEM_ID, S2E_NB_ELEM_ID
      IF(RecvMsg%SideToElem(iIndex,iSide).NE.0) THEN
        MPIGEO%SideToElem(iIndex,HALOoffsetSide+iSide) = &
        HALOoffsetElem + RecvMsg%SideToElem(iIndex,iSide) 
      ELSE 
        ! gets 0 if pointing to a non-Halo-cell.
        ! if 0 and no BC is associated this 
        ! is a Halo-bound-cell
        ! either way it stays at the initialized value of -1
      END IF
    END DO ! S2E_LOC_SIDE_ID, S2E_NB_LOC_SIDE_ID, S2E_FLIP
    MPIGEO%SideToElem(3:5,HALOoffsetSide+iSide) = &
        RecvMsg%SideToElem(3:5,iSide)
  END DO
  !--- BC Mapping ------------------------------------------------------!
  ! first identify the local NbProc of iProc
  DO iNbProc=1,nNbProcs
    IF(iProc.EQ.nbProc(iNbProc)) CurrentNbProc = iNbProc
  END DO

  DO iSide = 1,RecvMsg%nSides
    IF (RecvMsg%BC(1,iSide).GT.0) THEN
      MPIGEO%BC(1,HALOoffsetSide + iSide) = RecvMsg%BC(1,iSide)
    ELSE IF (RecvMsg%BC(1,iSide).EQ.-1) THEN    ! MPI-Side
      MPIGEO%BC(1,HALOoffsetSide+iSide) = -1 
      MPIGEO%BC(2:4,HALOoffsetSide+iSide) = RecvMsg%BC(2:4,iSide)
      MPIGEO%BC(5,HALOoffsetSide+iSide) = iProc  ! Where am I from? 
      IF (MPIGEO%BC(2,HALOoffsetSide+iSide).EQ.PMPIVAR%iProc) THEN ! case 1: my own mpi-side
        ! build mapping from mySides to MPIGEO-Sides
        IF (MPIGEO%BC(3,HALOoffsetSide+iSide).EQ.1) THEN ! MINE-Side for other proc > YOUR-Side for myproc
           MPIGEO%haloMPINbSide(offsetMPISides_YOUR(CurrentNbProc-1)&
                               +MPIGEO%BC(4,HALOoffsetSide+iSide)   &
                               -(nBCSides+nInnerSides)) =           & 
                               HALOoffsetSide+iSide
           ! do inverse mapping from halocells to myproc cells
           MPIGEO%BC(4,HALOoffsetSide+iSide) =-(offsetMPISides_YOUR(CurrentNbProc-1)& 
                                               +MPIGEO%BC(4,HALOoffsetSide+iSide))
        ELSE
           MPIGEO%haloMPINbSide(offsetMPISides_MINE(CurrentNbProc-1)&
                               +MPIGEO%BC(4,HALOoffsetSide+iSide)   &
                               -(nBCSides+nInnerSides)) =           & 
                               HALOoffsetSide+iSide
           ! do inverse mapping from halocells to myproc cells
           MPIGEO%BC(4,HALOoffsetSide+iSide) =-(offsetMPISides_MINE(CurrentNbProc-1)& 
                                               +MPIGEO%BC(4,HALOoffsetSide+iSide))
        END IF
      ELSE ! case 3: mpi-side in halo-cells
        ! make connection between halo-mpi-domains:
        ! only proceed if NbProc is lower than iProc because data for
        ! greater ProcID's are not available yet.
        IF (MPIGEO%BC(2,HALOoffsetSide+iSide).LT.iProc) THEN
          DO jSide = 1,HALOoffsetSide
            IF ((MPIGEO%BC(2,jSide).EQ.iProc) )THEN !.AND.&                             ! check proc
              IF  (MPIGEO%BC(3,jSide).NE.MPIGEO%BC(3,HALOoffsetSide+iSide)) THEN !.AND.& ! needs to be different (mine vs. your)
              IF  (MPIGEO%BC(4,jSide).EQ.MPIGEO%BC(4,HALOoffsetSide+iSide)) THEN ! MPI-Tag
              IF  (MPIGEO%BC(5,jSide).EQ.MPIGEO%BC(2,HALOoffsetSide+iSide)) THEN ! Check for Neighbor-Proc
                ! Link: delete MPI-Tag since not required anymore
                MPIGEO%BC(4,jSide) = -(HALOoffsetSide+iSide)
                MPIGEO%BC(4,HALOoffsetSide+iSide) = -jSide
            END IF
            END IF
            END IF
            END IF
          END DO
        END IF
      END IF
    ELSE
      ! innerSide
    END IF
  END DO

! TBD BCalphaInd  = recvMsg%Sides(iSide,4+4)
!    IF ((iBC.GT.0).AND.(jProc.EQ.-1)) THEN                                     ! Connected side - must find connection
!      IF (MPIGEO%SideToElem(S2E_NB_ELEM_ID,HALOoffsetSide+iSide).NE.0) CYCLE      ! with side%ind=iBC
!      connectionInd = recvMsg%Sides(iSide,4+2)                                 ! neighbor GlobalElemID
!      DO iElem=HALOoffsetElem+1,HALOoffsetElem+RecvMsg%nElems
!        IF (MPIGEO%ElemID(HALOoffsetElem+iElem).EQ.connectionInd) THEN
!          MPIGEO%SideToElem(S2E_NB_ELEM_ID,HALOoffsetSide+iSide)=iElem
!        END IF
!      END DO
!      IF (MPIGEO%SideToElem(S2E_NB_ELEM_ID,HALOoffsetSide+iSide).EQ.0) THEN
!        CALL abort(__STAMP__,&
!                             'Could not match connection of MPI Neighborhood side (inner side)',iSide,999.)
!      END IF
!    ELSE IF (iBC.LT.0) THEN                                                    ! Boundary Side
!      MPIGEO%BC(HALOoffsetSide+iSide) = -iBC
!      IF (Side%BC%PartBCtype.EQ.PIC%PartBound%PeriodicBC) THEN                 ! Periodic Boundary
!        Side%BC%BCalphaInd = BCalphaInd
!        IF ((ABS(BCalphaInd).GT.MESH%nVV).OR.(BCalphaind.EQ.0)) THEN
!          CALL abort(__STAMP__,&
!                               'Periodic Vector invalid!',BCalphaInd,999.)
!        END IF
!        IF (ASSOCIATED(Side%Connection)) CYCLE
!        connectionInd = recvMsg%Sides(iSide,4+3)
!        DO jSide=1,RecvMsg%nSides
!          IF (NeighborhoodData%Sides(jSide)%sp%ind.EQ.connectionInd) THEN
!            Side%Connection=>NeighborhoodData%Sides(jSide)%sp
!            NeighborhoodData%Sides(jSide)%sp%Connection=>Side
!            EXIT
!          END IF
!        END DO
!        IF (.NOT.ASSOCIATED(Side%Connection)) THEN
!          CALL abort(__STAMP__,&
!                               'Could not match connection of MPI Neighborhood side (periodic side)',iSide,999.)
!        END IF
!      END IF
!    ELSE IF ((iBC.GT.0).AND.(jProc.GE.0).AND.(jProc.LT.PMPIVAR%nProcs)) THEN ! MPI Boundary - must find connection
                                                                               ! from iProc=jProc with side%tag=iBC
!      IF (BCalphaInd.NE.0) THEN                                                ! Periodic MPI Boundary
!        IF (.NOT.ASSOCIATED(Side%BC)) THEN
!          CALL getNewBC(Side%BC)
!          Side%BC%PartBCtype=PIC%PartBound%PeriodicBC
!          Side%BC%BCalphaInd = BCalphaInd
!          IF ((ABS(BCalphaInd).GT.MESH%nVV).OR.(BCalphaind.EQ.0)) THEN
!            CALL abort(__STAMP__,&
!                                 'Periodic Vector invalid!',BCalphaInd,999.)
!          END IF
!        ELSE
!          CALL abort(__STAMP__,&
!                               'BC already associated, other than expected!',BCalphaInd,999.)
!        END IF
!      END IF
!      IF (ASSOCIATED(Side%Connection)) THEN
!        LOGWRITE(*,*)'MPI-connection already established with side ind ',Side%Connection%ind
!        CYCLE
!      END IF
!      !--- For local MPI sides, we cannot use side%connection directly to connect the sides.
!      !    Otherwise the DG code would be confused. So we use Side%MPIConnect%Connection
!      !    which exists exclusively for Particle applications
!      IF (jProc.EQ.PMPIVAR%iProc) THEN
!        Side%Connection=>PMPIVAR%MPIConnect(iProc)%tagToSide(iBC)%sp
!        PMPIVAR%MPIConnect(iProc)%tagToSide(iBC)%sp%MPIConnect%Connection=>Side
!        LOGWRITE(*,*)'connected side from proc ',iProc,' to my MPI-side with tag ',&
!                                                                 iBC,' (MPI-connection)!'
!        CYCLE
!      ELSE IF (.NOT.ASSOCIATED(FirstElem(jProc)%ep)) THEN
!!        IPWRITE(*,*)'WARNING: Could not make MPI connection from iproc ',iProc,&
!!                    ' to jproc ',jProc,' because no sides from jproc were communicated to me!'
!!        STOP
!      END IF
!      Elem=>FirstElem(jProc)%ep
!      DO WHILE (ASSOCIATED(Elem).AND.(.NOT.ASSOCIATED(Side%Connection)))
!        connectSide=>Elem%firstSide
!        DO WHILE (ASSOCIATED(connectSide))
!          IF ((connectSide%iLocElemSide.EQ.iBC).AND.(connectSide%nGP.EQ.iProc)) THEN  ! iLocElemSide and nGP are used to save
!                                                                                      ! MPI tag and neighbor proc
!            Side%Connection=>connectSide
!            connectSide%Connection=>Side
!            LOGWRITE(*,*)'connected side from proc ',iProc,' to MPI-side from proc ',&
!                                                     jProc,' tag ',iBC,' (MPI-connection)!'
!            EXIT
!          END IF
!          connectSide=>connectSide%nextElemSide
!        END DO
!        Elem=>Elem%nextElem
!      END DO
!      IF (.NOT.ASSOCIATED(Side%Connection)) THEN
!          CALL abort(__STAMP__,&
!                               'ERROR: Could not make connection of MPI Neighborhood side ',iSide,999.)
!        IPWRITE(*,*)'WARNING: Could not match connection of MPI Neighborhood side ',&
!          iSide,' from proc ',iProc,' to jProc ',jProc,' with tag ',iBC,' (MPI-connection)!'
!        STOP
!      END IF
!    ELSE                                                                       ! Error
!      IPWRITE(*,*)'ERROR: iProc=',iProc,' iSide=',iSide,' iBC=',iBC,', jProc=',jProc,'!'
!      CALL abort(__STAMP__,&
!                           'Could not connect MPI Neighborhood side',iSide,999.)
!    END IF
    
!    !--- assign nodes to side
!    DO iNode=1,Side%nNodes
!      jNode = RecvMsg%Sides(iSide,iNode)
!      IF ((jNode.LE.0).OR.(jNode.GT.RecvMsg%nNodes)) THEN
!        ERRWRITE(*,*)'======================================'
!        ERRWRITE(*,*)'iProc             =',iProc
!        ERRWRITE(*,*)'nNeighborhoodNodes=',nNeighborhoodNodes
!        ERRWRITE(*,*)'SendMsg%nElems    =',SendMsg%nElems
!        ERRWRITE(*,*)'RecvMsg%nElems    =',RecvMsg%nElems
!        ERRWRITE(*,*)'SendMsg%nSides    =',SendMsg%nSides
!        ERRWRITE(*,*)'RecvMsg%nSides    =',RecvMsg%nSides
!        ERRWRITE(*,*)'SendMsg%nNodes    =',SendMsg%nNodes
!        ERRWRITE(*,*)'RecvMsg%nNodes    =',RecvMsg%nNodes
!        ERRWRITE(*,'(A,I0,A,7(1X,I0))')'RecvMsg%Sides(',iSide,',1:7)=',RecvMsg%Sides(iSide,1:7)
!        ERRWRITE(*,*)'======================================'
!        CALL abort(__STAMP__,'jNode=',jNode,999.)
!      END IF
!      Side%node(iNode)%np => NeighborhoodData%Nodes(jNode)%np
!      Side%ind = recvMsg%Sides(iSide,4+1)
!      Side%iLocElemSide = recvMsg%Sides(iSide,4+2)   ! iLocElemSide is abused here to save the MPI tag (if present)
!      Side%nGP          = recvMsg%Sides(iSide,4+3)-1 ! nGP is abused here to save the neighbor proc (if present)
!      IF (Side%ind.LE.0) THEN
!        WRITE(*,*)'ERROR: NeighborhoodData(',iProc,')%Sides(',iSide,')%sp%ind=',Side%ind
!      END IF
!    END DO
    !--- compute side normal
!    t1=Side%Node(2)%np%x-Side%Node(1)%np%x
!    t2=Side%Node(3)%np%x-Side%Node(1)%np%x
!    ALLOCATE(side%normal(1,3,1))
!    Side%Normal(1,:,1)=cross(t1,t2)
!    Side%Normal(1,:,1)=Side%Normal(1,:,1)/sqrt(SUM(Side%Normal(1,:,1)*Side%Normal(1,:,1)))
!    LOGWRITE(*,*)'Creating NeighborhoodData(',iProc,')%Sides(',iSide,')%sp'
!  END DO

  END IF !(.NOT.RecvMsg%nNodes.EQ.0) 
  !--- Deallocate message buffers 
  IF (RecvMsg%nNodes.GT.0) THEN
    DEALLOCATE(RecvMsg%ElemToSide,       &
               RecvMsg%NodeCoords,       &          
               RecvMsg%ElemSideNodeID,   &          
               RecvMsg%SideToElem,       &          
               RecvMsg%BC,               &          
               RecvMsg%NativeElemID,     &          
               RecvMsg%PeriodicElemSide, &
               RecvMsg%ConcaveElemSide, &
               STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate RecvMsg%*, proc',iProc,999.)
  END IF
  IF (SendMsg%nNodes.GT.0) THEN
    DEALLOCATE(SendMsg%ElemToSide,       &
               SendMsg%NodeCoords,       &          
               SendMsg%ElemSideNodeID,   &          
               SendMsg%SideToElem,       &          
               SendMsg%BC,               &          
               SendMsg%NativeElemID,     &          
               SendMsg%PeriodicElemSide, &
               SendMsg%ConcaveElemSide,  &
               STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'Could not deallocate SendMsg%*, proc',iProc,999.)
  END IF
END SUBROUTINE ExchangeMPINeighborhoodGeometry

SUBROUTINE IdentifyNonImmediateMPINeighborhood(iProc,NodeIndex)
!===================================================================================================================================
!
! Searches for sides in the neighborhood of non-immediate MPI neighbors
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Particle_Vars,          ONLY : GEO
  USE MOD_part_MPI_Vars,          ONLY : PMPIVAR, halo_eps, FIBGMCellPadding
  USE MOD_mesh_vars,              ONLY : nNodes, nInnerSides, nBCSides, nElems, &
                                         ElemToSide, nSides
!  USE MOD_part_boundary_periodic, ONLY : periodicMove,getPeriodicVectorIndex,getNPVIndices
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)  :: iProc  ! MPI proc with which the local proc is to exchange boundary information
  INTEGER, INTENT(INOUT):: NodeIndex(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  TYPE tMPISideMessage
    REAL(KIND=8), POINTER :: Nodes(:,:) =>NULL()  ! nodes of boundary faces
    INTEGER(KIND=4)       :: nNodes               ! number of Nodes to send
  END TYPE
  TYPE(tMPISideMessage):: SendMsg
  TYPE(tMPISideMessage):: RecvMsg
  INTEGER              :: i
  INTEGER              :: iElem, jElem, iNode, iLocSide
  INTEGER              :: CellX ,CellY ,CellZ
  INTEGER              :: CellX2,CellY2,CellZ2
  INTEGER              :: PVindex
  INTEGER              :: nMoves
  INTEGER              :: ALLOCSTAT
  REAL                 :: xNode(1:3)
!===================================================================================================================================

  !--- Exchange MPI-Sides:
  !    Each proc receives all MPI sides that lie in FIBGM cells within an
  !    eps-distance of FIBGM cells containing any of my own MPI sides.
  !--- Go through each FIBGM cell and if there are MPI-Neighbors, search the
  !    nPaddingCells surrounding FIBGM cells for my own MPI sides.
  !--- Step1: find MPI side of myproc that are within padding distance to iProc

  NodeIndex(:)=0
  SendMsg%nNodes=0
  DO CellX=GEO%FIBGMimin,GEO%FIBGMimax
    DO CellY=GEO%FIBGMkmin,GEO%FIBGMkmax
      DO CellZ=GEO%FIBGMlmin,GEO%FIBGMlmax
        !IF (.NOT.ASSOCIATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) CYCLE
        IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%PaddingProcs)) CYCLE
        IF (GEO%FIBGM(CellX,CellY,CellZ)%PaddingProcs(1).LE.0) CYCLE       ! -1 if ???
        DO i = 2,GEO%FIBGM(CellX,CellY,CellZ)%PaddingProcs(1)+1
          IF (GEO%FIBGM(CellX,CellY,CellZ)%PaddingProcs(i).EQ.iProc) THEN
            DO iElem = 1, GEO%FIBGM(CellX,CellY,CellZ)%nElem
              jElem = GEO%FIBGM(CellX,CellY,CellZ)%Element(iElem)
              DO iLocSide=1,6
                IF (ElemToSide(E2S_SIDE_ID,iLocSide,jElem).GT.(nInnerSides+nBCSides)) THEN
                  DO iNode = 1,4
                    IF (NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,jElem)).EQ.0) THEN
                      SendMsg%nNodes=SendMsg%nNodes+1 
                      NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,jElem))=SendMsg%nNodes
                    END IF
                  END DO ! iNode
                END IF ! Side = MPI Side
              END DO ! Side
            END DO ! iElem
          END IF ! shapeProcs(i).EQ.iProc
        END DO ! i
      END DO ! CellZ
    END DO ! CellY
  END DO ! CellX

  !--- NOTE: IF SENDMSG%NNODES IS 0 AT THIS POINT, THEN I SHOULD BE ABLE TO RETURN HERE!!!
  !          This is not done yet because I'm not sure whether there are still inconsistencies in the code...

  !--- Debugging information
  LOGWRITE(*,*)'nNodes for iProc=',iProc,':',SendMsg%nNodes

  !--- Send number of MPI sides to MPI neighbor iProc and receive number of MPI
  !    sides from MPI neighbor iProc (immediate neighbor or not)
  IF (PMPIVAR%iProc.LT.iProc) THEN
    CALL MPI_SEND(SendMsg%nNodes,1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,IERROR)
    CALL MPI_RECV(RecvMsg%nNodes,1,MPI_INTEGER,iProc,1102,PMPIVAR%COMM,MPISTATUS,IERROR)
  ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
    CALL MPI_RECV(RecvMsg%nNodes,1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)
    CALL MPI_SEND(SendMsg%nNodes,1,MPI_INTEGER,iProc,1102,PMPIVAR%COMM,IERROR)
  END IF

  !--- Allocate Message
  IF (SendMsg%nNodes.GT.0) THEN
    ALLOCATE(SendMsg%Nodes(1:SendMsg%nNodes,3), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__,&
                           'Could not allocate SendMessage%Nodes ',SendMsg%nNodes,999.)
    END IF
    SendMsg%Nodes(:,:)=0.
  END IF
  IF (RecvMsg%nNodes.GT.0) THEN
    ALLOCATE(RecvMsg%Nodes(1:RecvMsg%nNodes,3), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__,&
                           'Could not allocate RecvMessage%Nodes ',RecvMsg%nNodes,999.)
    END IF
    RecvMsg%Nodes(:,:)=0.
  END IF

  !--- Send any (corner-) nodes from the MPI-sides to the MPI-neighbor iProc
  !    and receive iProc's (corner-) nodes in return
  !--- fill send buffers
  !--- Step 2: send myproc MPI-side-nodes to iProc
  DO iElem=1,nElems
    DO iLocSide=1,6
      IF (ElemToSide(E2S_SIDE_ID,iLocSide,iElem).GT.(nInnerSides+nBCSides)) THEN
        DO iNode = 1,4
          IF (NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,iElem)).NE.0) THEN
            xNode(1:3) = GEO%NodeCoords(1:3,GEO%ElemSideNodeID(iNode,iLocSide,iElem))
            SendMsg%Nodes(NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,iElem)),1:3) = xNode(1:3)
          END IF
        END DO ! iNode
      END IF ! Side = MPI Side
    END DO ! Side
  END DO

  ! CALL WriteDebugNodes(SendMsg%Nodes,SendMsg%nNodes,iProc)
  !--- send and receive data
  IF (PMPIVAR%iProc.LT.iProc) THEN
    IF (SendMsg%nNodes.GT.0) CALL MPI_SEND(SendMsg%Nodes,SendMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1103,PMPIVAR%COMM,IERROR)
    IF (RecvMsg%nNodes.GT.0) &
       CALL MPI_RECV(RecvMsg%Nodes,RecvMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1104,PMPIVAR%COMM,MPISTATUS,IERROR)
  ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
    IF (RecvMsg%nNodes.GT.0) &
       CALL MPI_RECV(RecvMsg%Nodes,RecvMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1103,PMPIVAR%COMM,MPISTATUS,IERROR)
    IF (SendMsg%nNodes.GT.0) CALL MPI_SEND(SendMsg%Nodes,SendMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1104,PMPIVAR%COMM,IERROR)
  END IF

  !--- For each node, identify the FIBGM cell(s) in which
  !    the node resides and search the surrounding nPaddingCells for neighboring
  !    elements.
  !--- For idiots: iProc tells me which nodes are on his MPI bound,
  !---             and I check which nodes are within eps range.
  NodeIndex(:)=0
  IF (RecvMsg%nNodes.GT.0) THEN
    CALL CheckMPINeighborhoodByFIBGM(RecvMsg%Nodes,RecvMsg%nNodes,NodeIndex)
  END IF

  !--- Deallocate Messages
  IF (SendMsg%nNodes.GT.0) THEN
    DEALLOCATE(SendMsg%Nodes, STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__,&
                           'Could not deallocate SendMessage%Nodes proc ',iProc,999.)
    END IF
  END IF
  IF (RecvMsg%nNodes.GT.0) THEN
    DEALLOCATE(RecvMsg%Nodes, STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__,&
                           'Could not deallocate RecvMessage%Nodes proc ',iProc,999.)
    END IF
  END IF
END SUBROUTINE IdentifyNonImmediateMPINeighborhood

SUBROUTINE CheckMPINeighborhoodByFIBGM(NodeData,nExternalNodes,NodeIndex)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars,          ONLY : GEO
  USE MOD_part_MPI_Vars,          ONLY : PMPIVAR, PMPIExchange, halo_eps, FIBGMCellPadding,&
                                         NbrOfCases,casematrix
  USE MOD_mesh_vars,              ONLY : nNodes
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,    INTENT(IN)      :: NodeData(1:nExternalNodes,1:3)
  INTEGER, INTENT(IN)      :: nExternalNodes
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(INOUT)   :: NodeIndex(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                  :: CellX, CellY, CellZ
  INTEGER                  :: CellX2,CellY2,CellZ2
  INTEGER                  :: iElem, iNode, jNode
  INTEGER                  :: nNbNodes, jElem, ind, ind2 , iCase
  REAL                     :: NodeX(1:3), Vec1(1:3), Vec2(1:3), Vec3(1:3)
!===================================================================================================================================

  !--- For each side built from 4 nodes, identify the FIBGM cell(s) in which
  !    the side resides and search the surrounding nPaddingCells for neighboring
  !    elements.
  !--- for idiots: get nodes of myProc that are within eps distance to MPI-bound 
  !                of iProc
  NodeIndex(:)=0
  nNbNodes=0
  DO iNode=1,nExternalNodes
    NodeX(1:3) = NodeData(iNode,1:3)
    CellX = INT((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
    CellY = INT((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
    CellZ = INT((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
    DO CellX2 = CellX-FIBGMCellPadding(1),CellX+FIBGMCellPadding(1)
      DO CellY2 = CellY-FIBGMCellPadding(2),CellY+FIBGMCellPadding(2)
        DO CellZ2 = CellZ-FIBGMCellPadding(3),CellZ+FIBGMCellPadding(3)
          IF ((CellX2.GT.GEO%FIBGMimax).OR.(CellX2.LT.GEO%FIBGMimin) .OR. &
              (CellY2.GT.GEO%FIBGMkmax).OR.(CellY2.LT.GEO%FIBGMkmin) .OR. &
              (CellZ2.GT.GEO%FIBGMlmax).OR.(CellZ2.LT.GEO%FIBGMlmin) ) CYCLE
          DO iElem = 1, GEO%FIBGM(CellX2,CellY2,CellZ2)%nElem
            jElem = GEO%FIBGM(CellX2,CellY2,CellZ2)%Element(iElem)
            DO jNode = 1,8
              IF ( (NodeIndex(GEO%ElemToNodeID(jNode,jElem)).EQ.0)                 .AND. &
                   (DOT_PRODUCT(GEO%NodeCoords(1:3,GEO%ElemToNodeID(jNode,jElem))-NodeX, &
                                GEO%NodeCoords(1:3,GEO%ElemToNodeID(jNode,jElem))-NodeX).LE.halo_eps*halo_eps) ) THEN
                nNbNodes=nNbNodes+1
                NodeIndex(GEO%ElemToNodeID(jNode,jElem))=nNbNodes
              END IF
            END DO ! jNode
          END DO ! iElem
        END DO ! CellZ
      END DO ! CellY
    END DO ! CellX2
  END DO ! iNode
  !--- if there are periodic boundaries, they need to be taken into account as well:
  IF (GEO%nPeriodicVectors.GT.0) THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (GEO%nPeriodicVectors.EQ.1) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    END IF
    IF (GEO%nPeriodicVectors.EQ.2) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    END IF
    IF (GEO%nPeriodicVectors.EQ.3) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
    END IF
    !--- check nodes shifted by periodic vectors, add to NodeIndex if match
    DO iNode=1,nExternalNodes
      DO iCase = 1, NbrOfCases
        IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
            (casematrix(iCase,2).EQ.0) .AND. &
            (casematrix(iCase,3).EQ.0)) CYCLE
        NodeX(1:3) = NodeData(iNode,1:3)           + &
                     casematrix(iCase,1)*Vec1(1:3) + &
                     casematrix(iCase,2)*Vec2(1:3) + &
                     casematrix(iCase,3)*Vec3(1:3) 
        CellX = INT((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
        CellY = INT((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
        CellZ = INT((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
        DO CellX2 = CellX-FIBGMCellPadding(1),CellX+FIBGMCellPadding(1)
          DO CellY2 = CellY-FIBGMCellPadding(2),CellY+FIBGMCellPadding(2)
            DO CellZ2 = CellZ-FIBGMCellPadding(3),CellZ+FIBGMCellPadding(3)
              IF ((CellX2.GT.GEO%FIBGMimax).OR.(CellX2.LT.GEO%FIBGMimin) .OR. &
                   (CellY2.GT.GEO%FIBGMkmax).OR.(CellY2.LT.GEO%FIBGMkmin) .OR. &
                   (CellZ2.GT.GEO%FIBGMlmax).OR.(CellZ2.LT.GEO%FIBGMlmin) ) CYCLE
                DO iElem = 1, GEO%FIBGM(CellX2,CellY2,CellZ2)%nElem
                  jElem = GEO%FIBGM(CellX2,CellY2,CellZ2)%Element(iElem)
                  DO jNode = 1,8
                     IF ( (NodeIndex(GEO%ElemToNodeID(jNode,jElem)).EQ.0)                 .AND. &
                          (DOT_PRODUCT(GEO%NodeCoords(1:3,GEO%ElemToNodeID(jNode,jElem))-NodeX, &
                                       GEO%NodeCoords(1:3,GEO%ElemToNodeID(jNode,jElem))-NodeX).LE.halo_eps*halo_eps) ) THEN
                      nNbNodes=nNbNodes+1
                      NodeIndex(GEO%ElemToNodeID(jNode,jElem))=nNbNodes
                    END IF
                  END DO ! jNode
                END DO ! iElem
            END DO ! CellZ
          END DO ! CellY
        END DO ! CellX2
      END DO ! iCase
    END DO ! iNode
  END IF  ! nperiodicvectors>0
END SUBROUTINE CheckMPINeighborhoodByFIBGM

SUBROUTINE Finalize()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_part_MPI_Vars, ONLY : PMPIVAR, PMPIExchange
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                :: ALLOCSTAT
  INTEGER                :: iProc
!===================================================================================================================================

  !--- deallocate message buffers and request buffers
  DO iProc=0,PMPIVAR%nProcs-1
    IF (ASSOCIATED(PMPIExchange%send_message)) THEN
      IF (ASSOCIATED(PMPIExchange%send_message(iProc)%content)) THEN
        DEALLOCATE( PMPIExchange%send_message(iProc)%content,STAT=ALLOCSTAT )
        IF (ALLOCSTAT .NE. 0) THEN
          WRITE(*,*)'ERROR in MPI Particle Finalize: cannot deallocate send_message(',iProc,')!'
        END IF
      END IF
    END IF
    !IF (ASSOCIATED(PMPIExchange%recv_message)) THEN
    !  IF (ASSOCIATED(PMPIExchange%recv_message(iProc)%content)) THEN
    !    DEALLOCATE( PMPIExchange%recv_message(iProc)%content,STAT=ALLOCSTAT )
    !    IF (ALLOCSTAT .NE. 0) THEN
    !      WRITE(*,*)'ERROR in MPI Particle Finalize: cannot deallocate recv_message(',iProc,')!'
    !    END IF
    !  END IF
    !END IF
  END DO
  IF (ASSOCIATED(PMPIExchange%send_message)) DEALLOCATE(PMPIExchange%send_message)
  IF (ASSOCIATED(PMPIExchange%send_request))DEALLOCATE(PMPIExchange%send_request)
  IF (ASSOCIATED(PMPIExchange%nbrOfSendParticles))DEALLOCATE(PMPIExchange%nbrOfSendParticles)
  !SDEALLOCATE(PMPIExchange%recv_message)
  IF (ASSOCIATED(PMPIExchange%recv_request))DEALLOCATE(PMPIExchange%recv_request)
  !IF (ASSOCIATED(PMPIExchange%nbrOfRecvParticles))DEALLOCATE(PMPIExchange%nbrOfRecvParticles)
END SUBROUTINE Finalize

#endif

END MODULE MOD_part_boundary_mpi
