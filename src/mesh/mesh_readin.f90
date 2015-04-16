#include "boltzplatz.h"

MODULE MOD_Mesh_ReadIn
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER    :: ElemInfoSize=6        !number of entry in each line of ElemInfo
INTEGER,PARAMETER    :: ELEM_Type=1           !entry position, 
INTEGER,PARAMETER    :: ELEM_Zone=2           
INTEGER,PARAMETER    :: ELEM_FirstSideInd=3
INTEGER,PARAMETER    :: ELEM_LastSideInd=4
INTEGER,PARAMETER    :: ELEM_FirstNodeInd=5
INTEGER,PARAMETER    :: ELEM_LastNodeInd=6

INTEGER,PARAMETER    :: SideInfoSize=4
INTEGER,PARAMETER    :: SIDE_Type=1           !entry position
INTEGER,PARAMETER    :: SIDE_ID=2
INTEGER,PARAMETER    :: SIDE_nbElemID=3
INTEGER,PARAMETER    :: SIDE_BCID=4

INTEGER,ALLOCATABLE  :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE     :: ElemWeight(:)
REAL,ALLOCATABLE     :: NodeCoords(:,:)
INTEGER,ALLOCATABLE  :: NodeMap(:)

INTEGER              :: BoundaryOrder_mesh
INTEGER              :: nNodeIDs,nSideIDs
INTEGER              :: offsetSideID,offsetNodeID
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMesh
  MODULE PROCEDURE ReadMesh
END INTERFACE

INTERFACE Qsort1Int
  MODULE PROCEDURE Qsort1Int
END INTERFACE

INTERFACE INVMAP
  MODULE PROCEDURE INVMAP
END INTERFACE

PUBLIC::ReadMesh,Qsort1Int,INVMAP
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
! Read boundary conditions from data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:BoundaryName,BoundaryType,nBCs,nUserBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE           :: BCMapping(:),BCType(:,:)
CHARACTER(LEN=255), ALLOCATABLE:: BCNames(:)
INTEGER                        :: iBC,iUserBC
INTEGER                        :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================
! Read boundary names from data file
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
nBCs=INT(HSize(1),4)
DEALLOCATE(HSize)
ALLOCATE(BCNames(nBCs))
ALLOCATE(BCMapping(nBCs))
CALL ReadArray('BCNames',1,(/nBCs/),Offset,1,StrArray=BCNames)  ! Type is a dummy type only
! User may have redefined boundaries in the ini file. So we have to create mappings for the boundaries.
BCMapping=0
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    DO iUserBC=1,nUserBCs
      IF(INDEX(TRIM(BCNames(iBC)),TRIM(BoundaryName(iUserBC))) .NE.0) BCMapping(iBC)=iUserBC
    END DO
  END DO
END IF

! Read boundary types from data file
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF(HSize(1).NE.nBCs) STOP 'Problem in readBC'
DEALLOCATE(HSize)
ALLOCATE(BCType(nBCs,4))
offset=0
CALL ReadArray('BCType',2,(/nBCs,4/),Offset,1,IntegerArray=BCType)
! Now apply boundary mappings
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    IF(BCMapping(iBC) .NE. 0)THEN
      SWRITE(Unit_StdOut,'(A,A)')    ' |     Boundary in HDF file found |  ',TRIM(BCNames(iBC))
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                            was | ',BCType(iBC,1),BCType(iBC,3)
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                      is set to | ',BoundaryType(BCMapping(iBC),:)
      BCType(iBC,1) = BoundaryType(BCMapping(iBC),BC_TYPE)
      BCType(iBC,3) = BoundaryType(BCMapping(iBC),BC_STATE)
    END IF
  END DO
END IF
IF(ALLOCATED(BoundaryName)) DEALLOCATE(BoundaryName)
IF(ALLOCATED(BoundaryType)) DEALLOCATE(BoundaryType)
ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,2))
BoundaryName = BCNames
BoundaryType(:,BC_TYPE)  = BCType(:,1)
BoundaryType(:,BC_STATE) = BCType(:,3)
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(Unit_StdOut,'(A,A16,A20,A10,A10)')'BOUNDARY CONDITIONS','|','Name','Type','State'
DO iBC=1,nBCs
  SWRITE(*,'(A,A33,A20,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:) 
END DO
SWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs


SUBROUTINE ReadMesh(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,          ONLY:NGeo
USE MOD_Mesh_Vars,          ONLY:offsetElem,nElems,nGlobalElems,nNodes
USE MOD_Mesh_Vars,          ONLY:nSides,nInnerSides,nBCSides,nMPISides
USE MOD_Mesh_Vars,          ONLY:useCurveds
USE MOD_Mesh_Vars,          ONLY:BoundaryType
USE MOD_Mesh_Vars,          ONLY:MeshInitIsDone
USE MOD_Mesh_Vars,          ONLY:Elems,Nodes
USE MOD_Mesh_Vars,          ONLY:aElem,aSide,bSide
USE MOD_Mesh_Vars,          ONLY:GETNEWELEM,GETNEWSIDE,createSides
#ifdef MPI
USE MOD_Mesh_Vars,          ONLY:ParticleMPIWeight
USE MOD_MPI_Vars,           ONLY:offsetElemMPI,nMPISides_Proc,nNbProcs,NbProc
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Restart_Vars,       ONLY:DoRestart,RestartFile
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: BCindex
INTEGER                        :: iElem,ElemID
INTEGER                        :: iNode,jNode,iNodeP,NodeID,SideID
INTEGER                        :: iLocSide,jLocSide
INTEGER                        :: iSide
INTEGER                        :: FirstNodeInd,LastNodeInd,FirstSideInd,LastSideInd,FirstElemInd,LastElemInd
LOGICAL                        :: oriented
INTEGER                        :: nPeriodicSides 
INTEGER                        :: ReduceData(7)
#ifdef MPI
INTEGER                        :: ReduceData_glob(7)
INTEGER                        :: iNbProc, locnPart
INTEGER                        :: iProc, curiElem
INTEGER,ALLOCATABLE            :: MPISideCount(:)
INTEGER,ALLOCATABLE            :: PartInt(:,:)
INTEGER,PARAMETER              :: ELEM_FirstPartInd=1
INTEGER,PARAMETER              :: ELEM_LastPartInd=2
REAL,ALLOCATABLE               :: ElemWeight(:)
REAL                           :: SumWeight, CurWeight
#endif
LOGICAL                        :: fileExists
LOGICAL                        :: doConnection
!===================================================================================================================================
IF(MESHInitIsDone) RETURN
IF(MPIRoot)THEN
  INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
  IF(.NOT.FileExists)  CALL abort(__STAMP__, &
          'readMesh from data file "'//TRIM(FileString)//'" does not exist',999,999.)
END IF

SWRITE(UNIT_stdOut,'(A)')'READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')

! Open data file
#ifdef MPI
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)
#else
CALL OpenDataFile(FileString,create=.FALSE.)
#endif

CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
nGlobalElems=INT(HSize(1),4) !global number of elements
DEALLOCATE(HSize)
#ifdef MPI
!simple partition: nGlobalelems/nprocs, do this on proc 0
IF(ALLOCATED(offsetElemMPI))DEALLOCATE(offsetElemMPI)
ALLOCATE(offsetElemMPI(0:nProcessors))
offsetElemMPI=0

IF (DoRestart) THEN 
  CALL CloseDataFile() 
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.)
  ALLOCATE(PartInt(1:nGlobalElems,2))
  PartInt(:,:)=0
  CALL ReadArray('PartInt',2,(/nGlobalElems,2/),0,1,IntegerArray=PartInt)
  CALL CloseDataFile() 
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)
  
  SumWeight = 0.0
  CurWeight = 0.0
  ParticleMPIWeight = GETREAL('Particles-MPIWeight','0.02')
  IF (ParticleMPIWeight.LT.0) THEN
    WRITE(*,*) "ERROR: Particle weight can't be negative!"
    STOP
  END IF
  ! load balancing for particles
  ! read in particle data

  ALLOCATE(ElemWeight(1:nGlobalElems))
  DO iElem = 1, nGlobalElems
    locnPart=PartInt(iElem,ELEM_LastPartInd)-PartInt(iElem,ELEM_FirstPartInd)
    ElemWeight(iElem) = locnPart*ParticleMPIWeight + 1.0
    SumWeight = SumWeight + ElemWeight(iElem)
  END DO

  SumWeight=SumWeight/nProcessors 
  curiElem = 1
  DO iProc=0, nProcessors-1
    offsetElemMPI(iProc)=curiElem - 1 
    DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc  
      CurWeight=CurWeight+ElemWeight(iElem)  
      IF (CurWeight.GE.SumWeight*(iProc+1)) THEN
        curiElem = iElem + 1 
        EXIT
      END IF
    END DO   
  END DO
  
ELSE
  nElems=nGlobalElems/nProcessors
  iElem=nGlobalElems-nElems*nProcessors
  DO iProc=0,nProcessors-1
    offsetElemMPI(iProc)=nElems*iProc+MIN(iProc,iElem)
  END DO
END IF
offsetElemMPI(nProcessors)=nGlobalElems
!local nElems and offset
nElems=offsetElemMPI(myRank+1)-offsetElemMPI(myRank)
offsetElem=offsetElemMPI(myRank)
LOGWRITE(*,*)'offset,nElems',offsetElem,nElems
#else /* MPI */
nElems=nGlobalElems   !local number of Elements 
offsetElem=0          ! offset is the index of first entry, hdf5 array starts at 0-.GT. -1 
#endif /* MPI */


CALL readBCs()

CALL ReadAttribute(File_ID,'BoundaryOrder',1,IntegerScalar=BoundaryOrder_mesh)
IF(useCurveds) THEN
  IF(NGeo+1.NE.BoundaryOrder_mesh) THEN
    CALL abort(__STAMP__, &
          ' NGeo does not correspond to boundary order in Meshfile, set NGeo to:',BoundaryOrder_mesh-1,999.)
    NGeo=BoundaryOrder_mesh-1
  END IF
END IF
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
FirstElemInd=offsetElem+1
LastElemInd=offsetElem+nElems
ALLOCATE(ElemInfo(FirstElemInd:LastElemInd,ElemInfoSize))
CALL ReadArray('ElemInfo',2,(/nElems,ElemInfoSize/),offsetElem,1,IntegerArray=ElemInfo)


ALLOCATE(Elems(FirstElemInd:LastElemInd))

DO iElem=FirstElemInd,LastElemInd
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  Elems(iElem)%ep=>GETNEWELEM()
  aElem=>Elems(iElem)%ep
  aElem%Ind    = iElem
  aElem%Type   = ElemInfo(iElem,ELEM_Type)
  aElem%Zone   = ElemInfo(iElem,ELEM_Zone)
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from data file 
offsetNodeID=ElemInfo(FirstElemInd,ELEM_FirstNodeInd) ! hdf5 array starts at 0-> -1
nNodeIDs=ElemInfo(LastElemInd,ELEM_LastNodeInd)-ElemInfo(FirstElemind,ELEM_FirstNodeInd)
FirstNodeInd=offsetNodeID+1
LastNodeInd=offsetNodeID+nNodeIDs
ALLOCATE(NodeInfo(FirstNodeInd:LastNodeInd))
CALL ReadArray('NodeInfo',1,(/nNodeIDs/),offsetNodeID,1,IntegerArray=NodeInfo)

!WRITE(*,*)'DEBUG, NodeInfo'
!DO i=FirstNodeInd,LastNodeInd
!  WRITE(*,*)i,':',NodeInfo(i)
!END DO
 
CALL GetNodeMap() !get nNodes and NodeMap from NodeInfo array
LOGWRITE(*,*)'MIN,MAX,SIZE of NodeMap',MINVAL(NodeMap),MAXVAL(NodeMap),SIZE(NodeMap,1)

ALLOCATE(Nodes(1:nNodes)) ! pointer list, entry is known by INVMAP(i,nNodes,NodeMap)
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
!assign nodes 
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  DO jNode=1,8
    iNode=iNode+1
    NodeID=ABS(NodeInfo(iNode))     !global, unique NodeID
    iNodeP=INVMAP(NodeID,nNodes,NodeMap)  ! index in local Nodes pointer array
    IF(iNodeP.LE.0) STOP 'Problem in INVMAP' 
    IF(.NOT.ASSOCIATED(Nodes(iNodeP)%np))THEN
      ALLOCATE(Nodes(iNodeP)%np) 
      Nodes(iNodeP)%np%ind=NodeID 
    END IF
    aElem%Node(jNode)%np=>Nodes(iNodeP)%np
  END DO
  CALL createSides(aElem)
  aElem%nCurvedNodes=0
  IF(useCurveds) THEN
    aElem%nCurvedNodes= ElemInfo(iElem,ELEM_LastNodeInd) - ElemInfo(iElem,ELEM_FirstNodeInd) - 14 ! corner + oriented nodes
    IF((aElem%nCurvedNodes.GT.0).OR.(NGeo.GT.1))THEN
      ALLOCATE(aElem%CurvedNode(aElem%nCurvedNodes))
      DO jNode=1,aElem%nCurvedNodes
        iNode=iNode+1
        NodeID=NodeInfo(iNode) !first oriented corner node
        iNodeP=INVMAP(NodeID,nNodes,NodeMap)  ! index in local Nodes pointer array
        IF(iNodeP.LE.0) STOP 'Problem in INVMAP' 
        IF(.NOT.ASSOCIATED(Nodes(iNodeP)%np))THEN
          ALLOCATE(Nodes(iNodeP)%np)
          Nodes(iNodeP)%np%ind=NodeID 
        END IF
        aElem%CurvedNode(jNode)%np=>Nodes(iNodeP)%np
      END DO !jNode=1,nCurvedNodes
    END IF
  END IF
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------


offsetSideID=ElemInfo(FirstElemInd,ELEM_FirstSideInd) ! hdf5 array starts at 0-> -1  
nSideIDs=ElemInfo(LastElemInd,ELEM_LastSideInd)-ElemInfo(FirstElemind,ELEM_FirstSideInd)
!read local SideInfo from data file 
FirstSideInd=offsetSideID+1
LastSideInd=offsetSideID+nSideIDs
ALLOCATE(SideInfo(FirstSideInd:LastSideInd,SideInfoSize))
CALL ReadArray('SideInfo',2,(/nSideIDs,SideinfoSize/),offsetSideID,1,IntegerArray=SideInfo)

!WRITE(*,*)'DEBUG, SideInfo'
!DO i=FirstSideInd,LastSideInd
!  WRITE(*,*)i,':',SideInfo(i,:)
!END DO


DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  !iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  !iNode=iNode+8
  iNode=ElemInfo(iElem,ELEM_LastNodeInd) !first index -1 in NodeInfo
  iNode=iNode-6
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  !build up sides of the element using element Nodes and CGNS standard
  ! assign flip
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    aSide%Elem=>aElem
    oriented=(Sideinfo(iSide,SIDE_ID).GT.0)
    
    aSide%Ind=ABS(SideInfo(iSide,SIDE_ID))
    iNode=iNode+1
    NodeID=NodeInfo(iNode) !first oriented corner node
    IF(oriented)THEN !oriented side
      aSide%flip=0
    ELSE !not oriented
      DO jNode=1,4
        IF(aSide%Node(jNode)%np%ind.EQ.ABS(NodeID)) EXIT
      END DO
      IF(jNode.GT.4) STOP 'NodeID doesnt belong to side'
      aSide%flip=jNode
    END IF
  END DO !i=1,locnSides
END DO !iElem

 
! build up side connection 
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    sideID  = ABS(SideInfo(iSide,SIDE_ID))
    elemID  = SideInfo(iSide,SIDE_nbElemID)
    BCindex = SideInfo(iSide,SIDE_BCID)
    doConnection=.TRUE. ! for periodic sides if BC is reassigned as non periodic
    IF(BCindex.NE.0)THEN !BC
      aSide%BCindex = BCindex
      IF(BoundaryType(aSide%BCindex,BC_TYPE).NE.1)THEN ! Reassignement from periodic to non-periodic
        doConnection=.FALSE.
        aSide%flip  =0
      END IF
    ELSE
      aSide%BCindex = 0
    END IF
    IF(.NOT.ASSOCIATED(aSide%connection))THEN
      IF((elemID.NE.0).AND.doConnection)THEN !connection 
        IF((elemID.LE.LastElemInd).AND.(elemID.GE.FirstElemInd))THEN !local connection
          DO jLocSide=1,6
            bSide=>Elems(elemID)%ep%Side(jLocSide)%sp
            IF(bSide%ind.EQ.aSide%ind)THEN
              aSide%connection=>bSide
              bSide%connection=>aSide
              EXIT
            END IF
          END DO
        ELSE !MPI connection
#ifdef MPI
          aSide%connection=>GETNEWSIDE()            
          aSide%connection%flip=aSide%flip
          aSide%connection%Elem=>GETNEWELEM()
          aSide%NbProc = ELEMIPROC(elemID)
#else
        CALL abort(__STAMP__, &
          ' elemID of neighbor not in global Elem list ',999,999.)
#endif
        END IF
      END IF
    END IF !connection associated
  END DO !iLocSide 
END DO !iElem

DEALLOCATE(ElemInfo,SideInfo,NodeInfo)

! get physical coordinates

ALLOCATE(NodeCoords(nNodes,3))

CALL ReadCoords(File_ID,'NodeCoords',(/nNodes,3/),NodeMap,NodeCoords)
!assign to pointer
DO iNode=1,nNodes
  IF(ASSOCIATED(Nodes(iNode)%np))Nodes(iNode)%np%x=NodeCoords(iNode,:)
END DO
DEALLOCATE(NodeCoords)
 
#ifdef MPI
!  DEALLOCATE(offsetElemMPI) ! DO NOT DEALLOCATE, will be used by readData/writedata (io_restart.f90)!
#endif /*MPI*/

DEALLOCATE(NodeMap)
CALL CloseDataFile() 

! COUNT SIDES

 
nBCSides=0
nSides=0
nPeriodicSides=0
nMPISides=0
#ifdef MPI
ALLOCATE(MPISideCount(0:nProcessors-1))
MPISideCount=0
#endif
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aElem%Side(iLocSide)%sp%tmp=0
  END DO
END DO
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%tmp.EQ.0)THEN
      nSides=nSides+1
      aSide%tmp=-1 !used as marker
      IF(ASSOCIATED(aSide%connection)) aSide%connection%tmp=-1
      IF(aSide%BCindex.NE.0)THEN !side is BC or periodic side
        IF(ASSOCIATED(aSide%connection))THEN
          nPeriodicSides=nPeriodicSides+1
        ELSE
          nBCSides=nBCSides+1
        END IF
      END IF
#ifdef MPI
      IF(aSide%NbProc.NE.-1) THEN
        nMPISides=nMPISides+1
        MPISideCount(aSide%NbProc)=MPISideCount(aSide%NbProc)+1
      END IF
#endif
    END IF
  END DO
END DO
nInnerSides=nSides-nBCSides-nMPISides !periodic side count to inner side!!!

LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nSides:',nSides
LOGWRITE(*,'(A22,I8)')'nBCSides:',nBCSides
LOGWRITE(*,'(A22,I8)')'nInnerSides:',nInnerSides
LOGWRITE(*,'(A22,I8)')'nMPISides:',nMPISides
LOGWRITE(*,*)'-------------------------------------------------------'
 !now MPI sides
#ifdef MPI
nNBProcs=0
DO iProc=0,nProcessors-1
  IF(iProc.EQ.myRank) CYCLE
  IF(MPISideCount(iProc).GT.0) nNBProcs=nNbProcs+1
END DO
IF(nNbProcs.EQ.0)THEN !MPI + 1Proc case !
  ALLOCATE(NbProc(1),nMPISides_Proc(1))
  nNbProcs=1
  NbProc=0
  nMPISides_Proc=0
ELSE
  ALLOCATE(NbProc(nNbProcs),nMPISides_Proc(1:nNbProcs))
  iNbProc=0
  DO iProc=0,nProcessors-1
    IF(iProc.EQ.myRank) CYCLE
    IF(MPISideCount(iProc).GT.0) THEN
      iNbProc=iNbProc+1
      NbProc(iNbProc)=iProc
      ! compute number of MPISides per neighbor proc and divide by two
      nMPISides_Proc(iNBProc)=MPISideCount(iProc)
    END IF
  END DO
END IF
DEALLOCATE(MPISideCount)

#endif /*MPI*/

ReduceData(1)=nElems
ReduceData(2)=nSides
ReduceData(3)=nNodes
ReduceData(4)=nInnerSides
ReduceData(5)=nPeriodicSides
ReduceData(6)=nBCSides
ReduceData(7)=nMPISides

#ifdef MPI
CALL MPI_REDUCE(ReduceData,ReduceData_glob,7,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
ReduceData=ReduceData_glob
#endif /*MPI*/

IF(MPIRoot)THEN
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nElems | ',ReduceData(1) !nElems
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides | ',ReduceData(2) !nSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nNodes | ',ReduceData(3) !nNodes
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nInnerSides,not periodic | ',ReduceData(4)-ReduceData(5) !nInnerSides-nPeriodicSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','                periodic | ',ReduceData(5) !nPeriodicSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nBCSides | ',ReduceData(6) !nBCSides
  WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nMPISides | ',ReduceData(7)/2 !nMPISides
  WRITE(UNIT_stdOut,'(132("."))')
END IF

SWRITE(UNIT_stdOut,'(132("."))')
END SUBROUTINE ReadMesh



SUBROUTINE GetNodeMap()
!===================================================================================================================================
! take NodeInfo array, sort it, eliminate mulitple IDs and return the Mapping 1->NodeID1, 2->NodeID2, ... 
! this is useful if the NodeID list of the mesh are not contiguous, essentially occuring when using domain decomposition (MPI)
!===================================================================================================================================
! MODULES
USE MOD_mesh_vars,ONLY:nNodes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: temp(nNodeIDs+1),i,nullpos
!===================================================================================================================================
temp(1)=0
temp(2:nNodeIDs+1)=NodeInfo
!sort
CALL Qsort1Int(temp)
nullpos=INVMAP(0,nNodeIDs+1,temp)
!count unique entries
nNodes=1
DO i=nullpos+2,nNodeIDs+1
  IF(temp(i).NE.temp(i-1)) nNodes = nNodes+1
END DO
!associate unique entries
ALLOCATE(NodeMap(nNodes))
nNodes=1
NodeMap(1)=temp(nullpos+1)
DO i=nullpos+2,nNodeIDs+1
  IF(temp(i).NE.temp(i-1)) THEN
    nNodes = nNodes+1
    NodeMap(nNodes)=temp(i)
  END IF
END DO
END SUBROUTINE GetNodeMap 


FUNCTION INVMAP(ID,nIDs,ArrID)
!===================================================================================================================================
! find the inverse Mapping p.e. NodeID-> entry in NodeMap (a sorted array of unique NodeIDs), using bisection 
! if Index is not in the range, -1 will be returned, if it is in the range, but is not found, 0 will be returned!!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ID            ! ID to search for
INTEGER, INTENT(IN)                :: nIDs          ! size of ArrID
INTEGER, INTENT(IN)                :: ArrID(nIDs)   ! 1D array of IDs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: INVMAP               ! index of ID in NodeMap array
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid
!===================================================================================================================================
INVMAP=0
maxSteps=INT(LOG(REAL(nIDs))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=1
up=nIDs
IF((ID.LT.ArrID(low)).OR.(ID.GT.ArrID(up))) THEN
  !WRITE(*,*)'WARNING, Node Index Not in local range -> set to -1'
  INVMAP=-1  ! not in the range!
  RETURN
END IF 
IF(ID.EQ.ArrID(low))THEN
  INVMAP=low
ELSEIF(ID.EQ.ArrID(up))THEN
  INVMAP=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF(ID .EQ. ArrID(mid))THEN
      INVMAP=mid                     !index found!
      EXIT
    ELSEIF(ID .GT. ArrID(mid))THEN ! seek in upper half
      low=mid
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION INVMAP 



#ifdef MPI
FUNCTION ELEMIPROC(ElemID)
!===================================================================================================================================
! find the proc on which Element with elemID lies, use offsetElemMPI array and bisection 
!===================================================================================================================================
! MODULES
USE MOD_Globals,   ONLY:nProcessors
USE MOD_MPI_vars,  ONLY:offsetElemMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ElemID            ! NodeID to search for
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: ELEMIPROC         
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid
!===================================================================================================================================
ELEMIPROC=0
maxSteps=INT(LOG(REAL(nProcessors))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=0
up=nProcessors-1
IF((ElemID.GT.offsetElemMPI(low)).AND.(ElemID.LE.offsetElemMPI(low+1)))THEN
  ELEMIPROC=low
ELSEIF((ElemID.GT.offsetElemMPI(up)).AND.(ElemID.LE.offsetElemMPI(up+1)))THEN
  ELEMIPROC=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF((ElemID.GT.offsetElemMPI(mid)).AND.(ElemID.LE.offsetElemMPI(mid+1)))THEN
      ELEMIPROC=mid                     !index found!
      EXIT
    ELSEIF(ElemID .GT. offsetElemMPI(mid+1))THEN ! seek in upper half
      low=mid+1
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION ELEMIPROC 
#endif /* MPI */



RECURSIVE SUBROUTINE Qsort1Int(A)
!===================================================================================================================================
! QuickSort for integer array A
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: A(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: marker
!===================================================================================================================================
IF(SIZE(A).GT.1) THEN
  CALL Partition1Int(A,marker)
  CALL Qsort1Int(A(:marker-1))
  CALL Qsort1Int(A(marker:))
END IF
RETURN
END SUBROUTINE Qsort1Int       



SUBROUTINE Partition1Int(A,marker)
!===================================================================================================================================
! Neeeded by QuickSort
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: A(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)              :: marker
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,j
INTEGER                          :: temp,x
!===================================================================================================================================
x= A(1)
i= 0
j= SIZE(A)+1
DO
  j=j-1
  DO
    IF(A(j).LE.x) EXIT
    j=j-1
  END DO
  i=i+1
  DO
    IF(A(i).GE.x) EXIT
    i=i+1
  END DO
  IF(i.LT.j)THEN
    ! exchange A(i) and A(j)
    temp=A(i)
    A(i)=A(j)
    A(j)=temp
  ELSEIF(i.EQ.j)THEN
    marker=i+1
    RETURN
  ELSE
    marker=i
    RETURN
  ENDIF                      
END DO                        
RETURN                        
END SUBROUTINE Partition1Int     

END MODULE MOD_Mesh_ReadIn
