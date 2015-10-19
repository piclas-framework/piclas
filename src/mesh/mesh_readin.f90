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

INTEGER,PARAMETER    :: SideInfoSize=5
INTEGER,PARAMETER    :: SIDE_Type=1           !entry position
INTEGER,PARAMETER    :: SIDE_ID=2
INTEGER,PARAMETER    :: SIDE_nbElemID=3
INTEGER,PARAMETER    :: SIDE_Flip=4
INTEGER,PARAMETER    :: SIDE_BCID=5

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMesh
  MODULE PROCEDURE ReadMesh
END INTERFACE

PUBLIC::ReadMesh
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
! Read boundary conditions from data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY:BoundaryName,BoundaryType,nBCs,nUserBCs
USE MOD_ReadInTools,ONLY:GETINTARRAY,CNTSTR,GETSTR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE            :: UserBCFound(:)
CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)
INTEGER, ALLOCATABLE           :: BCMapping(:),BCType(:,:)
INTEGER                        :: iBC,iUserBC
INTEGER                        :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================
! read in boundary conditions from ini file, will overwrite BCs from meshfile!
nUserBCs = CNTSTR('BoundaryName','0')
IF(nUserBCs.GT.0)THEN
  ALLOCATE(BoundaryName(1:nUserBCs))
  ALLOCATE(BoundaryType(1:nUserBCs,2))
  DO iBC=1,nUserBCs
    BoundaryName(iBC)   = GETSTR('BoundaryName')
    BoundaryType(iBC,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
  END DO
END IF !nUserBCs>0

! Read boundary names from data file
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
nBCs=INT(HSize(1),4)
DEALLOCATE(HSize)
ALLOCATE(BCNames(nBCs))
ALLOCATE(BCMapping(nBCs))
ALLOCATE(UserBCFound(nUserBCs))
CALL ReadArray('BCNames',1,(/nBCs/),Offset,1,StrArray=BCNames)  ! Type is a dummy type only
! User may have redefined boundaries in the ini file. So we have to create mappings for the boundaries.
BCMapping=0
UserBCFound=.FALSE.
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    DO iUserBC=1,nUserBCs
      IF(INDEX(TRIM(BCNames(iBC)),TRIM(BoundaryName(iUserBC))).NE.0)THEN
        BCMapping(iBC)=iUserBC
        UserBCFound(iUserBC)=.TRUE.
      END IF
    END DO
  END DO
END IF
DO iUserBC=1,nUserBCs
  IF(.NOT.UserBCFound(iUserBC)) CALL Abort(__STAMP__,&
    'Boundary condition specified in parameter file has not been found: '//TRIM(BoundaryName(iUserBC)))
END DO
DEALLOCATE(UserBCFound)

! Read boundary types from data file
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF((HSize(1).NE.4).OR.(HSize(2).NE.nBCs)) STOP 'Problem in readBC'
DEALLOCATE(HSize)
ALLOCATE(BCType(4,nBCs))
offset=0
CALL ReadArray('BCType',2,(/4,nBCs/),Offset,1,IntegerArray=BCType)
! Now apply boundary mappings
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    IF(BCMapping(iBC) .NE. 0)THEN
      IF((BoundaryType(BCMapping(iBC),1).EQ.1).AND.(BCType(1,iBC).NE.1)) &
        CALL abort(__STAMP__,'Remapping non-periodic to periodic BCs is not possible!')
      SWRITE(Unit_StdOut,'(A,A)')    ' |     Boundary in HDF file found |  ',TRIM(BCNames(iBC))
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                            was | ',BCType(1,iBC),BCType(3,iBC)
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                      is set to | ',BoundaryType(BCMapping(iBC),1:2)
      BCType(1,iBC) = BoundaryType(BCMapping(iBC),BC_TYPE)
      BCType(3,iBC) = BoundaryType(BCMapping(iBC),BC_STATE)
    END IF
  END DO
END IF
IF(ALLOCATED(BoundaryName)) DEALLOCATE(BoundaryName)
IF(ALLOCATED(BoundaryType)) DEALLOCATE(BoundaryType)
ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,3))
BoundaryName = BCNames
BoundaryType(:,BC_TYPE)  = BCType(1,:)
BoundaryType(:,BC_STATE) = BCType(3,:)
BoundaryType(:,BC_ALPHA) = BCType(4,:)
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(Unit_StdOut,'(A,A16,A20,A10,A10,A10)')'BOUNDARY CONDITIONS','|','Name','Type','State','Alpha'
DO iBC=1,nBCs
  SWRITE(*,'(A,A33,A20,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
END DO
SWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs


SUBROUTINE ReadMesh(FileString,NodeCoords)
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
USE MOD_Mesh_Vars,          ONLY:Elems
USE MOD_Mesh_Vars,          ONLY:aElem,aSide,bSide
USE MOD_Mesh_Vars,          ONLY:GETNEWELEM,GETNEWSIDE
#ifdef MPI
USE MOD_LoadBalance_Vars,   ONLY:nLoadBalance, LoadDistri, PartDistri,ParticleMPIWeight,WeightSum,TargetWeight
USE MOD_LoadDistribution,   ONLY:SingleStepOptimalPartition
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
REAL,ALLOCATABLE,INTENT(OUT)   :: NodeCoords(:,:,:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: NodeCoordsTmp(:,:,:,:,:)
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:)
INTEGER                        :: BCindex
INTEGER                        :: iElem,ElemID
INTEGER                        :: iLocSide,nbLocSide
INTEGER                        :: iSide
INTEGER                        :: FirstSideInd,LastSideInd,FirstElemInd,LastElemInd
INTEGER                        :: nPeriodicSides 
INTEGER                        :: ReduceData(7)
INTEGER                        :: nSideIDs,offsetSideID
#ifdef MPI
INTEGER                        :: ReduceData_glob(7)
INTEGER                        :: iNbProc, locnPart
INTEGER                        :: iProc, curiElem, MyElems, jProc,NewElems
INTEGER,ALLOCATABLE            :: MPISideCount(:)
INTEGER,ALLOCATABLE            :: PartInt(:,:)
INTEGER,PARAMETER              :: ELEM_FirstPartInd=1
INTEGER,PARAMETER              :: ELEM_LastPartInd=2
REAL,ALLOCATABLE               :: ElemWeight(:)
REAL                           :: CurWeight
! new weight distribution method
INTEGER                        :: BalanceMethod,iDistriIter
LOGICAL                        :: FoundDistribution
REAL                           :: LastProcDiff
REAL                           :: LoadDiff(0:nProcessors-1)
REAL                           :: MaxLoadDiff,LastLoadDiff
INTEGER                        :: ElemDistri(0:nProcessors-1),getElem
INTEGER,ALLOCATABLE            :: PartsInElem(:)
REAL                           :: diffLower,diffUpper
#endif
LOGICAL                        :: fileExists
LOGICAL                        :: doConnection
LOGICAL                        :: oriented
!===================================================================================================================================
IF(MESHInitIsDone) RETURN
IF(MPIRoot)THEN
  INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
  IF(.NOT.FileExists)  CALL abort(__STAMP__, &
          'readMesh from data file "'//TRIM(FileString)//'" does not exist')
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
IF(HSize(1).NE.6) &
  CALL abort(__STAMP__,'ERROR: Wrong size of ElemInfo, should be 6')
nGlobalElems=INT(HSize(2),4) !global number of elements
DEALLOCATE(HSize)
#ifdef MPI
!simple partition: nGlobalelems/nprocs, do this on proc 0
IF(ALLOCATED(offsetElemMPI)) DEALLOCATE(offsetElemMPI)
ALLOCATE(offsetElemMPI(0:nProcessors))
offsetElemMPI=0
SDEALLOCATE(LoadDistri)
ALLOCATE(LoadDistri(0:nProcessors-1))
LoadDistri(:)=0.
SDEALLOCATE(PartDistri)
ALLOCATE(PartDistri(0:nProcessors-1))
PartDistri(:)=0
IF (DoRestart) THEN 
  CALL CloseDataFile() 
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.)
  ALLOCATE(PartInt(1:nGlobalElems,2))
  PartInt(:,:)=0
  CALL ReadArray('PartInt',2,(/nGlobalElems,2/),0,1,IntegerArray=PartInt)
  CALL CloseDataFile() 
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)
  
  WeightSum = 0.0
  CurWeight = 0.0

  ! load balancing for particles
  ! read in particle data
  BalanceMethod     = GETINT('WeightDistributionMethod','0')

  ALLOCATE(ElemWeight(1:nGlobalElems))
  ALLOCATE(PartsInElem(1:nGlobalElems))

  DO iElem = 1, nGlobalElems
    locnPart=PartInt(iElem,ELEM_LastPartInd)-PartInt(iElem,ELEM_FirstPartInd)
    PartsInElem(iElem)=locnPart
    ElemWeight(iElem) = locnPart*ParticleMPIWeight + 1.0
    WeightSum = WeightSum + ElemWeight(iElem)
  END DO

  SELECT CASE(BalanceMethod)
  CASE(0) ! old scheme
    curiElem = 1
    WeightSum=WeightSum/REAL(nProcessors)
    DO iProc=0, nProcessors-1
      offsetElemMPI(iProc)=curiElem - 1 
      DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc  
        CurWeight=CurWeight+ElemWeight(iElem)  
        IF (CurWeight.GE.WeightSum*(iProc+1)) THEN
          curiElem = iElem + 1 
          EXIT
        END IF
      END DO   
    END DO
    offsetElemMPI(nProcessors)=nGlobalElems
  CASE(1)
    ! 1: last Proc receives the least load
    ! 2: Root receives the least load
    IF(MPIRoot)THEN
      FoundDistribution=.FALSE.
      targetWeight=WeightSum/REAL(nProcessors)
      LastProcDiff=0.
      iDistriIter=0
      DO WHILE(.NOT.FoundDistribution)
        iDistriIter=iDistriIter+1
        SWRITE(*,'(A19,I4,A19,G0)') '... LoadDistriIter ',iDistriIter,' with targetWeight=',targetWeight
        targetWeight=targetWeight+LastProcDiff/REAL(nProcessors)
        curiElem=1
        offSetElemMPI=0
        offsetElemMPI(nProcessors)=nGlobalElems
        LoadDistri=0.
        LoadDiff=0.
        DO iProc=0,nProcessors-1
          offSetElemMPI(iProc)=curiElem-1
          CurWeight=0.
          getElem=0
          DO iElem=curiElem, nGlobalElems - nProcessors +1 + iProc
            CurWeight=CurWeight+ElemWeight(iElem)
            getElem=getElem+1
            IF((CurWeight.GT.targetWeight) .OR. (iElem .EQ. nGlobalElems - nProcessors +1 + iProc))THEN
              diffLower=CurWeight-ElemWeight(iElem)-targetWeight
              diffUpper=Curweight-targetWeight
              IF(getElem.GT.1)THEN
                IF(iProc.EQ.nProcessors-1)THEN
                  LoadDistri(iProc)=CurWeight
                  LoadDiff(iProc)=diffUpper
                  curiElem=iElem+1
                  EXIT
                ELSE
                  IF(ABS(diffLower).LT.ABS(diffUpper) .AND. iElem.LT.nGlobalElems-nProcessors+1+iProc)THEN
                   LoadDiff(iProc)=diffLower
                   curiElem=iElem
                   LoadDistri(iProc)=CurWeight-ElemWeight(iElem)
                   EXIT
                  ELSE
                    LoadDiff(iProc)=diffUpper
                    curiElem=iElem+1
                    LoadDistri(iProc)=CurWeight
                    EXIT
                  END IF
                END IF
              ELSE
                LoadDiff(iProc)=diffUpper
                curiElem=iElem+1
                LoadDistri(iProc)=CurWeight
                EXIT
              END IF
            END IF
          END DO ! iElem
        END DO ! iProc
        ElemDistri=0
        DO iProc=0,nProcessors-1
          ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
          ! sanity check
          IF(ElemDistri(iProc).LE.0) CALL abort(&
          __STAMP__,' Process received zero elements during load distribution',iProc)
        END DO ! iPRoc
        LoadDistri(nProcessors-1)=ElemDistri(nProcessors-1) +&
                                  SUM(PartsInElem(offSetElemMPI(nProcessors-1)+1:nGlobalElems))*ParticleMPIWeight
        LastLoadDiff = LoadDistri(nProcessors-1)-targetWeight
        LoadDiff(nProcessors-1)=LastLoadDiff
        MaxLoadDiff=MAXVAL(LoadDiff(0:nProcessors-2))
        LastProcDiff=LastLoadDiff-MaxLoadDiff
        IF(LastProcDiff.LT.0.01*targetWeight)THEN
          FoundDistribution=.TRUE.
        END IF
        IF(iDistriIter.GT.nProcessors) THEN
         ! CALL abort(&
         ! __STAMP__,&
          SWRITE(UNIT_StdOut,'(A)') &
          'No valid load distribution throughout the processes found! Alter ParticleMPIWeight!'
          FoundDistribution=.TRUE.
        END IF         
        IF(ABS(WeightSum-SUM(LoadDistri)).GT.0.5) THEN
           CALL abort(&
          __STAMP__,' Lost Elements and/or Particles during load distribution!')
        END IF
      END DO
    END IF
    ! thend the ound distribution to all other procs
    CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  CASE(2)
     CALL abort(&
     __STAMP__,' error in load distritubion. please fix me!')
    ! 1: last Proc receives the least load
    ! 2: Root receives the least load
    IF(MPIRoot)THEN
      FoundDistribution=.FALSE.
      targetWeight=WeightSum/REAL(nProcessors)
      LastProcDiff=0.
      iDistriIter=0
      DO WHILE(.NOT.FoundDistribution)
        iDistriIter=iDistriIter+1
        SWRITE(*,*) 'LoadDistriIter', iDistriIter
        targetWeight=targetWeight+LastProcDiff/REAL(nProcessors)
        curiElem=nGlobalElems
        offSetElemMPI=0
        LoadDistri=0.
        LoadDiff=0.
        DO iProc=nProcessors-1,0,-1
          offSetElemMPI(iProc+1)=curiElem
          CurWeight=0.
          getElem=0
          DO iElem=curiElem,  iProc+1,-1
            CurWeight=CurWeight+ElemWeight(iElem)
            getElem=getElem+1
            IF((CurWeight.GT.targetWeight) .OR. (iElem .EQ. iProc+1))THEN ! take lower and upper is special case
              diffLower=CurWeight-ElemWeight(iElem)-targetWeight
              diffUpper=Curweight-targetWeight
              IF(getElem.GT.1)THEN
                IF(iProc.EQ.0)THEN
                  LoadDistri(iProc)=CurWeight
                  LoadDiff(iProc)=diffLower
                  curiElem=iElem
                  EXIT
                ELSE
                  IF(ABS(diffLower).GT.ABS(diffUpper) .AND. iElem.GT.iProc+1)THEN
                    ! take upper
                    LoadDiff(iProc)=diffUpper
                    curiElem=iElem+1
                    LoadDistri(iProc)=CurWeight-ElemWeight(iElem)
                    EXIT
                  ELSE
                    LoadDiff(iProc)=diffLower
                    curiElem=iElem
                    LoadDistri(iProc)=CurWeight!+ElemWeight(iElem)
                    EXIT
                  END IF
                END IF
              ELSE
                LoadDiff(iProc)=diffLower
                curiElem=iElem
                LoadDistri(iProc)=CurWeight
                EXIT
              END IF
            END IF
          END DO ! iElem
        END DO ! iProc
        ElemDistri=0
        DO iProc=0,nProcessors-1
          ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
          ! sanity check
          IF(ElemDistri(iProc).LE.0) CALL abort(&
          __STAMP__,' Process received zero elements during load distribution',iProc)
        END DO ! iPRoc
        LoadDistri(0)=ElemDistri(0) +&
                                  SUM(PartsInElem(1:offSetElemMPI(1)))*ParticleMPIWeight
        LastLoadDiff = LoadDistri(0)-targetWeight
        LoadDiff(0)=LastLoadDiff
        MaxLoadDiff=MAXVAL(LoadDiff(1:nProcessors-1))
        LastProcDiff=LastLoadDiff-MaxLoadDiff
        IF(LastProcDiff.GT.0)THEN
          FoundDistribution=.TRUE.
        END IF
        IF(iDistriIter.EQ.nProcessors) CALL abort(&
          __STAMP__,&
          'No valid load distribution throughout the processes found! Alter ParticleMPIWeight!')
        IF(ABS(WeightSum-SUM(LoadDistri)).GT.0.5) THEN
           WRITE(*,*) WeightSum-SUM(LoadDistri)
           WRITE(*,*) OffSetElemMPI(1)
           WRITE(*,*) ElemDistri
           WRITE(*,*) LoadDistri
           CALL abort(&
          __STAMP__,' Lost Elements and/or Particles during load distribution!')
        END IF
      END DO
    END IF
    ! thend the ound distribution to all other procs
    CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  CASE(3)
    ! 1: last Proc receives the least load
    targetWeight=WeightSum/REAL(nProcessors)
    LastProcDiff=0.
    curiElem=1
    offSetElemMPI=0
    offsetElemMPI(nProcessors)=nGlobalElems
    LoadDistri=0.
    LoadDiff=0.
    DO iProc=0,nProcessors-1
      offSetElemMPI(iProc)=curiElem-1
      CurWeight=0.
      getElem=0
      DO iElem=curiElem, nGlobalElems - nProcessors +1 + iProc
        CurWeight=CurWeight+ElemWeight(iElem)
        getElem=getElem+1
        IF((CurWeight.GT.targetWeight) .OR. (iElem .EQ. nGlobalElems - nProcessors +1 + iProc))THEN
          diffLower=CurWeight-ElemWeight(iElem)-targetWeight
          diffUpper=Curweight-targetWeight
          IF(getElem.GT.1)THEN
            IF(iProc.EQ.nProcessors-1)THEN
              LoadDistri(iProc)=CurWeight
              LoadDiff(iProc)=diffUpper
              curiElem=iElem+1
              EXIT
            ELSE
              IF(ABS(diffLower).LT.ABS(diffUpper) .AND. iElem.LT.nGlobalElems-nProcessors+1+iProc)THEN
               LoadDiff(iProc)=diffLower
               curiElem=iElem
               LoadDistri(iProc)=CurWeight-ElemWeight(iElem)
               EXIT
              ELSE
                LoadDiff(iProc)=diffUpper
                curiElem=iElem+1
                LoadDistri(iProc)=CurWeight
                EXIT
              END IF
            END IF
          ELSE
            LoadDiff(iProc)=diffUpper
            curiElem=iElem+1
            LoadDistri(iProc)=CurWeight
            EXIT
          END IF
        END IF
      END DO ! iElem
    END DO ! iProc
    ElemDistri=0
    DO iProc=0,nProcessors-1
      ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
      ! sanity check
      IF(ElemDistri(iProc).LE.0) CALL abort(&
      __STAMP__,' Process received zero elements during load distribution',iProc)
    END DO ! iPRoc
    ! redistribute element weight
    DO iProc=1,nProcessors
      FirstElemInd=OffSetElemMPI(MyRank)+1
      LastElemInd =OffSetElemMPI(MyRank+1)
      MyElems=ElemDistri(MyRank)
      CALL SingleStepOptimalPartition(MyElems,NewElems,ElemWeight(FirstElemInd:LastElemInd))
      ElemDistri=0
      CALL MPI_ALLGATHER(NewElems,1,MPI_INTEGER,ElemDistri(:),1,MPI_INTEGER,MPI_COMM_WORLD,iERROR)
      ! calculate proc offset
      OffSetElemMPI(0)=0
      DO jProc=0,nProcessors-1
        OffSetElemMPI(jProc+1) = OffsetElemMPI(jProc) + ElemDistri(jProc)
      END DO ! jProc=0,nProcessors-1
    END DO ! iProc=1,nProcessors
    ! compute load distri
    DO iProc=0,nProcessors-1
      FirstElemInd=OffSetElemMPI(iProc)+1
      LastElemInd =OffSetElemMPI(iProc+1)
      LoadDistri(iProc) = LastElemInd-OffSetElemMPI(iProc) &
                        + SUM(PartsInElem(FirstElemInd:LastElemInd))*ParticleMPIWeight
    !  SWRITE(*,*) FirstElemInd,LastElemInd,LoadDistri(iProc),SUM(PartsInElem(FirstElemInd:LastElemInd))
    END DO ! iPRoc
  CASE(4)
    ! predistribute elements
    curiElem = 1
    targetWeight=WeightSum/REAL(nProcessors)
    offsetElemMPI(nProcessors)=nGlobalElems
    DO iProc=0, nProcessors-1
      offsetElemMPI(iProc)=curiElem - 1 
      DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc  
        CurWeight=CurWeight+ElemWeight(iElem)  
        IF (CurWeight.GE.targetWeight*(iProc+1)) THEN
          curiElem = iElem + 1 
          EXIT
        END IF
      END DO   
    END DO
    ElemDistri=0
    DO iProc=0,nProcessors-1
      ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
      ! sanity check
      IF(ElemDistri(iProc).LE.0) CALL abort(&
      __STAMP__,' Process received zero elements during load distribution',iProc)
    END DO ! iPRoc
    ! redistribute element weight
    DO iProc=1,nProcessors
      FirstElemInd=OffSetElemMPI(MyRank)+1
      LastElemInd =OffSetElemMPI(MyRank+1)
      MyElems=ElemDistri(MyRank)
      CALL SingleStepOptimalPartition(MyElems,NewElems,ElemWeight(FirstElemInd:LastElemInd))
      ElemDistri=0
      CALL MPI_ALLGATHER(NewElems,1,MPI_INTEGER,ElemDistri(:),1,MPI_INTEGER,MPI_COMM_WORLD,iERROR)
      ! calculate proc offset
      OffSetElemMPI(0)=0
      DO jProc=0,nProcessors-1
        OffSetElemMPI(jProc+1) = OffsetElemMPI(jProc) + ElemDistri(jProc)
      END DO ! jProc=1,nProcessors
    END DO ! jProc=0,nProcessors
  ! compute load distri
    LoadDistri=0.
    DO iProc=0,nProcessors-1
      FirstElemInd=OffSetElemMPI(iProc)+1
      LastElemInd =OffSetElemMPI(iProc+1)
      LoadDistri(iProc) = LastElemInd-OffSetElemMPI(iProc) &
                        + SUM(PartsInElem(FirstElemInd:LastElemInd))*ParticleMPIWeight
    !  SWRITE(*,*) FirstElemInd,LastElemInd,LoadDistri(iProc),SUM(PartsInElem(FirstElemInd:LastElemInd))
    END DO ! iPRoc
  END SELECT
  offsetElemMPI(nProcessors)=nGlobalElems
  DO iProc=0,nProcessors-1
    DO iElem=offSetElemMPI(iProc)+1,offSetElemMPI(iProc+1)
      PartDistri(iProc)=PartDistri(iProc)+PartsInElem(iElem)
    END DO
  END DO ! iProc
  DEALLOCATE(PartsInElem)
ELSE
  nElems=nGlobalElems/nProcessors
  iElem=nGlobalElems-nElems*nProcessors
  DO iProc=0,nProcessors-1
    offsetElemMPI(iProc)=nElems*iProc+MIN(iProc,iElem)
  END DO
  offsetElemMPI(nProcessors)=nGlobalElems
END IF
!local nElems and offset
nElems=offsetElemMPI(myRank+1)-offsetElemMPI(myRank)
offsetElem=offsetElemMPI(myRank)
LOGWRITE(*,*)'offset,nElems',offsetElem,nElems
#else /* MPI */
nElems=nGlobalElems   !local number of Elements 
offsetElem=0          ! offset is the index of first entry, hdf5 array starts at 0-.GT. -1 
#endif /* MPI */

CALL readBCs()
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
FirstElemInd=offsetElem+1
LastElemInd=offsetElem+nElems
ALLOCATE(Elems(                FirstElemInd:LastElemInd))
ALLOCATE(ElemInfo(ElemInfoSize,FirstElemInd:LastElemInd))
CALL ReadArray('ElemInfo',2,(/ElemInfoSize,nElems/),offsetElem,2,IntegerArray=ElemInfo)

DO iElem=FirstElemInd,LastElemInd
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  Elems(iElem)%ep=>GETNEWELEM()
  aElem=>Elems(iElem)%ep
  aElem%Ind    = iElem
  aElem%Type   = ElemInfo(ELEM_Type,iElem)
  aElem%Zone   = ElemInfo(ELEM_Zone,iElem)
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------

offsetSideID=ElemInfo(ELEM_FirstSideInd,FirstElemInd) ! hdf5 array starts at 0-> -1  
nSideIDs=ElemInfo(ELEM_LastSideInd,LastElemInd)-ElemInfo(ELEM_FirstSideInd,FirstElemInd)
!read local SideInfo from data file 
FirstSideInd=offsetSideID+1
LastSideInd=offsetSideID+nSideIDs
ALLOCATE(SideInfo(SideInfoSize,FirstSideInd:LastSideInd))
CALL ReadArray('SideInfo',2,(/SideInfoSize,nSideIDs/),offsetSideID,2,IntegerArray=SideInfo)


DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  !build up sides of the element according to CGNS standard
  ! assign flip
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    aSide%Elem=>aElem
    oriented=(Sideinfo(SIDE_ID,iSide).GT.0)
    aSide%Ind=ABS(SideInfo(SIDE_ID,iSide))
    IF(oriented)THEN !oriented side
      aSide%flip=0
    ELSE !not oriented
      aSide%flip=MOD(Sideinfo(SIDE_Flip,iSide),10)
      IF((aSide%flip.LT.0).OR.(aSide%flip.GT.4)) STOP 'NodeID doesnt belong to side'
    END IF
  END DO !i=1,locnSides
END DO !iElem

 
! build up side connection 
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    elemID  = SideInfo(SIDE_nbElemID,iSide)
    BCindex = SideInfo(SIDE_BCID,iSide)
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

    IF(.NOT.doConnection) CYCLE
    IF(ASSOCIATED(aSide%connection)) CYCLE

    ! check if neighbor on local proc or MPI connection
    IF((elemID.LE.LastElemInd).AND.(elemID.GE.FirstElemInd))THEN !local
      nbLocSide=Sideinfo(SIDE_Flip,iSide)/10
      IF((nbLocSide.LT.1).OR.(nbLocSide.GT.6))&
        CALL abort(__STAMP__,'SideInfo: Index of local side must be between 1 and 6!')
      bSide=>Elems(elemID)%ep%Side(nbLocSide)%sp
      aSide%connection=>bSide
      bSide%connection=>aSide
      IF(bSide%ind.NE.aSide%ind)&
        CALL abort(__STAMP__,&
        'SideInfo: Index of side and neighbor side have to be identical!')
    ELSE !MPI
#ifdef MPI
      aSide%connection=>GETNEWSIDE()            
      aSide%connection%flip=aSide%flip
      aSide%connection%Elem=>GETNEWELEM()
      aSide%NbProc = ELEMIPROC(elemID)
#else
      CALL abort(__STAMP__, &
        ' elemID of neighbor not in global Elem list ')
#endif
    END IF
  END DO !iLocSide 
END DO !iElem

DEALLOCATE(ElemInfo,SideInfo)

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------
! get physical coordinates
IF(useCurveds)THEN
  ALLOCATE(NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),offsetElem*(NGeo+1)**3,2,RealArray=NodeCoords)
ELSE
  ALLOCATE(NodeCoords(   3,0:1,   0:1,   0:1,   nElems))
  ALLOCATE(NodeCoordsTmp(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),offsetElem*(NGeo+1)**3,2,RealArray=NodeCoordsTmp)
  NodeCoords(:,0,0,0,:)=NodeCoordsTmp(:,0,   0,   0,   :)
  NodeCoords(:,1,0,0,:)=NodeCoordsTmp(:,NGeo,0,   0,   :)
  NodeCoords(:,0,1,0,:)=NodeCoordsTmp(:,0,   NGeo,0,   :)
  NodeCoords(:,1,1,0,:)=NodeCoordsTmp(:,NGeo,NGeo,0,   :)
  NodeCoords(:,0,0,1,:)=NodeCoordsTmp(:,0,   0,   NGeo,:)
  NodeCoords(:,1,0,1,:)=NodeCoordsTmp(:,NGeo,0,   NGeo,:)
  NodeCoords(:,0,1,1,:)=NodeCoordsTmp(:,0,   NGeo,NGeo,:)
  NodeCoords(:,1,1,1,:)=NodeCoordsTmp(:,NGeo,NGeo,NGeo,:)
  DEALLOCATE(NodeCoordsTmp)
  NGeo=1
ENDIF
nNodes=nElems*(NGeo+1)**3

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
          !nPeriodicSides=nPeriodicSides+1
          IF(BoundaryType(aSide%BCindex,BC_TYPE).EQ.1) nPeriodicSides=nPeriodicSides+1
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


END MODULE MOD_Mesh_ReadIn
