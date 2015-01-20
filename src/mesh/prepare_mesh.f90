#include "boltzplatz.h"

MODULE MOD_Prepare_Mesh
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

INTERFACE setLocalSideIDs
  MODULE PROCEDURE setLocalSideIDs
END INTERFACE

INTERFACE fillMeshInfo
  MODULE PROCEDURE fillMeshInfo
END INTERFACE

INTERFACE fillElemGeo
  MODULE PROCEDURE fillElemGeo
END INTERFACE

INTERFACE getVolumeMapping
  MODULE PROCEDURE getVolumeMapping
END INTERFACE

PUBLIC::setLocalSideIDs,fillMeshInfo,fillElemGeo,getVolumeMapping

#ifdef MPI
INTERFACE exchangeFlip
  MODULE PROCEDURE exchangeFlip
END INTERFACE

PUBLIC::exchangeFlip 
#endif
!===================================================================================================================================

CONTAINS


SUBROUTINE setLocalSideIDs()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY: nElems,nInnerSides,nBCSides,offsetElem
USE MOD_Mesh_Vars,  ONLY: aElem,aSide
USE MOD_Mesh_Vars,  ONLY: Elems,nMPISides_MINE,nMPISides_YOUR
#ifdef MPI
USE MOD_Mesh_Vars,  ONLY: nSides
USE MOD_MPI_Vars,   ONLY: nNbProcs,NbProc,nMPISides_Proc,nMPISides_MINE_Proc,nMPISides_YOUR_Proc
USE MOD_MPI_Vars,   ONLY: offsetElemMPI,offsetMPISides_MINE,offsetMPISides_YOUR
USE MOD_Mesh_ReadIn,ONLY: Qsort1Int,INVMAP
#endif
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: iElem,FirstElemInd,LastElemInd
INTEGER   :: iLocSide,iSide,iInnerSide,iBCSide
#ifdef MPI
INTEGER               :: iNbProc,ioUnit
INTEGER,ALLOCATABLE   :: SideIDMap(:)
CHARACTER(LEN=10)     :: formatstr
INTEGER,ALLOCATABLE   :: NBinfo(:,:),NBinfo_glob(:,:,:),nNBProcs_glob(:),Procinfo_glob(:,:),tmparray(:,:)  !for output only
REAL,ALLOCATABLE      :: tmpreal(:,:)
INTEGER               :: i,j,ProcInfo(4),nNBmax      !for output only
#endif
!===================================================================================================================================
FirstElemInd= offsetElem+1
LastElemInd = offsetElem+nElems
! ----------------------------------------
! Set side IDs to arrange sides:
! 1. BC sides
! 2. inner sides
! 3. MPI sides
! MPI Sides are not included here!
! ----------------------------------------

DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
      aElem%Side(iLocSide)%sp%sideID=-1
  END DO
END DO

iSide=0
iBCSide=0
iInnerSide=nBCSides
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%sideID.EQ.-1)THEN
      IF(aSide%NbProc.EQ.-1)THEN ! no MPI Sides
        IF(ASSOCIATED(aSide%connection))THEN
          iInnerSide=iInnerSide+1
          iSide=iSide+1
          aSide%SideID=iInnerSide
          aSide%connection%SideID=iInnerSide
        ELSE
          iBCSide=iBCSide+1
          iSide=iSide+1
          aSide%SideID=iBCSide
        END IF !associated connection
      END IF ! .NOT. MPISide
    END IF !sideID NE -1
  END DO ! iLocSide=1,6
END DO !iElem
IF(iSide.NE.(nInnerSides+nBCSides)) STOP 'not all SideIDs are set!'

nMPISides_MINE=0
nMPISides_YOUR=0
#ifdef MPI
! SPLITTING MPISides in MINE and YOURS
ALLOCATE(nMPISides_MINE_Proc(1:nNbProcs),nMPISides_YOUR_Proc(1:nNbProcs))
nMPISides_MINE_Proc=0
nMPISides_YOUR_Proc=0
DO iNbProc=1,nNbProcs
  IF(myRank.LT.NbProc(iNbProc)) THEN
    nMPISides_MINE_Proc(iNbProc)=nMPISides_Proc(iNbProc)/2
  ELSE
    nMPISides_MINE_Proc(iNbProc)=nMPISides_Proc(iNbProc)-nMPISides_Proc(iNbProc)/2
  END IF    
  nMPISides_YOUR_Proc(iNbProc)=nMPISides_Proc(iNbProc)-nMPISides_MINE_Proc(iNbProc)
END DO
nMPISides_MINE=SUM(nMPISides_MINE_Proc)
nMPISides_YOUR=SUM(nMPISides_YOUR_Proc)

ALLOCATE(offsetMPISides_YOUR(0:nNbProcs),offsetMPISides_MINE(0:nNbProcs))
offsetMPISides_MINE=0
offsetMPISides_YOUR=0
! compute offset, first all MINE , then all YOUR MPISides
offsetMPISides_MINE(0)=nInnerSides+nBCSides
DO iNbProc=1,nNbProcs
  offsetMPISides_MINE(iNbProc)=offsetMPISides_MINE(iNbProc-1)+nMPISides_MINE_Proc(iNbProc)
END DO
offsetMPISides_YOUR(0)=offsetMPISides_MINE(nNbProcs)
DO iNbProc=1,nNbProcs
  offsetMPISides_YOUR(iNbProc)=offsetMPISides_YOUR(iNbProc-1)+nMPISides_YOUR_Proc(iNbProc)
END DO
IF(nProcessors.EQ.1) RETURN
DO iNbProc=1,nNbProcs
  ALLOCATE(SideIDMap(nMPISides_Proc(iNbProc)))
  iSide=0
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      IF(aSide%NbProc.NE.NbProc(iNbProc))CYCLE
      iSide=iSide+1
      SideIDMap(iSide)=aSide%ind !global Side Index 
    END DO !iLocSide
  END DO !iElem
  CALL Qsort1Int(SideIDMap) !sort by global side index
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      IF(aSide%NbProc.NE.NbProc(iNbProc))CYCLE
      aSide%SideID=INVMAP(aSide%ind,nMPISides_Proc(iNbProc),SideIDMap) ! get sorted iSide
      IF(myRank.LT.aSide%NbProc)THEN
        IF(aSide%SideID.LE.nMPISides_MINE_Proc(iNbProc))THEN !MINE
          aSide%SideID=aSide%SideID +offsetMPISides_MINE(iNbProc-1)
        ELSE !YOUR
          aSide%SideID=(aSide%SideID-nMPISides_MINE_Proc(iNbProc))+offsetMPISides_YOUR(iNbProc-1)
        END IF
      ELSE
        IF(aSide%SideID.LE.nMPISides_YOUR_Proc(iNbProc))THEN !MINE
          aSide%SideID=aSide%SideID +offsetMPISides_YOUR(iNbProc-1)
        ELSE !YOUR
          aSide%SideID=(aSide%SideID-nMPISides_YOUR_Proc(iNbProc))+offsetMPISides_MINE(iNbProc-1)
        END IF
      END IF !myrank<NbProc
    END DO !iLocSide
  END DO !iElem
  DEALLOCATE(SideIDMap)
END DO !nbProc(i)


WRITE(formatstr,'(a5,I2,a3)')'(A22,',nNBProcs,'I8)'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nNbProcs:',nNbProcs
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'NbProc:'   ,NbProc
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'nMPISides_Proc:',nMPISides_Proc
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'nMPISides_MINE_Proc:',nMPISides_MINE_Proc
LOGWRITE(*,formatstr)'nMPISides_YOUR_Proc:',nMPISides_YOUR_Proc
WRITE(formatstr,'(a5,I2,a3)')'(A22,',nNBProcs+1,'I8)'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'offsetMPISides_MINE:',offsetMPISides_MINE
LOGWRITE(*,formatstr)'offsetMPISides_YOUR:',offsetMPISides_YOUR
LOGWRITE(*,*)'-------------------------------------------------------'

!output partitioning info
ProcInfo(1)=nElems
ProcInfo(2)=nSides
ProcInfo(3)=nInnerSides
ProcInfo(4)=nBCSides
IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:nProcessors-1))
  ALLOCATE(ProcInfo_glob(4,0:nProcessors-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot 
CALL MPI_GATHER(nNBProcs,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
CALL MPI_GATHER(ProcInfo,4,MPI_INTEGER,ProcInfo_glob,4,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(6,nNBmax,0:nProcessors))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError) 
ALLOCATE(NBinfo(6,nNbmax))
NBinfo=0
NBinfo(1,1:nNBProcs)=NBProc
NBinfo(2,1:nNBProcs)=nMPISides_Proc
NBinfo(3,1:nNBProcs)=nMPISides_MINE_Proc
NBinfo(4,1:nNBProcs)=nMPISides_YOUR_Proc
NBinfo(5,1:nNBProcs)=offsetMPISides_MINE(0:nNBProcs-1)
NBinfo(6,1:nNBProcs)=offsetMPISides_YOUR(0:nNBProcs-1)
CALL MPI_GATHER(NBinfo,6*nNBmax,MPI_INTEGER,NBinfo_glob,6*nNBmax,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  OPEN(UNIT=ioUnit,FILE='partitionInfo.out',STATUS='REPLACE')
  WRITE(ioUnit,*)'Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',nProcessors
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  WRITE(ioUnit,'(8(A15))')'Rank','nElems','nSides','nInnerSides','nBCSides','nMPISides','nMPISides_MINE','nNBProcs'
  WRITE(ioUnit,'(A120)')&
      '======================================================================================================================='
  !statistics
  ALLOCATE(tmparray(7,0:3),tmpreal(7,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=HUGE(-1)  !max
  tmparray(:,3)=HUGE(1)   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i)
    tmparray(2,0)=Procinfo_glob(2,i)
    tmparray(3,0)=Procinfo_glob(3,i)
    tmparray(4,0)=Procinfo_glob(4,i)
    tmparray(5,0)=SUM(NBinfo_glob(2,:,i))
    tmparray(6,0)=SUM(NBinfo_glob(3,:,i))
    tmparray(7,0)=nNBProcs_glob(i)
    DO j=1,7
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO
  END DO
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(nProcessors) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i)
    tmparray(2,0)=Procinfo_glob(2,i)
    tmparray(3,0)=Procinfo_glob(3,i)
    tmparray(4,0)=Procinfo_glob(4,i)
    tmparray(5,0)=SUM(NBinfo_glob(2,:,i))
    tmparray(6,0)=SUM(NBinfo_glob(3,:,i))
    tmparray(7,0)=nNBProcs_glob(i)
    DO j=1,7
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2 
    END DO
  END DO
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(nProcessors))
  WRITE(ioUnit,'(A15,7(5X,F10.2))')'   MEAN        ',tmpreal(:,1)
  WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A15,7(5X,F10.2))')'   RMS         ',tmpreal(:,2)
  WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A15,7(5X,I10))')'   MIN         ',tmparray(:,3)
  WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A15,7(5X,I10))')'   MAX         ',tmparray(:,2)
  WRITE(ioUnit,'(A120)')&
      '======================================================================================================================='
  DO i=0,nProcessors-1
    WRITE(ioUnit,'(8(5X,I10))')i,Procinfo_glob(:,i),SUM(NBinfo_glob(2,:,i)),SUM(NBinfo_glob(3,:,i)),nNBProcs_glob(i)
    WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  END DO
  WRITE(ioUnit,*)' '
  WRITE(ioUnit,*)'Information per neighbor processor'
  WRITE(ioUnit,*)' '
  WRITE(ioUnit,'(7(A15))')'Rank','NBProc','nMPISides_Proc','nMPISides_MINE','nMPISides_YOUR','offset_MINE','offset_YOUR'
  WRITE(ioUnit,'(A120)')&
      '======================================================================================================================='
  DO i=0,nProcessors-1
    WRITE(ioUnit,'(7(5X,I10))')i,NBinfo_glob(:,1,i)
    DO j=2,nNBProcs_glob(i)
      WRITE(ioUnit,'(A15,6(5X,I10))')' ',NBinfo_glob(:,j,i)
    END DO
    WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  END DO
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit) 
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)
#endif /*MPI*/  
END SUBROUTINE setLocalSideIDs

SUBROUTINE fillMeshInfo()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems,offsetElem,nInnerSides,nBCSides
USE MOD_Mesh_Vars,ONLY:nMPISides_MINE
USE MOD_Mesh_Vars,ONLY:ElemToSide,BC,SideToElem,SideToElem2
USE MOD_Mesh_Vars,ONLY:aElem,aSide
USE MOD_Mesh_Vars,ONLY:Elems
#ifdef MPI
USE MOD_Mesh_Vars,ONLY:nMPISides
USE MOD_MPI_vars
#endif
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iSide,LocSideID,nSides_flip(0:4)
#ifdef MPI
INTEGER             :: dummy(0:4)
#endif
!===================================================================================================================================
!  LOGWRITE(*,'(4A8)')'SideID', 'globID','NbProc','Flip'
!  DO iElem=1,nElems
!    aElem=>Elems(iElem+offsetElem)%ep
!    DO LocSideID=1,6
!      aSide=>aElem%Side(LocSideID)%sp
!      LOGWRITE(*,'(4I8)') aSide%SideID,aSide%ind,aSide%NbProc,aSide%flip
!    END DO !iLocSide
!  END DO !iElem
! ELement to Side mapping
nSides_flip=0
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    ElemToSide(E2S_SIDE_ID,LocSideID,iElem)=aSide%SideID
    ElemToSide(E2S_FLIP,LocSideID,iElem)   =aSide%Flip
    nSides_flip(aSide%flip)=nSides_flip(aSide%flip)+1
  END DO ! LocSideID
END DO ! iElem

! Side to Element mapping, sorted by SideID

DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%Flip.EQ.0)THEN !root side
      SideToElem(S2E_ELEM_ID,aSide%SideID)         = iElem !root Element
      SideToElem(S2E_LOC_SIDE_ID,aSide%SideID)     = LocSideID
    ELSE
      SideToElem(S2E_NB_ELEM_ID,aSide%SideID)      = iElem ! element with flipped side
      SideToElem(S2E_NB_LOC_SIDE_ID,aSide%SideID)  = LocSideID
      SideToElem(S2E_FLIP,aSide%SideID)            = aSide%Flip
    END IF
    IF(aSide%sideID .LE. nBCSides)THEN
      BC(aSide%sideID)=aSide%BCIndex
    END IF
  END DO ! LocSideID
END DO ! iElem
! Side to Element mapping, sorted by iElem, only MINE are added

iSide=0
! first BCSides, InnerSides and MINE MPISides, attention, all mixed, since element-wise ordering
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%SideID.GT.nBCSides+nInnerSides+nMPISides_MINE)CYCLE
    iSide=iSide+1
    SideToElem2(S2E2_ELEM_ID,iSide)     = iElem
    SideToElem2(S2E2_SIDE_ID,iSide)     = aSide%SideID
    SideToElem2(S2E2_LOC_SIDE_ID,iSide) = LocSideID
    SideToElem2(S2E2_FLIP,iSide)        = aSide%flip
  END DO ! LocSideID
END DO ! iElem
IF(iSide.NE.(2*nInnerSides+nBCSides+nMPISides_MINE)) STOP 'wrong number of nInnerSides+nBCSides+nMPISides_MINE'
#ifdef MPI
! now the last Sides, only YOUR MPISides
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%SideID.LE.nBCSides+nInnerSides+nMPISides_MINE)CYCLE
    iSide=iSide+1
    SideToElem2(S2E2_ELEM_ID,iSide)     = iElem
    SideToElem2(S2E2_SIDE_ID,iSide)     = aSide%SideID
    SideToElem2(S2E2_LOC_SIDE_ID,iSide) = LocSideID
    SideToElem2(S2E2_FLIP,iSide)        = aSide%flip
  END DO ! LocSideID
END DO ! iElem
IF(iSide.NE.(2*nInnerSides+nBCSides+nMPISides)) STOP 'wrong number of nInnerSides+nBCSides+nMPISides_MINE'
#endif

#ifdef MPI
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nSides_flip,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(nSides_flip,dummy,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=0 | ',nSides_flip(0)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=1 | ',nSides_flip(1)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=2 | ',nSides_flip(2)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=3 | ',nSides_flip(3)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=4 | ',nSides_flip(4)
SWRITE(UNIT_StdOut,'(132("."))')
END SUBROUTINE fillMeshInfo


#ifdef MPI
SUBROUTINE exchangeFlip()
!===================================================================================================================================
! set flip of MINE sides to zero, therefore send flip of MINE to other processor, so that YOUR sides get their corresponding flip>0
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems,offsetElem
USE MOD_Mesh_Vars,ONLY:aElem,aSide
USE MOD_Mesh_Vars,ONLY:Elems,nBCSides,nInnerSides,nMPISides_MINE
USE MOD_MPI_vars
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,LocSideID
INTEGER             :: Flip_MINE(offsetMPISides_MINE(0)+1:offsetMPISides_MINE(nNBProcs))
INTEGER             :: Flip_YOUR(offsetMPISides_YOUR(0)+1:offsetMPISides_YOUR(nNBProcs))
INTEGER             :: SendRequest(nNbProcs),RecRequest(nNbProcs)
!===================================================================================================================================
IF(nProcessors.EQ.1) RETURN
!fill MINE flip info
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF((aSide%SideID.GT.nBCSides+nInnerSides).AND.&
       (aSide%SideID.LE.nBCSides+nInnerSides+nMPISides_MINE))THEN
      Flip_MINE(aSide%sideID)=aSide%flip
    END IF
  END DO ! LocSideID
END DO ! iElem
DO iNbProc=1,nNbProcs
  ! Start send flip from MINE 
  IF(nMPISides_MINE_Proc(iNbProc).GT.0)THEN
    nSendVal    =nMPISides_MINE_Proc(iNbProc)
    SideID_start=OffsetMPISides_MINE(iNbProc-1)+1
    SideID_end  =OffsetMPISides_MINE(iNbProc)
    CALL MPI_ISEND(Flip_MINE(SideID_start:SideID_end),nSendVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,SendRequest(iNbProc),iError)
  END IF
  ! Start receive flip to YOUR
  IF(nMPISides_YOUR_Proc(iNbProc).GT.0)THEN
    nRecVal     =nMPISides_YOUR_Proc(iNbProc)
    SideID_start=OffsetMPISides_YOUR(iNbProc-1)+1
    SideID_end  =OffsetMPISides_YOUR(iNbProc)
    CALL MPI_IRECV(Flip_YOUR(SideID_start:SideID_end),nRecVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,RecRequest(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs
DO iNbProc=1,nNbProcs
  IF(nMPISides_YOUR_Proc(iNbProc).GT.0)CALL MPI_WAIT(RecRequest(iNbProc) ,MPIStatus,iError)
  IF(nMPISides_MINE_Proc(iNBProc).GT.0)CALL MPI_WAIT(SendRequest(iNbProc),MPIStatus,iError)
END DO !iProc=1,nNBProcs
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%NbProc.EQ.-1) CYCLE !no MPISide
    IF(aSide%SideID.GT.nBCSides+nInnerSides+nMPISides_MINE)THEN
      IF(aSide%flip.EQ.0)THEN
        IF(Flip_YOUR(aSide%SideID).EQ.0) STOP 'problem in exchangeflip'
        aSide%flip=Flip_YOUR(aSide%sideID)
      END IF
    ELSE
      aSide%flip=0 !MINE MPISides flip=0
    END IF
  END DO ! LocSideID
END DO ! iElem
  
END SUBROUTINE exchangeFlip
#endif


SUBROUTINE fillElemGeo(XCL_NGeo)
!===================================================================================================================================
! Fill XCL_NGeo with interpolation points and transform those to Chebyshev-Lobatto points
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems,offsetElem
USE MOD_Mesh_Vars,ONLY:Elems
USE MOD_Mesh_Vars,ONLY:aElem
USE MOD_Mesh_Vars,ONLY:NGeo,Vdm_NGeo_CLNGeo
USE MOD_ChangeBasis,ONLY:changeBasis3D
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iNode,iGeo(3)
INTEGER            :: NodeMap((NGeo+1)**3,3)
!===================================================================================================================================
IF(NGeo.EQ.1)THEN
  NodeMap(1,:)=(/   0,   0,   0/)
  NodeMap(2,:)=(/NGeo,   0,   0/)
  NodeMap(3,:)=(/NGeo,NGeo,   0/)
  NodeMap(4,:)=(/   0,NGeo,   0/)
  NodeMap(5,:)=(/   0,   0,NGeo/)
  NodeMap(6,:)=(/NGeo,   0,NGeo/)
  NodeMap(7,:)=(/NGeo,NGeo,NGeo/)
  NodeMap(8,:)=(/   0,NGeo,NGeo/)
ELSE
  CALL getBasisMappingHexa(NGeo,NodeMap)
END IF

DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  IF(aElem%nCurvedNodes.EQ.0)THEN !Only use corner nodes
    DO iNode=1,8
      iGeo=NodeMap(iNode,:)
      XCL_NGeo(:,iGeo(1),iGeo(2),iGeo(3),iElem) = aElem%Node(iNode)%np%x
    END DO !iNode
  ELSE !NGeo > 1, so curved elements
    DO iNode=1,aElem%nCurvedNodes
      iGeo=NodeMap(iNode,:)
      XCL_NGeo(:,iGeo(1),iGeo(2),iGeo(3),iElem) = aElem%curvedNode(iNode)%np%x
    END DO !iNode
  END IF
  !from equidistant to CL
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_NGeo_CLNGeo,XCL_NGeo(:,:,:,:,iElem),XCL_NGeo(:,:,:,:,iElem))
END DO


END SUBROUTINE fillElemGeo

SUBROUTINE getBasisMappingHexa(N,bMap,bMapInv)
!===================================================================================================================================
! mapping from iAns -> i,j  in [0,Deg], can be used for nodeMap too: CALL getBasisMapping(nNodes1D-1,nodeMap) 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: N 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)          :: bMap((N+1)**3,3)
INTEGER,INTENT(OUT),OPTIONAL :: bMapInv(0:N,0:N,0:N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: iAns,i,j,k
!===================================================================================================================================
iAns=0
DO k=0,N; DO j=0,N; DO i=0,N
  iAns=iAns+1
  bMap(iAns,:)=(/i,j,k/)
END DO; END DO; END DO
IF(PRESENT(bMapInv))THEN
  bMapInv=0
  iAns=0
  DO k=0,N; DO j=0,N; DO i=0,N
    iAns=iAns+1
    bMapInv(i,j,k)=iAns
  END DO; END DO; END DO
END IF
END SUBROUTINE getBasisMappingHexa

SUBROUTINE getVolumeMapping(XCL_NGeo)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,   ONLY:NGeo,Vdm_NGeo_CLNGeo,Xi_NGeo
USE MOD_Mesh_Vars,   ONLY:nElems
USE MOD_ChangeBasis, ONLY:changeBasis3D
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,i,j,k
REAL               :: xi(3),rxi(3)
REAL               :: X_loc(3,0:NGeo,0:NGeo,0:NGeo)
!===================================================================================================================================
DO iElem=1,nElems
  X_loc=XCL_NGeo(:,:,:,:,iElem)
! Test
!  X_loc=0.
!  DO k=0,NGeo  
!    DO j=0,NGeo  
!      DO i=0,NGeo  
!        X_loc(:,i,j,k)=REAL((/i,j,k/))/REAL(NGeo)
!      END DO  !i
!    END DO  !j
!  END DO  !k
  DO k=1,NGeo-1  
    DO j=1,NGeo-1  
      DO i=1,NGeo-1
!        xi  = REAL((/i,j,k/))/REAL(NGeo) ![0,1]
        xi  = 0.5*((/Xi_NGeo(i),Xi_NGeo(j),Xi_NGeo(k)/)+1.)
        rxi = 1.-xi
        !FACES-EDGES+CORNERS
        X_loc(:,i,j,k) =  (  rxi(1)*X_loc(:,0,j,k) + xi(1)*X_loc(:,NGeo,j,k)   &
                           + rxi(2)*X_loc(:,i,0,k) + xi(2)*X_loc(:,i,NGeo,k)   &
                           + rxi(3)*X_loc(:,i,j,0) + xi(3)*X_loc(:,i,j,NGeo) ) &
                        - (  rxi(1)*rxi(2)*X_loc(:,   0,   0,k) + rxi(1)*xi(2)*X_loc(:,   0,NGeo,   k)   &
                           +  xi(1)*rxi(2)*X_loc(:,NGeo,   0,k) +  xi(1)*xi(2)*X_loc(:,NGeo,NGeo,   k)   &
                           + rxi(1)*rxi(3)*X_loc(:,   0,   j,0) + rxi(1)*xi(3)*X_loc(:,   0,   j,NGeo)   &
                           +  xi(1)*rxi(3)*X_loc(:,NGeo,   j,0) +  xi(1)*xi(3)*X_loc(:,NGeo,   j,NGeo)   &
                           + rxi(2)*rxi(3)*X_loc(:,   i,   0,0) + rxi(2)*xi(3)*X_loc(:,   i,   0,NGeo)   &
                           +  xi(2)*rxi(3)*X_loc(:,   i,NGeo,0) +  xi(2)*xi(3)*X_loc(:,   i,NGeo,NGeo) ) &
                        + (  rxi(1)*rxi(2)*rxi(3)*X_loc(:,   0,   0,   0) + rxi(1)*rxi(2)*xi(3)*X_loc(:,   0,   0,NGeo) &
                           + rxi(1)* xi(2)*rxi(3)*X_loc(:,   0,NGeo,   0) + rxi(1)* xi(2)*xi(3)*X_loc(:,   0,NGeo,NGeo) &
                           +  xi(1)*rxi(2)*rxi(3)*X_loc(:,NGeo,   0,   0) +  xi(1)*rxi(2)*xi(3)*X_loc(:,NGeo,   0,NGeo) &
                           +  xi(1)* xi(2)*rxi(3)*X_loc(:,NGeo,NGeo,   0) +  xi(1)* xi(2)*xi(3)*X_loc(:,NGeo,NGeo,NGeo) ) 
      END DO  !i
    END DO  !j
  END DO  !k
  !XCL_NGeo(:,:,:,:,iElem)=X_loc
  !from equidistant to CL
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_NGeo_CLNGeo,X_loc,XCL_NGeo(:,:,:,:,iElem))
END DO  !iElem
END SUBROUTINE getVolumeMapping


END MODULE MOD_Prepare_Mesh
