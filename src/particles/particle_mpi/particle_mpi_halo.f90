#include "boltzplatz.h"

MODULE MOD_Particle_MPI_Halo
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

#ifdef MPI
INTERFACE IdentifyHaloMPINeighborhood
  MODULE PROCEDURE IdentifyHaloMPINeighborhood
END INTERFACE

INTERFACE ExchangeHaloGeometry
  MODULE PROCEDURE ExchangeHaloGeometry
END INTERFACE

INTERFACE WriteParticlePartitionInformation
  MODULE PROCEDURE WriteParticlePartitionInformation
END INTERFACE

INTERFACE IdentifySimpleHaloMPINeighborhood
  MODULE PROCEDURE IdentifySimpleHaloMPINeighborhood
END INTERFACE


!INTERFACE WriteParticleMappingPartitionInformation
!  MODULE PROCEDURE WriteParticleMappingPartitionInformation
!END INTERFACE


PUBLIC :: IdentifyHaloMPINeighborhood,ExchangeHaloGeometry
PUBLIC :: IdentifySimpleHaloMPINeighborhood
PUBLIC :: WriteParticlePartitionInformation!,WriteParticleMappingPartitionInformation

!===================================================================================================================================

CONTAINS


SUBROUTINE IdentifyHaloMPINeighborhood(iProc,SideIndex,ElemIndex)
!===================================================================================================================================
! Searches for sides in the neighborhood of MPI-neighbors
! mark all Sides for communication which are in range of halo_eps of own MPI_Sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,         ONLY:GEO
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI
USE MOD_Particle_Surfaces_Vars,     ONLY:BezierControlPoints3D
!USE MOD_Particle_Tracking_Vars,     ONLY:DoRefMapping
USE MOD_Mesh_Vars,                  ONLY:NGeo,ElemToSide,nElems,nSides,firstMPISide_MINE,MortarSlave2MasterInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)   :: iProc  ! MPI proc with which the local proc has to exchange boundary information
INTEGER, INTENT(INOUT):: SideIndex(nSides)
INTEGER, INTENT(INOUT):: ElemIndex(PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iBGM,jBGM,kBGM,iNBProc,iElem,ElemID,ilocSide,SideID,iSide
INTEGER                     :: DataSize
! MPI exchange
TYPE tMPISideMessage
  REAL(KIND=8), ALLOCATABLE :: BezierSides3D(:,:,:,:) ! BezierSides of boundary faces
  INTEGER(KIND=4)           :: nMPISides              ! number of Sides to send
END TYPE
TYPE(tMPISideMessage)       :: SendMsg
TYPE(tMPISideMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT
!=================================================================================================================================

! 1) Exchange Sides:
!    Each proc receives all sides that lie in FIBGM cells within an  eps-distance of FIBGM cells containing 
!    any of my own MPI sides.
! 2) Go through each FIBGM cell and if there are MPI-Neighbors, search the neighbor  surrounding FIBGM
!     cells for my own MPI sides.

! Step1: find Sides of myProc that are within halo_eps distance to iProc

! here only MPI Sides.... why not INNER Sides, too????????/
! PO: correction: we send only the MPI sides, because only the MPI-Faces has to be checked, not the interior faces for
!     distance calulation. if MPI-Side is in range, then all the other sides are in range, too.
!     caution: maybe there can be a special case, in which the nBCSides are in HaloRange but not the MPI-Side
SideIndex(:)=0
SendMsg%nMPISides=0
DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
  DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
    DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
      IF (.NOT.ALLOCATED(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)) CYCLE
      IF (GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs(1).LE.0) CYCLE       ! -1 if ???
      DO iNBProc = 2,GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs(1)+1
        IF (GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs(iNBProc).EQ.iProc) THEN
          DO iElem = 1, GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
            ElemID = GEO%FIBGM(iBGM,jBGM,kBGM)%Element(iElem)
            DO iLocSide=1,6
              SideID=ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
              IF(SideID.LE.0) CYCLE
              IF(MortarSlave2MasterInfo(SideID).NE.-1) CYCLE
              IF(SideID.GE.firstMPISide_MINE)THEN
                !IF(SideID.GT.(nInnerSides+nBCSides).AND.(SideIndex(SideID).EQ.0))THEN
                ! because of implicit, but here I send for checking, other process sends the required halo region
                IF(SideIndex(SideID).EQ.0)THEN
                  SendMsg%nMPISides=SendMsg%nMPISides+1
                  SideIndex(SideID)=SendMsg%nMPISides
                END IF
              END IF
            END DO ! ilocSide
          END DO ! iElem
        END IF ! shapeProcs(i).EQ.iProc
      END DO ! iNBProc
    END DO ! iBGM
  END DO ! jBGM
END DO ! kBGM

!IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Number of Sides-To Send:   ', SendMsg%nMPISides

!--- NOTE: IF SENDMSG%NNODES IS 0 AT THIS POINT, THEN I SHOULD BE ABLE TO RETURN HERE!!!
!          This is not done yet because I'm not sure whether there are still inconsistencies in the code...

!--- Debugging information
LOGWRITE(*,*)' nMPISides for iProc=',iProc,':',SendMsg%nMPISides
!IPWRITE(UNIT_stdOut,*)' nMPISides for iProc=',iProc,':',SendMsg%nMPISides

!--- Send number of MPI sides to MPI neighbor iProc and receive number of MPI
!    sides from MPI neighbor iProc (immediate neighbor or not)
IF (PartMPI%MyRank.LT.iProc) THEN
  CALL MPI_SEND(SendMsg%nMPISides,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_RECV(RecvMsg%nMPISides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  CALL MPI_RECV(RecvMsg%nMPISides,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_SEND(SendMsg%nMPISides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
END IF

!IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Number of Sides-To Receive:', RecvMsg%nMPISides

!--- Allocate Message
IF (SendMsg%nMPISides.GT.0) THEN
  ALLOCATE(SendMsg%BezierSides3D(1:3,0:NGeo,0:nGEO,1:SendMsg%nMPISides), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMessage%BezierSides3D ',SendMsg%nMPISides)
  END IF
  SendMsg%BezierSides3D=0.
END IF
IF (RecvMsg%nMPISides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nMPISides), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMessage%BezierSides3D ',RecvMsg%nMPISides)
  END IF
  RecvMsg%BezierSides3D=0.
END IF

!--- Send any (corner-) nodes from the MPI-sides to the MPI-neighbor iProc
!    and receive iProc's (corner-) nodes in return
!--- fill send buffers
!--- Step 2: send myproc MPI-side-nodes to iProc

DO iSide=1,nSides
  IF(SideIndex(iSide).NE.0)THEN
     SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(iSide))=BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
  END IF 
END DO ! iSide=1,nSides

! CALL WriteDebugNodes(SendMsg%Nodes,SendMsg%nNodes,iProc)
!--- send and receive data
DataSize=3*(NGeo+1)*(NGeo+1)
IF(PartMPI%MyRank.LT.iProc)THEN
  IF (SendMsg%nMPISides.GT.0) CALL MPI_SEND(SendMsg%BezierSides3D       &
                                           ,SendMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1103                        &
                                           ,PartMPI%COMM                &
                                           ,IERROR                      )
  IF (RecvMsg%nMPISides.GT.0) CALL MPI_RECV(RecvMsg%BezierSides3D       &
                                           ,RecvMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1104                        &
                                           ,PartMPI%COMM                &
                                           ,MPISTATUS                   &
                                           ,IERROR                      )
ELSE IF(PartMPI%MyRank.GT.iProc)THEN
  IF (RecvMsg%nMPISides.GT.0) CALL MPI_RECV(RecvMsg%BezierSides3D       &
                                           ,RecvMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1103                        &
                                           ,PartMPI%COMM                &
                                           ,MPISTATUS                   &
                                           ,IERROR                      )
  IF (SendMsg%nMPISides.GT.0) CALL MPI_SEND(SendMsg%BezierSides3D       &
                                           ,SendMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1104                        &
                                           ,PartMPI%COMM                &
                                           ,IERROR                      )
END IF



! PO
! For each side, identifz the FIBGM cell(s) in which the side resides and
! search the surrounding nPaddingCells for neighboring elements
! For idiots: iProc and tells me which sides are on his MPI bound,
!             and I check which sides are within eps range

SideIndex(:)=0
ElemIndex(:)=0
IF (RecvMsg%nMPISides.GT.0) THEN
  CALL CheckMPINeighborhoodByFIBGM(RecvMsg%BezierSides3D,RecvMsg%nMPISides,SideIndex,ElemIndex)
END IF

!--- Deallocate Messages
IF (SendMsg%nMPISides.GT.0) THEN
  DEALLOCATE(SendMsg%BezierSides3D, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not deallocate SendMessage%BezierSides3D proc ',iProc)
  END IF
END IF
IF (RecvMsg%nMPISides.GT.0) THEN
  DEALLOCATE(RecvMsg%BezierSides3D, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not deallocate RecvMessage%BezierSides3D proc ',iProc)
  END IF
END IF

END SUBROUTINE IdentifyHaloMPINeighborhood


SUBROUTINE IdentifySimpleHaloMPINeighborhood(iProc,ElemIndex)
!===================================================================================================================================
! Searches for sides in the neighborhood of MPI-neighbors
! mark all Sides for communication which are in range of halo_eps of own MPI_Sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_MPI_Vars,                   ONLY:OffsetMPISides_send,OffsetMPISides_rec,nNbProcs,NbProc &
                                        ,nMPISides_send,nMPISides_rec,nMPISides_send,nMPISides_rec
USE MOD_Particle_Mesh_Vars,         ONLY:GEO
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,         ONLY:ElemBaryNGeo,ElemRadiusNGeo
!USE MOD_Particle_Tracking_Vars,     ONLY:DoRefMapping
USE MOD_Mesh_Vars,                  ONLY:ElemToSide,nElems,nInnerSides,nBCSides,nSides,SideToElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)   :: iProc  ! MPI proc with which the local proc has to exchange boundary information
INTEGER, INTENT(INOUT):: ElemIndex(PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iBGM,jBGM,kBGM,iNBProc,iElem,ElemID,ilocSide,SideID,SendID
INTEGER                     :: DataSize, firstSide,lastSide,iSide
! MPI exchange
TYPE tMPIElemMessage
  REAL(KIND=8), ALLOCATABLE :: ElemBaryAndRadius(:,:) ! BezierSides of boundary faces
  INTEGER(KIND=4)           :: nMPIElems              ! number of Sides to send
END TYPE
TYPE(tMPIElemMessage)       :: SendMsg
TYPE(tMPIElemMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT
!=================================================================================================================================

! 1) Exchange Sides:
!    Each proc receives all sides that lie in FIBGM cells within an  eps-distance of FIBGM cells containing 
!    any of my own MPI sides.
! 2) Go through each FIBGM cell and if there are MPI-Neighbors, search the neighbor  surrounding FIBGM
!     cells for my own MPI sides.

! Step1: find Sides of myProc that are within halo_eps distance to iProc

! here only MPI Sides.... why not INNER Sides, too????????/
! PO: correction: we send only the MPI sides, because only the MPI-Faces has to be checked, not the interior faces for
!     distance calulation. if MPI-Side is in range, then all the other sides are in range, too.
!     caution: maybe there can be a special case, in which the nBCSides are in HaloRange but not the MPI-Side
ElemIndex(:)=0
SendMsg%nMPIElems=0
! first, check if iproc is a DG-neighbor and mark direct MPI interface elements
DO iNbProc=1,nNbProcs
  IF(nbProc(iNbProc).NE.iProc) CYCLE
  ! SendID from flux 
  ! receive neighbor master sides
  ! send my Master sides
  SendID=1
  ! from MPI-Send, master sides
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    firstSide=OffsetMPISides_send(iNbProc-1,SendID)+1
    lastSide =OffsetMPISides_send(iNbProc,SendID)
    DO iSide=firstSide,lastSide
      ElemID=SideToElem(S2E_ELEM_ID,iSide)
      IF(ElemID.LT.1) STOP 'wrong'
      IF(ElemIndex(ElemID).EQ.0)THEN
        SendMsg%nMPIElems=SendMsg%nMPIElems+1
        ElemIndex(ElemID)=SendMsg%nMPIElems
      END IF
    END DO !  iSide=firstSide,lastSide
  END IF
  ! from MPI-recv, slave sides
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    firstSide=OffsetMPISides_rec(iNbProc-1,SendID)+1
    lastSide =OffsetMPISides_rec(iNbProc,SendID)
    DO iSide=firstSide,lastSide
      ElemID=SideToElem(S2E_NB_ELEM_ID,iSide)
      IF(ElemID.LT.1) STOP 'wrong'
      IF(ElemIndex(ElemID).EQ.0)THEN
        SendMsg%nMPIElems=SendMsg%nMPIElems+1
        ElemIndex(ElemID)=SendMsg%nMPIElems
      END IF
    END DO !  iSide=firstSide,lastSide
  END IF
END DO ! iNbProc=1,nNbProcs

! next, get informations from halo mesh
! get over halo sides
DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
  DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
    DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
      IF (.NOT.ALLOCATED(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)) CYCLE
      IF (GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs(1).LE.0) CYCLE       ! -1 if ???
      DO iNBProc = 2,GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs(1)+1
        IF (GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs(iNBProc).EQ.iProc) THEN
          DO iElem = 1, GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
            ElemID = GEO%FIBGM(iBGM,jBGM,kBGM)%Element(iElem)
            IF((ElemID.LT.1).OR.(ElemID.GT.PP_nElems)) CYCLE
            DO iLocSide=1,6
              SideID=ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
              IF(SideID.GT.0)THEN
                IF(SideID.GT.nSides) CYCLE
                !IF((SideID.LE.nBCSides).OR.(SideID.GT.(nBCSides+nInnerSides)))THEN
                ! PO bc sides does not have to be checked
                IF(SideID.GT.(nBCSides+nInnerSides))THEN
                  !IF(SideID.GT.(nInnerSides+nBCSides).AND.(SideIndex(SideID).EQ.0))THEN
                  ! because of implicit, but here I send for checking, other process sends the required halo region
                  IF(ElemIndex(ElemID).EQ.0)THEN
                    SendMsg%nMPIElems=SendMsg%nMPIElems+1
                    ElemIndex(ElemID)=SendMsg%nMPIElems
                  END IF
                END IF
              END IF
            END DO ! ilocSide
          END DO ! iElem
        END IF ! shapeProcs(i).EQ.iProc
      END DO ! iNBProc
    END DO ! iBGM
  END DO ! jBGM
END DO ! kBGM

!IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Number of Sides-To Send:   ', SendMsg%nMPIElems

!--- NOTE: IF SENDMSG%NNODES IS 0 AT THIS POINT, THEN I SHOULD BE ABLE TO RETURN HERE!!!
!          This is not done yet because I'm not sure whether there are still inconsistencies in the code...

!--- Debugging information
LOGWRITE(*,*)' nMPIElems for iProc=',iProc,':',SendMsg%nMPIElems
!IPWRITE(UNIT_stdOut,*)' nMPIElems for iProc=',iProc,':',SendMsg%nMPIElems

!--- Send number of MPI sides to MPI neighbor iProc and receive number of MPI
!    sides from MPI neighbor iProc (immediate neighbor or not)
IF (PartMPI%MyRank.LT.iProc) THEN
  CALL MPI_SEND(SendMsg%nMPIElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_RECV(RecvMsg%nMPIElems,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  CALL MPI_RECV(RecvMsg%nMPIElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_SEND(SendMsg%nMPIElems,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
END IF

!IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Number of Sides-To Receive:', RecvMsg%nMPIElems

!--- Allocate Message
IF (SendMsg%nMPIElems.GT.0) THEN
  ALLOCATE(SendMsg%ElemBaryAndRadius(1:4,1:SendMsg%nMPIElems), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMessage%ElemBaryAndRadius ',SendMsg%nMPIElems)
  END IF
  SendMsg%ElemBaryAndRadius=0.
END IF
IF (RecvMsg%nMPIElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemBaryAndRadius(1:4,1:RecvMsg%nMPIElems), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMessage%ElemBaryAndRadius ',RecvMsg%nMPIElems)
  END IF
  RecvMsg%ElemBaryAndRadius=0.
END IF

!--- Send any (corner-) nodes from the MPI-sides to the MPI-neighbor iProc
!    and receive iProc's (corner-) nodes in return
!--- fill send buffers
!--- Step 2: send myproc MPI-side-nodes to iProc
DO iElem=1,nElems
  IF(ElemIndex(iElem).NE.0)THEN
     SendMsg%ElemBaryAndRadius(1:3,ElemIndex(iElem))=ElemBaryNGeo(1:3,iElem)
     SendMsg%ElemBaryAndRadius( 4 ,ElemIndex(iElem))=ElemRadiusNGeo(iElem)
  END IF ! SideIndex NE.0
END DO ! iElem

! CALL WriteDebugNodes(SendMsg%Nodes,SendMsg%nNodes,iProc)
!--- send and receive data
DataSize=4
IF(PartMPI%MyRank.LT.iProc)THEN
  IF (SendMsg%nMPIElems.GT.0) CALL MPI_SEND(SendMsg%ElemBaryAndRadius   & 
                                           ,SendMsg%nMPIElems*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1103                        &
                                           ,PartMPI%COMM                &
                                           ,IERROR                      )
  IF (RecvMsg%nMPIElems.GT.0) CALL MPI_RECV(RecvMsg%ElemBaryAndRadius   &
                                           ,RecvMsg%nMPIElems*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1104                        &
                                           ,PartMPI%COMM                &
                                           ,MPISTATUS                   &
                                           ,IERROR                      )
ELSE IF(PartMPI%MyRank.GT.iProc)THEN
  IF (RecvMsg%nMPIElems.GT.0) CALL MPI_RECV(RecvMsg%ElemBaryAndRadius   &
                                           ,RecvMsg%nMPIElems*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1103                        &
                                           ,PartMPI%COMM                &
                                           ,MPISTATUS                   &
                                           ,IERROR                      )
  IF (SendMsg%nMPIElems.GT.0) CALL MPI_SEND(SendMsg%ElemBaryAndRadius   &
                                           ,SendMsg%nMPIElems*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1104                        &
                                           ,PartMPI%COMM                &
                                           ,IERROR                      )
END IF




! PO
! For each side, identifz the FIBGM cell(s) in which the side resides and
! search the surrounding nPaddingCells for neighboring elements
! For idiots: iProc and tells me which sides are on his MPI bound,
!             and I check which sides are within eps range

ElemIndex(:)=0
IF (RecvMsg%nMPIElems.GT.0) THEN
  CALL CheckSimpleMPINeighborhoodByFIBGM(RecvMsg%ElemBaryAndRadius,RecvMsg%nMPIElems,ElemIndex)
END IF

!--- Deallocate Messages
IF (SendMsg%nMPIElems.GT.0) THEN
  DEALLOCATE(SendMsg%ElemBaryAndRadius, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not deallocate SendMessage%ElemBaryAndRadius proc ',iProc)
  END IF
END IF
IF (RecvMsg%nMPIElems.GT.0) THEN
  DEALLOCATE(RecvMsg%ElemBaryAndRadius, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not deallocate RecvMessage%ElemBaryAndRadius proc ',iProc)
  END IF
END IF

END SUBROUTINE IdentifySimpleHaloMPINeighborhood


SUBROUTINE CheckSimpleMPINeighborhoodByFIBGM(ElemBaryAndRadius,nExternalElems,ElemIndex)
!===================================================================================================================================
! Compute distance of MPI-Side to my-Sides, if distance below helo_eps2 distance, mark size for MPI-Exchange
! Question: Why does one does not mark the INNER Sides, too??
!           Then, the halo region would only require the elements inside of itself and not only the MPI sides...??
!           Or is it sufficient?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,        ONLY:GEO, FIBGMCellPadding,NbrOfCases,casematrix
USE MOD_Particle_MPI_Vars,         ONLY:halo_eps
USE MOD_Particle_Mesh_Vars,        ONLY:ElemBaryNGeo,ElemRadiusNGeo
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,    INTENT(IN)      :: ElemBaryAndRadius(1:4,1:nExternalElems)
INTEGER, INTENT(IN)      :: nExternalElems
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)   :: ElemIndex(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem, ElemID,iBGMElem,iCase,NbOfElems
INTEGER                  :: iBGMmin,iBGMmax,jBGMmin,jBGMmax,kBGMmin,kBGMmax,iPBGM,jPBGM,kPBGM
REAL                     :: Vec1(1:3),Vec2(1:3),Vec3(1:3), NodeX(1:3), Radius,Distance,Vec0(1:3)
REAL                     :: xmin, xmax, ymin,ymax,zmin,zmax, Radius2
LOGICAL                  :: IsChecked(1:PP_nElems)
!===================================================================================================================================

! For each (NGeo+1)^2 BezierControlPoint of each side, the FIBGM cell(s) in which the side 
! resides is identified and the surrouding nPaddingCells for each neighboring element
! are searched
!--- for idiots: get BezierControlPoints of myProc that are within eps distance to MPI-bound 
!                of iProc
ElemIndex=0
NBOfElems=0

DO iElem=1,nExternalElems
  IsChecked=.FALSE.
  ! check only the min-max values
  ! get elem extension based on barycenter and radius
  NodeX(1:3) = ElemBaryAndRadius(1:3,iElem) 
  Radius     = ElemBaryAndRadius( 4 ,iElem)

  xmin = NodeX(1) -Radius
  ymin = NodeX(2) -Radius
  zmin = NodeX(3) -Radius
  xmax = NodeX(1) +Radius
  ymax = NodeX(2) +Radius
  zmax = NodeX(3) +Radius
  Radius2=Radius*Radius

  ! BGM mesh cells
  iBGMmin = CEILING((xMin-GEO%xminglob)/GEO%FIBGMdeltas(1))-FIBGMCellPadding(1)
  iBGMmax = CEILING((xMax-GEO%xminglob)/GEO%FIBGMdeltas(1))+FIBGMCellPadding(1)

  jBGMmin = CEILING((yMin-GEO%yminglob)/GEO%FIBGMdeltas(2))-FIBGMCellPadding(2)
  jBGMmax = CEILING((yMax-GEO%yminglob)/GEO%FIBGMdeltas(2))+FIBGMCellPadding(2)

  kBGMmin = CEILING((zMin-GEO%zminglob)/GEO%FIBGMdeltas(3))-FIBGMCellPadding(3)
  kBGMmax = CEILING((zMax-GEO%zminglob)/GEO%FIBGMdeltas(3))+FIBGMCellPadding(3)

  iBGMmin  = MAX(GEO%FIBGMimin,iBGMmin)
  jBGMmin  = MAX(GEO%FIBGMjmin,jBGMmin)
  kBGMmin  = MAX(GEO%FIBGMkmin,kBGMmin)

  iBGMmax  = MIN(GEO%FIBGMimax,iBGMmax)
  jBGMmax  = MIN(GEO%FIBGMjmax,jBGMmax)
  kBGMmax  = MIN(GEO%FIBGMkmax,kBGMmax)

  ! loop over BGM cells    
  DO iPBGM = iBGMmin,iBGMmax
    DO jPBGM = jBGMmin,jBGMmax
      DO kPBGM = kBGMmin,kBGMmax
        IF(.NOT.ALLOCATED(GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element))CYCLE
        !IF((iPBGM.GT.GEO%FIBGMimax).OR.(iPBGM.LT.GEO%FIBGMimin) .OR. &
        !   (jPBGM.GT.GEO%FIBGMjmax).OR.(jPBGM.LT.GEO%FIBGMjmin) .OR. &
        !   (kPBGM.GT.GEO%FIBGMkmax).OR.(kPBGM.LT.GEO%FIBGMkmin) ) CYCLE
        DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
          ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
          !IF((ElemID.LE.0.).OR.(ElemID.GT.PP_nElems))CYCLE
          IF(IsChecked(ElemID)) THEN
            CYCLE
          ELSE
            IF(ElemIndex(ElemID).EQ.0)THEN
              Vec0=ElemBaryNGeo(1:3,ElemID)-ElemBaryAndRadius(1:3,iElem)
              !Distance=SQRT(DOT_PRODUCT(Vec0,Vec0)) &
              !        -Radius-ElemRadiusNGeo(ElemID)
              Distance=SQRT(DOT_PRODUCT(Vec0,Vec0)) &
                      -Radius-ElemRadiusNGeo(ElemID)
              IF((Distance.LE.halo_eps).OR.(Distance.LT.0.))THEN
                NbOfElems=NbOfElems+1
                ElemIndex(ElemID)=NbofElems
              END IF ! in range
              IsChecked(ElemID)=.TRUE.
            END IF ! ElemIndex(ElemID).EQ.0
          END IF
        END DO ! iBGMElem
      END DO ! kPBGM
    END DO ! jPBGM
  END DO ! i PBGM
END DO ! iElem


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
  !--- check sides shifted by periodic vectors, add to SideIndex if match
  DO iElem=1,nExternalElems
    DO iCase = 1, NbrOfCases
      IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
          (casematrix(iCase,2).EQ.0) .AND. &
          (casematrix(iCase,3).EQ.0)) CYCLE
     ! DO iNode=1,8
        NodeX =ElemBaryAndRadius(1:3,iElem)
        Radius=ElemBaryAndRadius( 4 ,iElem)
        Radius2=Radius*Radius
        !SELECT CASE(iNode)
        NodeX(:) = NodeX(:) + &
                   casematrix(iCase,1)*Vec1(1:3) + &
                   casematrix(iCase,2)*Vec2(1:3) + &
                   casematrix(iCase,3)*Vec3(1:3) 


        xmin = NodeX(1) -Radius
        ymin = NodeX(2) -Radius
        zmin = NodeX(3) -Radius
        xmax = NodeX(1) +Radius
        ymax = NodeX(2) +Radius
        zmax = NodeX(3) +Radius

        ! BGM mesh cells
        iBGMmin = CEILING((xMin-GEO%xminglob)/GEO%FIBGMdeltas(1))-FIBGMCellPadding(1)
        iBGMmax = CEILING((xMax-GEO%xminglob)/GEO%FIBGMdeltas(1))+FIBGMCellPadding(1)
      
        jBGMmin = CEILING((yMin-GEO%yminglob)/GEO%FIBGMdeltas(2))-FIBGMCellPadding(2)
        jBGMmax = CEILING((yMax-GEO%yminglob)/GEO%FIBGMdeltas(2))+FIBGMCellPadding(2)
      
        kBGMmin = CEILING((zMin-GEO%zminglob)/GEO%FIBGMdeltas(3))-FIBGMCellPadding(3)
        kBGMmax = CEILING((zMax-GEO%zminglob)/GEO%FIBGMdeltas(3))+FIBGMCellPadding(3)
      
        iBGMmin  = MIN(MAX(GEO%FIBGMimin,iBGMmin),GEO%FIBGMimax)
        jBGMmin  = MIN(MAX(GEO%FIBGMjmin,jBGMmin),GEO%FIBGMjmax)
        kBGMmin  = MIN(MAX(GEO%FIBGMkmin,kBGMmin),GEO%FIBGMkmax)
      
        iBGMmax  = MAX(MIN(GEO%FIBGMimax,iBGMmax),GEO%FIBGMimin)
        jBGMmax  = MAX(MIN(GEO%FIBGMjmax,jBGMmax),GEO%FIBGMjmin)
        kBGMmax  = MAX(MIN(GEO%FIBGMkmax,kBGMmax),GEO%FIBGMkmin)

        ! loop over BGM cells    
        DO iPBGM = iBGMmin,iBGMmax
          DO jPBGM = jBGMmin,jBGMmax
            DO kPBGM = kBGMmin,kBGMmax
              IF(.NOT.ALLOCATED(GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element))CYCLE
              DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
                ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
                IF((ElemID.LT.0).OR.(ElemID.GT.PP_nElems))CYCLE
                IF(ElemIndex(ElemID).EQ.0)THEN
                  Vec0=ElemBaryNGeo(1:3,ElemID)-NodeX(1:3)
                  Distance=SQRT(DOT_PRODUCT(Vec0,Vec0)) &
                          -Radius-ElemRadiusNGeo(ElemID)
                  IF(Distance.LE.halo_eps)THEN
                    NbOfElems=NbOfElems+1
                    ElemIndex(ElemID)=NbofElems
                  END IF ! in range
                END IF ! ElemIndex(ElemID).EQ.0
              END DO ! iBGMElem
            END DO ! kPBGM
          END DO ! jPBGM
        END DO ! iPBGM
     ! END DO ! Node=1,8
    END DO ! iCase
  END DO ! iElem
END IF  ! nperiodicvectors>0

END SUBROUTINE CheckSimpleMPINeighborhoodByFIBGM


SUBROUTINE CheckMPINeighborhoodByFIBGM(BezierSides3D,nExternalSides,SideIndex,ElemIndex)
!===================================================================================================================================
! Compute distance of MPI-Side to my-Sides, if distance below helo_eps2 distance, mark size for MPI-Exchange
! Question: Why does one does not mark the INNER Sides, too??
!           Then, the halo region would only require the elements inside of itself and not only the MPI sides...??
!           Or is it sufficient?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,                 ONLY:NGeo,ElemToSide,nSides,SideToElem,MortarSlave2MasterInfo,ElemToElemGlob,OffSetElem
USE MOD_Mesh_Vars,                 ONLY:MortarType,MortarInfo
USE MOD_Particle_Mesh_Vars,        ONLY:GEO, FIBGMCellPadding,NbrOfCases,casematrix
USE MOD_Particle_MPI_Vars,         ONLY:halo_eps
USE MOD_Particle_Surfaces_Vars,    ONLY:BezierControlPoints3D
!USE MOD_Particle_Tracking_Vars,    ONLY:DoRefMapping

!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,    INTENT(IN)      :: BezierSides3D(1:3,0:NGeo,0:NGeo,1:nExternalSides)
INTEGER, INTENT(IN)      :: nExternalSides
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)   :: SideIndex(nSides)
INTEGER, INTENT(INOUT)   :: ElemIndex(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iSide, NbOfSides,p,q,ElemID,ilocSide,SideID,r,s,iBGMElem,iCase,NbOfElems,NBElemID,iMortar,locMortarType
REAL                     :: NodeX(1:3)!,xNodes(1:3,0:NGeo,0:NGeo)
INTEGER                  :: iBGM,jBGM,kBGM,iPBGM,jPBGM,kPBGM
LOGICAL                  :: leave
REAL                     :: Vec1(1:3),Vec2(1:3),Vec3(1:3)
!===================================================================================================================================

! For each (NGeo+1)^2 BezierControlPoint of each side, the FIBGM cell(s) in which the side 
! resides is identified and the surrouding nPaddingCells for each neighboring element
! are searched
!--- for idiots: get BezierControlPoints of myProc that are within eps distance to MPI-bound 
!                of iProc
SideIndex(:)=0
ElemIndex=0
NbOfSides=0
NBOfElems=0

DO iSide=1,nExternalSides
  !IF(MortarSlave2MasterInfo(iSide).NE.-1) CYCLE
  DO q=0,NGeo
    DO p=0,NGeo
      NodeX(:) = BezierSides3D(:,p,q,iSide)
      !iBGM = INT((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      !jBGM = INT((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      !kBGM = INT((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      iBGM = CEILING((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
      iBGM = MIN(GEO%FIBGMimax,iBGM)                             
      jBGM = CEILING((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
      jBGM = MIN(GEO%FIBGMjmax,jBGM) 
      kBGM = CEILING((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
      kBGM = MIN(GEO%FIBGMkmax,kBGM)
      DO iPBGM = iBGM-FIBGMCellPadding(1),iBGM+FIBGMCellPadding(1)
        DO jPBGM = jBGM-FIBGMCellPadding(2),jBGM+FIBGMCellPadding(2)
          DO kPBGM = kBGM-FIBGMCellPadding(3),kBGM+FIBGMCellPadding(3)
            IF((iPBGM.GT.GEO%FIBGMimax).OR.(iPBGM.LT.GEO%FIBGMimin) .OR. &
               (jPBGM.GT.GEO%FIBGMjmax).OR.(jPBGM.LT.GEO%FIBGMjmin) .OR. &
               (kPBGM.GT.GEO%FIBGMkmax).OR.(kPBGM.LT.GEO%FIBGMkmin) ) CYCLE
            DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
              ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
              IF(ElemIndex(ElemID).GT.0) CYCLE ! element is already marked
              DO ilocSide=1,6
                SideID=ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
                ! caution, not save if corect
                leave=.FALSE.
                IF(SideIndex(SideID).EQ.0)THEN
                  DO s=0,NGeo
                    DO r=0,NGeo
                      IF(SQRT(DOT_Product(BezierControlPoints3D(:,r,s,SideID)-NodeX &
                                         ,BezierControlPoints3D(:,r,s,SideID)-NodeX )).LE.halo_eps)THEN
                        NbOfSides=NbOfSides+1
                        SideIndex(SideID)=NbOfSides
                        IF(ElemIndex(ElemID).EQ.0)THEN
                          NbOfElems=NbOfElems+1
                          ElemIndex(ElemID)=NbofElems
                        END IF
                        leave=.TRUE.
                        ! mark potential Inner elements on the other side
                        DO iMortar=1,4
                          CHECKSAFEINT(ElemToElemGlob(iMortar,ilocSide,offSetElem+ElemID),4)
                          NBElemID=ElemToElemGlob(iMortar,ilocSide,offSetElem+ElemID)-offSetElem
                          IF(NBElemID.LE.0) CYCLE
                          IF(NBElemID.GT.PP_nElems) CYCLE
                          ! check if NBElem is already marked, if not, mark it!
                          IF(ElemIndex(NbElemID).EQ.0)THEN
                            NbOfElems=NbOfElems+1
                            ElemIndex(NbElemID)=NbofElems
                          END IF
                        END DO ! iMortar=1,4
                        EXIT
                      END IF
                    END DO ! r
                    IF(leave) EXIT
                  END DO ! s
                END IF ! SideIndex(SideID).EQ.0
                IF(leave) EXIT ! Element is marked, all other sides of THIS element have not to be checked
              END DO ! ilocSide
            END DO ! iBGMElem
          END DO ! kPBGM
        END DO ! jPBGM
      END DO ! i PBGM
    END DO ! p
  END DO ! q
END DO ! iSide

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
  !--- check sides shifted by periodic vectors, add to SideIndex if match
  DO iSide=1,nExternalSides
    DO iCase = 1, NbrOfCases
      IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
          (casematrix(iCase,2).EQ.0) .AND. &
          (casematrix(iCase,3).EQ.0)) CYCLE
      DO q=0,NGeo
        DO p=0,NGeo
          NodeX(:) = BezierSides3d(:,p,q,iSide) + &
                     casematrix(iCase,1)*Vec1(1:3) + &
                     casematrix(iCase,2)*Vec2(1:3) + &
                     casematrix(iCase,3)*Vec3(1:3) 
          iBGM = INT((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          jBGM = INT((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          kBGM = INT((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          DO iPBGM = iBGM-FIBGMCellPadding(1),iBGM+FIBGMCellPadding(1)
            DO jPBGM = jBGM-FIBGMCellPadding(2),jBGM+FIBGMCellPadding(2)
              DO kPBGM = kBGM-FIBGMCellPadding(3),kBGM+FIBGMCellPadding(3)
                IF((iPBGM.GT.GEO%FIBGMimax).OR.(iPBGM.LT.GEO%FIBGMimin) .OR. &
                   (jPBGM.GT.GEO%FIBGMjmax).OR.(jPBGM.LT.GEO%FIBGMjmin) .OR. &
                   (kPBGM.GT.GEO%FIBGMkmax).OR.(kPBGM.LT.GEO%FIBGMkmin) ) CYCLE
                DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
                  ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
                  IF(ElemIndex(ElemID).GT.0) CYCLE ! element is already marked
                  DO ilocSide=1,6
                    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
                    ! caution, not save if corect
                    leave=.FALSE.
                    IF(SideIndex(SideID).EQ.0)THEN
                      DO s=0,NGeo
                        DO r=0,NGeo
                          IF(SQRT(DOT_Product(BezierControlPoints3D(:,r,s,SideID)-NodeX &
                                             ,BezierControlPoints3D(:,r,s,SideID)-NodeX )).LE.halo_eps)THEN
                            NbOfSides=NbOfSides+1
                            SideIndex(SideID)=NbOfSides
                            IF(ElemIndex(ElemID).EQ.0)THEN
                              NbOfElems=NbOfElems+1
                              ElemIndex(ElemID)=NbofElems
                            END IF
                            leave=.TRUE.
                            EXIT
                          END IF
                        END DO ! r
                        IF(leave) EXIT
                      END DO ! s
                    END IF ! SideIndex(SideID).EQ.0
                    IF(leave) EXIT
                  END DO ! ilocSide
                END DO ! iElem
              END DO ! kPBGM
            END DO ! jPBGM
          END DO ! i PBGM
        END DO ! p
      END DO ! q
    END DO ! iCase
  END DO ! iSide
END IF  ! nperiodicvectors>0

! finally, all elements connected to this side have to be marked
! this is a sanity step and could be omitted
DO iSide=1,nSides
  IF(SideIndex(iSide).GT.0)THEN
    ! check both elements connected to side if they are marked
    ! master
    ElemID=SideToElem(S2E_ELEM_ID,iSide)
    IF((ElemID.GT.0).AND.(ElemID.LE.PP_nElems))THEN
      IF(ElemIndex(ElemID).EQ.0)THEN
        NbOfElems=NbOfElems+1
        ElemIndex(ElemID)=NbofElems
      END IF
    END IF
    ! slave
    ElemID=SideToElem(S2E_NB_ELEM_ID,iSide)
    IF((ElemID.GT.0).AND.(ElemID.LE.PP_nElems))THEN
      IF(ElemIndex(ElemID).EQ.0)THEN
        NbOfElems=NbOfElems+1
        ElemIndex(ElemID)=NbofElems
      END IF
    END IF
  END IF
END DO ! iSide=1,nSides

!IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Number of marked sides:   ', NbOfSides

END SUBROUTINE CheckMPINeighborhoodByFIBGM


!SUBROUTINE ExchangeHaloGeometry(iProc,SideList,ElemList)
SUBROUTINE ExchangeHaloGeometry(iProc,ElemList)
!===================================================================================================================================
! exchange of halo geometry
! including:
! BezierControlPoints3D
! ElemToSide
! BC-Type and State
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartHaloElemToProc
USE MOD_Mesh_Vars,              ONLY:nElems, nSides, nBCSides, nInnerSides, ElemToSide, BC,nGeo,SideToElem
USE MOD_Mesh_Vars,              ONLY:MortarSlave2MasterInfo,OffSetElem
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalSides,nTotalElems,SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElemGlob,GEO,nTotalBCSides,ElemBaryNGeo
!USE MOD_Particle_Surfaces_Vars, ONLY:ElemSlabNormals,ElemSlabIntervals  
USE MOD_Mesh_Vars,              ONLY:XCL_NGeo,dXCL_NGeo,MortarType
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicDisplacement,PartElemToElemGlob
! should not be needed annymore
!USE MOD_Particle_MPI_Vars,      ONLY:nNbProcs,offsetMPISides_MINE, offsetMPISides_YOUR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iProc       ! MPI proc with which the local proc is to exchange boundary information
!INTEGER, INTENT(INOUT)          :: SideList(nSides)
INTEGER, INTENT(INOUT)          :: ElemList(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE tMPISideMessage
  REAL,ALLOCATABLE          :: BezierControlPoints3D(:,:,:,:)
  REAL,ALLOCATABLE          :: XCL_NGeo (:,:,:,:,:)
  REAL,ALLOCATABLE          :: dXCL_NGeo(:,:,:,:,:,:)
!  REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ElemSlabNormals    
!  REAL,ALLOCATABLE,DIMENSION(:,:)     :: ElemSlabIntervals   
  REAL,ALLOCATABLE,DIMENSION(:,:)     :: ElemBaryNGeo   
  INTEGER,ALLOCATABLE       :: ElemToSide(:,:,:) 
  INTEGER(KIND=8),ALLOCATABLE ::ElemToElemGlob(:,:,:) 
  INTEGER,ALLOCATABLE       :: SideToElem(:,:)
  INTEGER,ALLOCATABLE       :: MortarType(:,:)
  INTEGER,ALLOCATABLE       :: SideBCType(:)
  INTEGER,ALLOCATABLE       :: MortarSlave2MasterInfo(:)
  INTEGER,ALLOCATABLE       :: BC(:)
  INTEGER,ALLOCATABLE       :: NativeElemID(:)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: SideSlabNormals                  ! normal vectors of bounding slab box
  REAL,ALLOCATABLE,DIMENSION(:,:)    :: SideSlabIntervals               ! intervalls beta1, beta2, beta3
  LOGICAL,ALLOCATABLE,DIMENSION(:)   :: BoundingBoxIsEmpty
  !INTEGER,ALLOCATABLE       :: PeriodicElemSide(:,:)
  INTEGER                   :: nSides                 ! number of sides to send
  INTEGER                   :: nElems                 ! number of elems to send
END TYPE
TYPE(tMPISideMessage)       :: SendMsg
TYPE(tMPISideMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT
INTEGER                     :: newSideID,haloSideID,ioldSide,oldElemID,newElemID
REAL                        :: xNodes(1:3,0:NGeo,0:NGeo)
REAL                        :: xNodes2(1:3,0:NGeo,0:NGeo)
LOGICAL                     :: isDoubleSide
LOGICAL,ALLOCATABLE         :: isElem(:),isSide(:),isDone(:)
INTEGER, ALLOCATABLE        :: ElemIndex(:), SideIndex(:),HaloInc(:)
INTEGER                     :: iElem, ilocSide,SideID,iSide,iIndex,iHaloSide,flip
INTEGER                     :: nDoubleSides,nDoubleBezier,tmpnSides,tmpnElems
INTEGER                     :: datasize,datasize2,datasize3,ElemIDGlob
INTEGER                     :: p,q,tmpbcsides
INTEGER                     :: ElemID,ElemID2,hostElemId,idisplace,locsideid,newbcsideid
!===================================================================================================================================

ALLOCATE(isElem(1:nElems))
IF (.NOT.ALLOCATED(isElem)) CALL abort(&
  __STAMP__&
  ,'Could not allocate isElem')
isElem(:) = .FALSE.

ALLOCATE(isSide(1:nSides))
IF (.NOT.ALLOCATED(isSide)) CALL abort(&
  __STAMP__&
  ,'Could not allocate isSide')
isSide(:) = .FALSE.

ALLOCATE(ElemIndex(1:nElems))
IF (.NOT.ALLOCATED(ElemIndex)) CALL abort(&
  __STAMP__&
  ,'Could not allocate ElemIndex')
ElemIndex(:) = 0

ALLOCATE(SideIndex(1:nSides))
IF (.NOT.ALLOCATED(SideIndex)) CALL abort(&
  __STAMP__&
  ,'Could not allocate SideIndex')
SideIndex(:) = 0

!--- First, count marker node indices (nNeighborhoodNodes are within eps distance of at least one MPI-neighbor's node)
!--- For each MPI neighbor, identify the number of sides and elements to be sent
SendMsg%nElems=0
SendMsg%nSides=0
!LOGWRITE(*,*)'nNeighborhoodNodes=',nNeighborhoodNodes

! 1) get number of elements and sides
DO iElem=1,nElems
  IF(ElemList(iElem).NE.0)THEN
    !IF(.NOT.isElem(iElem)) THEN
    SendMsg%nElems=SendMsg%nElems+1
    ElemIndex(iElem) = SendMsg%nElems
    isElem(iElem)=.TRUE.
    !END IF ! NOT isElem
  END IF
END DO ! iElem

! 2) mark all required sides and get number of send sides
DO iElem=1,nElems
  IF(isElem(iElem))THEN
    DO ilocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      IF(DoRefMapping)THEN
        IF(.NOT.isSide(SideID).AND.SideID.LE.nBCSides) THEN
          ! missing: what do do with BC sides??"
          SendMsg%nSides=SendMsg%nSides+1
          SideIndex(SideID) = SendMsg%nSides
          isSide(SideID)=.TRUE.
        END IF ! not isSide
      ELSE
        IF(.NOT.isSide(SideID)) THEN
          SendMsg%nSides=SendMsg%nSides+1
          SideIndex(SideID) = SendMsg%nSides
          isSide(SideID)=.TRUE.
        END IF ! not isSide
      END IF
    END DO ! ilocSide
  END IF ! Element is marked to send
END DO ! iElem

!IPWRITE(*,*) 'iproc,sideindex,nsides',SideIndex(:),SendMsg%nSides
!IPWRITE(*,*) 'iproc,nsides',iProc,SendMsg%nSides,nTotalSides,nTotalBCSides

!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!WRITE(*,*) "Nodes:", SendMsg%nNodes,"Sides:",SendMsg%nSides,"elems:",SendMsg%nElems,"iProc",PMPIVAR%iProc
!--- Communicate number of sides (trias,quads), elems (tets,hexas) and nodes to each MPI proc
  
IF (PartMPI%MyRank.LT.iProc) THEN
  CALL MPI_SEND(SendMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_SEND(SendMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
  CALL MPI_RECV(RecvMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_RECV(RecvMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  CALL MPI_RECV(RecvMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_RECV(RecvMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_SEND(SendMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_SEND(SendMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
END IF

! allocate send buffers for nodes, sides and elements for each MPI neighbor
! ElemToSide Mapping 
IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%ElemToSide(1:2,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate SendMsg%ElemToSide',SendMsg%nElems)
  SendMsg%ElemToSide(:,:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemToSide(1:2,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
  RecvMsg%ElemToSide(:,:,:)=0
END IF

! ElemToElemGlob Mapping 
IF (SendMsg%nElems.GT.0) THEN       ! ElemToElem(1:4,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%ElemToElemGlob(1:4,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate SendMsg%ElemToElemGlob',SendMsg%nElems)
  SendMsg%ElemToElemGlob(:,:,:)=-1
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemToElemGlob(1:4,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%ElemToElemGlob',RecvMsg%nElems)
  RecvMsg%ElemToElemGlob(:,:,:)=-1
END IF
! BezierControlPoints3D for exchange
IF (SendMsg%nSides.GT.0) THEN       ! Beziercontrolpoints3d
  ALLOCATE(SendMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%BezierControlPoints3D',SendMsg%nSides)
  SendMsg%BezierControlPoints3D=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%BezierControlPoints3D',RecvMsg%nSides)
  RecvMsg%BezierControlPoints3D=0.
END IF

IF(DoRefMapping)THEN
  ! XCL_NGeo for exchange
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%XCL_NGeo',SendMsg%nElems)
    SendMsg%XCL_NGeo(:,:,:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%XCL_NGeo',RecvMsg%nElems)
    RecvMsg%XCL_NGeo(:,:,:,:,:)=0
  END IF
  
  ! DXCL_NGeo for exchange
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%DXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%DXCL_NGeo',SendMsg%nElems)
    SendMsg%DXCL_NGeo(:,:,:,:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%DXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
    RecvMsg%DXCL_NGeo(:,:,:,:,:,:)=0
  END IF

!  ! ElemSlabNormals for exchange
!  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
!    ALLOCATE(SendMsg%ElemSlabNormals(1:3,0:3,1:SendMsg%nElems),STAT=ALLOCSTAT) 
!    IF (ALLOCSTAT.NE.0) CALL abort(&
!     __STAMP__&
!      ,'Could not allocate SendMsg%ElemSlabNormals',SendMsg%nElems)
!    SendMsg%ElemSlabNormals(:,:,:)=0
!  END IF
!  IF (RecvMsg%nElems.GT.0) THEN
!    ALLOCATE(RecvMsg%ElemSlabNormals(1:3,0:3,1:RecvMsg%nElems),STAT=ALLOCSTAT)  
!    IF (ALLOCSTAT.NE.0) CALL abort(&
!      __STAMP__&
!      ,'Could not allocate RecvMsg%ElemSlabNormals',RecvMsg%nElems)
!    RecvMsg%ElemSlabNormals(:,:,:)=0
!  END IF
  
!  ! ElemSlabIntervals for exchange
!  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
!    ALLOCATE(SendMsg%ElemSlabIntervals(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT) 
!    IF (ALLOCSTAT.NE.0) CALL abort(&
!      __STAMP__&
!      ,'Could not allocate SendMsg%ElemSlabIntervals',SendMsg%nElems)
!    SendMsg%DXCL_NGeo(:,:,:,:,:,:)=0
!  END IF
!  IF (RecvMsg%nElems.GT.0) THEN
!    ALLOCATE(RecvMsg%ElemSlabIntervals(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT) 
!    IF (ALLOCSTAT.NE.0) CALL abort(&
!      __STAMP__&
!      ,'Could not allocate RecvMsg%ElemSlabIntervals',RecvMsg%nElems)
!    RecvMsg%DXCL_NGeo(:,:,:,:,:,:)=0
!  END IF
END IF

! ElemBaryNGeo
IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%ElemBaryNGeo(1:3,1:SendMsg%nElems),STAT=ALLOCSTAT) 
  IF (ALLOCSTAT.NE.0) CALL abort(&
   __STAMP__&
    ,'Could not allocate SendMsg%ElemBaryNGeo',SendMsg%nElems)
  SendMsg%ElemBaryNGeo(:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemBaryNGeo(1:3,1:RecvMsg%nElems),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%ElemBary',RecvMsg%nElems)
  RecvMsg%ElemBaryNGeo(:,:)=0
END IF

! SideToElem Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SideToElem(1:2,1:nSides) 
  ALLOCATE(SendMsg%SideToElem(1:5,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideToElem',SendMsg%nSides)
  SendMsg%SideToElem(:,:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideToElem(1:5,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideToElem',RecvMsg%nSides)
  RecvMsg%SideToElem(:,:)=0
END IF

! BC Mapping
! BC(1:4,1:nSides) 
! 1:BC, 
! 2:NBProc, 
! 3:1=Mine/2=Yours 
! 4:SIDE_ID-MPI_Offset(NBProc)  
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%BC(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%BC',SendMsg%nSides,999.)
  SendMsg%BC(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BC(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%BC',RecvMsg%nSides,999.)
  RecvMsg%BC(:)=0
END IF
! mortartype
IF(.NOT.DoRefMapping)THEN
  IF (SendMsg%nSides.GT.0) THEN       
    ALLOCATE(SendMsg%MortarType(1:2,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%MortarType',SendMsg%nSides,999.)
    SendMsg%MortarType(:,:)=0
  END IF
  IF (RecvMsg%nSides.GT.0) THEN
    ALLOCATE(RecvMsg%MortarType(1:2,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%MortarType',RecvMsg%nSides,999.)
    RecvMsg%MortarType(:,:)=0
  END IF
END IF
! SideBCType
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%SideBCType(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideBCType',SendMsg%nSides,999.)
  SendMsg%SideBCType(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideBCType(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideBCType',RecvMsg%nSides,999.)
  RecvMsg%SideBCType(:)=0
END IF
! NativeElemID 
IF (SendMsg%nElems.GT.0) THEN 
  ALLOCATE(SendMsg%NativeElemID(1:SendMsg%nElems),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%NativeElemID',SendMsg%nElems,999.)
  SendMsg%NativeElemID(:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%NativeElemID(1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%NativeElemID',RecvMsg%nElems,999.)
  RecvMsg%NativeElemID(:)=0
END IF
! SideSlabNormals Mapping
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%SideSlabNormals(1:3,1:3,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideSlabNormals',SendMsg%nSides)
  SendMsg%SideSlabNormals(:,:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideSlabNormals(1:3,1:3,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideSlabNormals',RecvMsg%nSides)
  RecvMsg%SideSlabNormals(:,:,:)=0.
END IF
! SideSlabIntervals Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SideSlabIntervals(1:2,1:nSides) 
  ALLOCATE(SendMsg%SideSlabIntervals(1:6,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideSlabIntervals',SendMsg%nSides)
  SendMsg%SideSlabIntervals(:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideSlabIntervals(1:6,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideSlabIntervals',RecvMsg%nSides)
  RecvMsg%SideSlabIntervals(:,:)=0.
END IF
! BoundingBoxIsEmpty Mapping
IF (SendMsg%nSides.GT.0) THEN       ! BoundingBoxIsEmpty(1:2,1:nSides) 
  ALLOCATE(SendMsg%BoundingBoxIsEmpty(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%BoundingBoxIsEmpty',SendMsg%nSides)
  SendMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BoundingBoxIsEmpty(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%BoundingBoxIsEmpty',RecvMsg%nSides)
  RecvMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
!! PeriodicElemSide Mapping 
!IF (SendMsg%nElems.GT.0) THEN       ! PeriodicElemSide(1:iLocSide,1:nElems)
!  ALLOCATE(SendMsg%PeriodicElemSide(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
!   ,'Could not allocate SendMsg%PeriodicElemSide',SendMsg%nElems,999.)
!  SendMsg%PeriodicElemSide(:,:)=0
!END IF
!IF (RecvMsg%nElems.GT.0) THEN
!  ALLOCATE(RecvMsg%PeriodicElemSide(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
!   ,'Could not allocate RecvMsg%PeriodicElemSide',RecvMsg%nElems,999.)
!  RecvMsg%PeriodicElemSide(:,:)=0
!END IF


! fill send buffers with node, side and element data (including connectivity!)
! ElemtoSide 
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
    IF(DoRefMapping)THEN
      SendMsg%XCL_NGeo(:,:,:,:,ElemIndex(iElem))=XCL_NGeo(:,:,:,:,iElem)
      SendMsg%dXCL_NGeo(:,:,:,:,:,ElemIndex(iElem))=dXCL_NGeo(:,:,:,:,:,iElem)
      !SendMsg%ElemSlabNormals(:,:,ElemIndex(iElem))=ElemSlabNormals(:,:,iElem)
      !SendMsg%ElemSlabIntervals(:,ElemIndex(iElem))=ElemSlabIntervals(:,iElem)
    END IF
    SendMsg%ElemBaryNGeo(:,ElemIndex(iElem))=ElemBaryNGeo(:,iElem)
    DO iLocSide = 1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      IF(SideIndex(SideID).GT.0)THEN
        SendMsg%ElemToSide(1,iLocSide,ElemIndex(iElem)) = &
                 SideIndex(ElemToSide(E2S_SIDE_ID,iLocSide,iElem))
             ! CAUTION DEBUG correct sideid????
        SendMsg%ElemToSide(2,iLocSide,ElemIndex(iElem)) = &
                ElemToSide(2,iLocSide,iElem)
      END IF
    END DO
  END IF
END DO

! part of elemtoelemglob which is needed
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
    SendMsg%ElemToElemGlob(1:4,1:6,ElemIndex(iElem)) = &
            PartElemToElemGlob(1:4,1:6,iElem)
        IF(Myrank.EQ.1) WRITE(*,*) 'sending...',PartElemToElemGlob(1:4,1:6,iElem)
  END IF
END DO

! SideToElem Mapping & BezierControloints3D
DO iSide = 1,nSides
  IF (SideIndex(iSide).GT.0) THEN
  !IF(isSide(iSide))THEN
    DO iIndex = 1,2   ! S2E_ELEM_ID, S2E_NB_ELEM_ID
      IF (SideToElem(iIndex,iSide).GT.0) THEN
         SendMsg%SideToElem(iIndex,SideIndex(iSide)) = &
                 ElemIndex(SideToElem(iIndex,iSide))     
         SendMsg%SideBCType(SideIndex(iSide)) = SidePeriodicType(iSide)
      END IF
    END DO ! S2E_LOC_SIDE_ID, S2E_NB_LOC_SIDE_ID, S2E_FLIP
    SendMsg%SideToElem(3:5,SideIndex(iSide)) = &
        SideToElem(3:5,iSide)
    SendMsg%BezierControlPoints3D(:,:,:,SideIndex(iSide)) = &
        BezierControlPoints3D(:,:,:,iSide)
    ! slabnormals
    SendMsg%SideSlabNormals(:,:,SideIndex(iSide)) = &
        SideSlabNormals(:,:,iSide)
    ! slabintervalls
    SendMsg%SideSlabIntervals(:,SideIndex(iSide)) = &
        SideSlabIntervals(:,iSide)
    ! BoundingBoxIsEmpty
    SendMsg%BoundingBoxIsEmpty(SideIndex(iSide)) = &
        BoundingBoxIsEmpty(iSide)
  END IF
END DO
!--- BC Mapping ------------------------------------------------------!
DO iSide = 1,nSides  ! no need to go through all side since BC(1:nBCSides)
  IF (SideIndex(iSide).NE.0) THEN
    SendMsg%BC(SideIndex(iSide)) = BC(iSide)
    IF(.NOT.DoRefMapping) SendMsg%MortarType(:,SideIndex(iSide))=MortarType(:,iSide)
  END IF
END DO

!DO iSide=1,SendMsg%nSides
!  IF(SUM(ABS(SendMsg%SideSlabIntervals(:,iSide))).EQ.0)THEN
!    IPWRITE(*,*) 'error sending slabnormal',iSide
!  END IF
!END DO


!! CAUTION & DEBUG Do we need this????
!! do we have to do this???
!DO iSide = nBCSides+nInnerSides+1,nSides ! only go through MPI-Sides 
!  IF(SideIndex(iSide).NE.0) THEN
!    DO iNbProc =1,nNbProcs  ! ??? iNbProc=1 start at nBCSides+nInnerSides ???
!      ! first check "mine"-sides:
!      SendMsg%BC(1,SideIndex(iSide)) = -1 
!      IF ((iSide.GE.offsetMPISides_MINE(iNbProc-1)+1).AND. &
!          (iSide.LE.offsetMPISides_MINE(iNbProc))) THEN
!          SendMsg%BC(2,SideIndex(iSide)) = NbProc(iNbProc)  ! global NbProcID
!          SendMsg%BC(3,SideIndex(iSide)) = 1                ! 1=Mine/2=Yours  
!          SendMsg%BC(4,SideIndex(iSide)) = iSide-offsetMPISides_MINE(iNbProc-1) ! MPITag 
!      ! then check "your"-sides:
!      ELSE IF ((iSide.GE.offsetMPISides_YOUR(iNbProc-1)+1).AND. &
!          (iSide.LE.offsetMPISides_YOUR(iNbProc))) THEN
!          SendMsg%BC(2,SideIndex(iSide)) = NbProc(iNbProc)  ! global NbProcID
!          SendMsg%BC(3,SideIndex(iSide)) = 2                ! 1=Mine/2=Yours  
!          SendMsg%BC(4,SideIndex(iSide)) = iSide-offsetMPISides_YOUR(iNbProc-1) ! MPITag 
!      END IF
!    END DO
!  END IF
!END DO
!! proceed MPI sides
!DO iSide = nBCSides+1,nBCSides+nInnerSides
!  IF(SideIndex(iSide).NE.0) THEN
!    IF ((SendMsg%SideToElem(1,SideIndex(iSide)).EQ.0).OR. &
!        (SendMsg%SideToElem(2,SideIndex(iSide)).EQ.0)) THEN
!      SendMsg%BC(1,SideIndex(iSide)) = 424242        ! EPS-Region-Bound in innersides
!    END IF
!  END IF
!END DO
! NativeElemID 
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
    SendMsg%NativeElemID(ElemIndex(iElem)) = iElem
  END IF
END DO 
!! PeriodicElemSide 
!! DEBUG CAUTION should be altered to this
! IF(GEO%nPeriodicVectors.GT.0)THEN
!   DO iElem = 1,nElems
!     IF (ElemIndex(iElem).NE.0) THEN
!       DO iLocSide = 1,6
!         SendMsg%PeriodicElemSide(iLocSide,ElemIndex(iElem)) = &
!              NeighborElemID(ilocSide,iElem)
!              !GEO%PeriodicElemSide(iLocSide,iElem)
!       END DO
!     END IF
!   END DO
! END IF


!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!IPWRITE(UNIT_stdOut,*) " Now MPI exchange"

dataSize=3*(NGeo+1)*(NGeo+1)
dataSize2=3*(NGeo+1)*(NGeo+1)*(NGeo+1)
dataSize3=9*(NGeo+1)*(NGeo+1)*(NGeo+1)

IF (PartMPI%MyRank.LT.iProc) THEN
  ! Send:
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BezierControlPoints3D,SendMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,IERROR)
!  IF(GEO%nPeriodicVectors.GT.0)THEN
!    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%PeriodicElemSide,SendMsg%nElems*6,MPI_INTEGER ,iProc,1109,PartMPI%COMM,IERROR)
!  END IF
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideBCType,SendMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabNormals,SendMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabIntervals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
  IF(DoRefMapping)THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)
!    IF (SendMsg%nElems.GT.0) &
!        CALL MPI_SEND(SendMsg%ElemSlabNormals,SendMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,IERROR)
!    IF (SendMsg%nElems.GT.0) &
!        CALL MPI_SEND(SendMsg%ElemSlabIntervals,SendMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemBaryNGeo,SendMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToElemGlob,SendMsg%nElems*24,MPI_LONG       ,iProc,1119,PartMPI%COMM,IERROR)
  IF(.NOT.DoRefMapping) THEN
    IF (SendMsg%nSides.GT.0) &
        CALL MPI_SEND(SendMsg%MortarType,SendMsg%nSides*2,MPI_INTEGER       ,iProc,1120,PartMPI%COMM,IERROR)
  END IF

  ! Receive:
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
  CALL MPI_RECV(RecvMsg%BezierControlPoints3D,RecvMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,MPISTATUS,IERROR)

  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,MPISTATUS,IERROR)
  !IF(GEO%nPeriodicVectors.GT.0)THEN
  !  IF (RecvMsg%nElems.GT.0) &
  !       CALL MPI_RECV(RecvMsg%PeriodicElemSide,RecvMsg%nElems*6,MPI_INTEGER,iProc,1109,PartMPI%COMM,MPISTATUS,IERROR)
  !END IF
  IF (RecvMsg%nSides.GT.0) CALL MPI_RECV(RecvMsg%SideBCType,RecvMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabNormals,RecvMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabIntervals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  IF(DoRefMapping)THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
!    IF (RecvMsg%nElems.GT.0) &
!        CALL MPI_RECV(RecvMsg%ElemSlabNormals,RecvMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,MPISTATUS,IERROR)
!    IF (RecvMsg%nElems.GT.0) &
!        CALL MPI_RECV(RecvMsg%ElemSlabIntervals,RecvMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemBaryNGeo,RecvMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemToElemGlob,RecvMsg%nElems*24,MPI_LONG,iProc,1119,PartMPI%COMM,MPISTATUS,IERROR)
  IF(.NOT.DoRefMapping) THEN
    IF (RecvMsg%nSides.GT.0) &
        CALL MPI_RECV(RecvMsg%MortarType,RecvMsg%nSides*2,MPI_INTEGER,iProc,1120,PartMPI%COMM,MPISTATUS,IERROR)
  END IF


ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  ! Receive:
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
  CALL MPI_RECV(RecvMsg%BezierControlPoints3D,RecvMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,MPISTATUS,IERROR)
  !IF(GEO%nPeriodicVectors.GT.0)THEN
  !  IF (RecvMsg%nElems.GT.0) &
  !       CALL MPI_RECV(RecvMsg%PeriodicElemSide,RecvMsg%nElems*6,MPI_INTEGER,iProc,1109,PartMPI%COMM,MPISTATUS,IERROR)
  !END IF
  IF (RecvMsg%nSides.GT.0) CALL MPI_RECV(RecvMsg%SideBCType,RecvMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabNormals,RecvMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabIntervals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  IF(DoRefMapping)THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
!    IF (RecvMsg%nElems.GT.0) &
!        CALL MPI_RECV(RecvMsg%ElemSlabNormals,RecvMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,MPISTATUS,IERROR)
!    IF (RecvMsg%nElems.GT.0) &
!        CALL MPI_RECV(RecvMsg%ElemSlabIntervals,RecvMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemBaryNGeo,RecvMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemToElemGlob,RecvMsg%nElems*24,MPI_LONG,iProc,1119,PartMPI%COMM,MPISTATUS,IERROR)
  IF(.NOT.DoRefMapping) THEN
    IF (RecvMsg%nSides.GT.0) &
        CALL MPI_RECV(RecvMsg%MortarType,RecvMsg%nSides*2,MPI_INTEGER,iProc,1120,PartMPI%COMM,MPISTATUS,IERROR)
  END IF

  ! Send:
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BezierControlPoints3D,SendMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,IERROR)
  !IF(GEO%nPeriodicVectors.GT.0)THEN
  !  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%PeriodicElemSide,SendMsg%nElems*6,MPI_INTEGER ,iProc,1109,PartMPI%COMM,IERROR)
  !END IF
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideBCType,SendMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabNormals,SendMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabIntervals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)

  IF(DoRefMapping)THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)
!    IF (SendMsg%nElems.GT.0) &
!        CALL MPI_SEND(SendMsg%ElemSlabNormals,SendMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,IERROR)
!    IF (SendMsg%nElems.GT.0) &
!        CALL MPI_SEND(SendMsg%ElemSlabIntervals,SendMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemBaryNGeo,SendMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemToElemGlob,SendMsg%nElems*24,MPI_LONG,iProc,1119,PartMPI%COMM,IERROR)
  IF(.NOT.DoRefMapping) THEN
    IF (SendMsg%nSides.GT.0) &
        CALL MPI_SEND(SendMsg%MortarType,SendMsg%nSides*2,MPI_INTEGER       ,iProc,1120,PartMPI%COMM,IERROR)
  END IF
END IF

IF ((RecvMsg%nElems.EQ.0) .AND. (RecvMsg%nSides.GT.0))THEN
    ERRWRITE(*,*)'ERROR: nElems=0 when nSides=',RecvMsg%nSides,' and nSides=',RecvMsg%nSides,'!'
    CALL abort(&
     __STAMP__&
     ,'nElems=0 while nSides=',RecvMsg%nSides)
END IF

DEALLOCATE(isElem,isSide,ElemIndex,SideIndex)

!DO iSide=1,RecvMsg%nSides
!  IF(SUM(ABS(RecvMsg%SideSlabIntervals(:,iSide))).EQ.0)THEN
!    IPWRITE(*,*) 'error recv slabnormal',iSide
!  END IF
!END DO


!IPWRITE(UNIT_stdOut,*) " Recive stuff"
IF(DoRefMapping)THEN
  IF (RecvMsg%nElems.GT.0) THEN
    ! now, the famous reconstruction of geometry
    ! add the halo region to the existing geometry
    ! therefore, the PartElemToSide,... has to be extended
    ! multiple sides ( MPI-Sides) should be ignored
    
    ! get number of double sides
    ! BezierControlPoints are build in master system, therefore, the node indicies are unique and 0,0 and NGeo,NGeo are on diagonal,
    ! opposide sides of the side
    nDoubleSides=0
    ALLOCATE(isSide(1:RecvMsg%nSides))
    ALLOCATE(isDone(1:RecvMsg%nSides))
    isDone=.FALSE.
    isSide=.TRUE.
    ! get increament for each halo side
    ! 1) get increment of halo side id
    ! HaloSideID is increment of new side
    ALLOCATE(HaloInc(1:RecvMsg%nSides))
    HaloInc=0
    HaloSideID=0
    DO iHaloSide=1,RecvMsg%nSides
      IF(isSide(iHaloSide))THEN
        HaloSideID=HaloSideID+1
        HaloInc(iHaloSide)=HaloSideID
      END IF
    END DO ! iHaloSide
    
    ! new number of sides
    !IPWRITE(*,*) 'nTotalSides,ntotBCSides,ntotalelems',nTotalSides,nTotalBCSides,nTotalElems
    tmpnSides    =nTotalSides
    tmpnElems    =nTotalElems
    tmpBCSides   =nTotalBCSides
    nTotalSides  =nTotalSides+RecvMsg%nSides-nDoubleSides
    nTotalBCSides=nTotalBCSides+RecvMsg%nSides-nDoubleSides
    nTotalElems  =nTotalElems+RecvMsg%nElems
    !tmpnSides    =nTotalSides
    !IPWRITE(*,*) 'NewnTotalSides,Newntotalbcsides,Newntotalelesm',nTotalSides,nTotalBCSides,nTotalElems
    CALL ResizeParticleMeshData(tmpnSides,tmpnElems,nTotalSides,nTotalElems,tmpBCSides,nTotalBCSides)

    ! loop over all elements and add them
    DO iElem=1,RecvMsg%nElems
      newElemID=tmpnElems+iElem
      PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
      PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
      XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
      dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
      ElemBaryNGeo(1:3,newElemID) = RecvMsg%ElemBaryNGeo(1:3,iElem)
    END DO

    ! loop over all sides and add them
    DO iSide=1,RecvMsg%nSides
      newSideID  =tmpnSides+iSide
      ! BezierControlPoints and BoundingBox stuff
      BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,iSide)
      ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
      SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,iSide)
      SideSlabIntervals(1:6,newSideID) =RecvMsg%SideSlabIntervals(1:6 ,iSide) 
      BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( iSide) 
      ! BC type, etc. 
      BC(newSideID)=RecvMsg%BC(iSide)
      SidePeriodicType(newSideID)=RecvMsg%SideBCType(iSide)
      PartBCSideList(newSideID)=tmpBCSides+iSide 
    END DO 

    ! fill lists 
    ! caution: PartSideToElem is only filled for BC sides
    DO iElem=1,RecvMsg%nElems
      newElemID=tmpnElems+iElem
      DO ilocSide=1,6
        SideID=RecvMsg%ElemToSide(1,iLocSide,iElem)
        IF(SideID.GT.0)THEN
          ! fill PartElemToSide
          newSideID=tmpnSides+SideID
          PartElemToSide(E2S_SIDE_ID,iLocSide,NewElemID)=newSideID
          PartElemToSide(E2S_FLIP,iLocSide,NewElemID)=0
          ! and SideToElem
          PartSideToElem(S2E_ELEM_ID,newSideID)=newElemID
        END IF
      END DO 
      ! list from ElemToElemGlob mapped to process local element
      ! new list points from local-elem-id to global
      PartElemToElemGlob(1:4,1:6,newElemID) = RecvMsg%ElemToElemGlob(1:4,1:6,iElem)
    END DO

    IF(.NOT.PartMPI%isMPINeighbor(iProc))THEN
      PartMPI%isMPINeighbor(iProc) = .true.
      PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
    END IF
    DEALLOCATE(isSide)
    DEALLOCATE(isDone)
    DEALLOCATE(HaloInc)
    DO iSide=1,nTotalSides
      IF(SUM(ABS(SideSlabNormals(:,:,iSide))).EQ.0)THEN
        CALL abort(&
__STAMP__&
          ,' SideSlabNormals is zero!,iSide',iSide)
      END IF
    END DO 
  END IF ! RecvMsg%nSides>0
ELSE ! DoRefMappping=F
 IF (RecvMsg%nSides.GT.0) THEN
    ! now, the famous reconstruction of geometry
    ! add the halo region to the existing geometry
    ! therefore, the PartElemToSide,... has to be extended
    ! multiple sides ( MPI-Sides) should are duplicated to deal easier with 
    ! mortar faces, the elemtoelemandside lists are reconstructed via the global element list
    ! this loop is not optimized, hence, the MPI-faces exists twice, 
    ! once for the master and the slave element
    
    nDoubleSides=0
    ALLOCATE(isSide(1:RecvMsg%nSides))
    ALLOCATE(isDone(1:RecvMsg%nSides))
    isDone=.FALSE.
    isSide=.TRUE.
    ! get increament for each halo side
    ! 1) get increment of halo side id
    ! HaloSideID is increment of new side
    ALLOCATE(HaloInc(1:RecvMsg%nSides))
    HaloInc=0
    HaloSideID=0
    DO iHaloSide=1,RecvMsg%nSides
      IF(isSide(iHaloSide))THEN
        HaloSideID=HaloSideID+1
        HaloInc(iHaloSide)=HaloSideID
      END IF
    END DO ! iHaloSide
    
    tmpnSides =nTotalSides
    tmpnElems=nTotalElems
    nTotalSides=nTotalSides+RecvMsg%nSides-nDoubleSides
    nTotalElems=nTotalElems+RecvMsg%nElems
    CALL ResizeParticleMeshData(tmpnSides,tmpnElems,nTotalSides,nTotalElems)
  
    !DO iElem=tmpnElems+1,nTotalElems
    DO iElem=1,RecvMsg%nElems
      !print*,'iElem',iElem
      ! first, new SideID=entry of RecvMsg+tmpnSides
      newElemID=tmpnElems+iElem
      DO ilocSide=1,6
        haloSideID=RecvMsg%ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
        ! first, set new sideid
      !  print*,'haloSideId',haloSideID
        IF(isSide(haloSideID)) THEN
          IF(HaloInc(haloSideID).EQ.0) IPWRITE(UNIT_stdOut,*) ' Warning: wrong halo inc'
          newSideID=tmpnSides+haloinc(haloSideID)
          IF(newSideID.LT.tmpnSides) IPWRITE(UNIT_stdOut,*) 'Warning: wrong new sideid', newsideid
        END IF
        ! build PartSideToElem, so much as possible
        ! get correct side out of RecvMsg%SideToElem
        IF(iElem.EQ.RecvMsg%SideToElem(S2E_ELEM_ID,haloSideID))THEN
          PartSideToElem(S2E_ELEM_ID,newSideID)    =newElemID
          PartSideToElem(S2E_LOC_SIDE_ID,newSideID)=ilocSide
          PartSideToElem(S2E_FLIP       ,newSideID)=RecvMsg%SideToElem(S2E_FLIP,haloSideID)
        ELSE IF(iElem.EQ.RecvMsg%SideToElem(S2E_NB_ELEM_ID,haloSideID))THEN
          PartSideToElem(S2E_NB_ELEM_ID,newSideID) =newElemID
          PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) =ilocSide
          PartSideToElem(S2E_FLIP      ,newSideID) =RecvMsg%SideToElem(S2E_FLIP,haloSideID)
        ELSE ! should be found, because there should be halo sides without any connection
            CALL abort(&
 __STAMP__&
          ,'Non-Critical error in domain reconstrution. IF NOT encountered, something is terrible wrong.')
        END IF
        !BC(1:4,newSideID)=RecvMsg%BC(1:4,haloSideID)
        IF(.NOT.isDone(haloSideID))THEN
          BC(newSideID)=RecvMsg%BC(haloSideID)
          SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
  
          ! copy Bezier to new side id
          BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
          ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
          SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
          SideSlabIntervals(1:6,newSideID)  =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
          BoundingBoxIsEmpty(newSideID)     =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
          MortarType(1:2,newSideID)         =RecvMsg%MortarType(1:2,haloSideID)
          isDone(haloSideID)=.TRUE.
        END IF
        ! build entry to PartElemToSide
        PartElemToSide(1,iLocSide,newElemId)=newSideID
        PartElemToSide(2,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
      END DO ! ilocSide
      ! set native elemID
      PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
      PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
      ElemBaryNGeo(1:3,newElemID) = RecvMsg%ElemBaryNGeo(1:3,iElem)
      ! list from ElemToElemGlob mapped to process local element
      ! new list points from local-elem-id to global
      SWRITE(*,*) '--------------------------newElemID,iElem',newElemID,iElem
      SWRITE(*,*) '--------------------------',RecvMsg%ElemToElemGlob(1:4,1:6,iElem)
      PartElemToElemGlob(1:4,1:6,newElemID) = RecvMsg%ElemToElemGlob(1:4,1:6,iElem)
    END DO ! iElem
    ! build rest: PartElemToElem, PartLocSideID
    DO iElem=PP_nElems+1,nTotalElems
      DO ilocSide=1,6
        flip   = PartElemToSide(E2S_FLIP,ilocSide,iElem)
        SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        ! check of sideid
      END DO ! ilocSide
    END DO ! Elem
    IF(.NOT.PartMPI%isMPINeighbor(iProc))THEN
      PartMPI%isMPINeighbor(iProc) = .true.
      PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
    END IF
    DEALLOCATE(isSide)
    DEALLOCATE(isDone)
    DEALLOCATE(HaloInc)
  END IF ! RecvMsg%nSides>0
END IF

END SUBROUTINE ExchangeHaloGeometry


SUBROUTINE ResizeParticleMeshData(nOldSides,nOldElems,nTotalSides,nTotalElems,nOldBCSides,nTotalBCSides)
!===================================================================================================================================
! resize the partilce mesh data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartHaloElemToProc
USE MOD_Mesh_Vars,              ONLY:BC,nGeo,nElems,XCL_NGeo,DXCL_NGEO,OffSetElem,MortarType
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElemGlob,ElemBaryNGeo
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
!USE MOD_Particle_Surfaces_Vars, ONLY:ElemSlabNormals,ElemSlabIntervals  
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: nOldSides,nOldElems,nTotalSides,nTotalElems
INTEGER,INTENT(IN),OPTIONAL        :: nOldBCSides,nTotalBCSides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: ALLOCSTAT,nLower
INTEGER,ALLOCATABLE                :: DummyElemToSide(:,:,:)                                
INTEGER,ALLOCATABLE                :: DummyBC(:)                                
REAL,ALLOCATABLE                   :: DummyBezierControlPoints3D(:,:,:,:)                                
REAL,ALLOCATABLE                   :: DummyXCL_NGEO (:,:,:,:,:)                                
REAL,ALLOCATABLE                   :: DummydXCL_NGEO(:,:,:,:,:,:)                                
REAL,ALLOCATABLE                   :: DummyElemBaryNGeo(:,:)                                
INTEGER,ALLOCATABLE                :: DummyHaloToProc(:,:)                                 
INTEGER,ALLOCATABLE                :: DummyMortarType(:,:)                                 
INTEGER,ALLOCATABLE                :: DummySideToElem(:,:)
INTEGER,ALLOCATABLE                :: DummySideBCType(:),DummyPartBCSideList(:)
INTEGER(KIND=8),ALLOCATABLE        :: DummyElemToElem(:,:,:)
REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummySideSlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummySideSlabIntervals               ! intervalls beta1, beta2, beta3
REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummyElemSlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummyElemSlabIntervals               ! intervalls beta1, beta2, beta3
LOGICAL,ALLOCATABLE,DIMENSION(:)   :: DummyBoundingBoxIsEmpty
!===================================================================================================================================

! reallocate shapes
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!IPWRITE(*,*) 'oldSides,oldElems,totalSides,TotalElems,oldBCSides,bcsides',nOldSides,nOldElems,nTotalSides &
!                                                                          ,ntotalElems,noldBCSides,nTotalBCSides
!print*,'rank, in reshape',myrank
!print*,'rank, in reshape',myrank, nOldSides,nOldElems,nTotalSides,nTotalElems
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
! PartElemToSide
ALLOCATE(DummyElemToSide(1:2,1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemToSide)) CALL abort(&
  __STAMP__&
  ,'Could not allocate ElemIndex')
DummyElemToSide=PartElemToSide
!IPWRITE(UNIT_stdOut,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
DEALLOCATE(PartElemToSide)
ALLOCATE(PartElemToSide(1:2,1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'Could not allocate PartElemToSide')
PartElemToSide=-1
PartElemToSide(:,:,1:nOldElems) =DummyElemToSide(:,:,1:nOldElems)
DEALLOCATE(DummyElemToSide)

!! ElemToElemGlob
!ALLOCATE(DummyElemToElemGlob(1:4,1:6,1:nOldElems))
!IF (.NOT.ALLOCATED(DummyElemToSide)) CALL abort(&
!  __STAMP__&
!  ,'Could not allocate ElemIndex')
!DummyElemToElemGlob(:,:,1:nOldElems)=ElemToElemGlob(:,:,offsetElem+1:offSetElem+nOldElems)
!!IPWRITE(UNIT_stdOut,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
!DEALLOCATE(ElemToElemGlob)
!ALLOCATE(ElemToElemGlob(1:4,1:6,offSetElem+1:offSetElem+nTotalElems),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!  __STAMP__&
!  ,'Could not allocate PartElemToSide')
!ElemToElemGlob=-1
!ElemToElemGlob(:,:,offSetElem+1:offSetElem+nOldElems) =DummyElemToElemGlob(:,:,1:nOldElems)
!DEALLOCATE(DummyElemToSide)

IF(DoRefMapping)THEN
  ! XCL_NGeo
  ALLOCATE(DummyXCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyXCL_NGeo)) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  DummyXCL_NGeo=XCL_NGeo
  !IPWRITE(UNIT_stdOut,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
  DEALLOCATE(XCL_NGeo)
  ALLOCATE(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate XCL_NGeo')
  XCL_NGeo=0.
  XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems) =DummyXCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems)
  DEALLOCATE(DummyXCL_NGeo)
  ! dXCL_NGeo
  ALLOCATE(DummydXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems))
  IF (.NOT.ALLOCATED(DummydXCL_NGeo)) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  DummydXCL_NGeo=dXCL_NGeo
  !IPWRITE(UNIT_stdOut,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
  DEALLOCATE(dXCL_NGeo)
  ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate dXCL_NGeo')
  dXCL_NGeo=0.
  dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems) =DummydXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems)
  DEALLOCATE(DummydXCL_NGeo)
  ! PartBCSideList
  ALLOCATE(DummyPartBCSideList(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyPartBCSideList)) CALL abort(&
    __STAMP__&
    ,'Could not allocate PartBCSideList')
  DummyPartBCSideList(1:nOldSides)=PartBCSideList(1:nOldSides)
  DEALLOCATE(PartBCSideList)
  ALLOCATE(PartBCSideList(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  PartBCSideList=-1 !HUGE(1)
  PartBCSideList(1:nOldSides) =DummyPartBCSideList(1:nOldSides)
  DEALLOCATE(DummyPartBCSideList)
!  ! ElemSlabNormals
!  ALLOCATE(DummyElemSlabNormals(1:3,0:3,1:nOldElems))
!  IF (.NOT.ALLOCATED(DummyElemSlabNormals)) CALL abort(&
!    __STAMP__&
!    ,'Could not allocate ElemIndex')
!  DummyElemSlabNormals=ElemSlabNormals
!  DEALLOCATE(ElemSlabNormals)
!  ALLOCATE(ElemSlabNormals(1:3,0:3,1:nTotalElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(&
!    __STAMP__&
!    ,'Could not allocate ElemIndex')
!  ElemSlabNormals=0
!  ElemSlabNormals(1:3,0:3,1:nOldElems) =DummyElemSlabNormals(1:3,0:3,1:nOldElems)
!  DEALLOCATE(DummyElemSlabNormals)
!  ! ElemSlabIntervals
!  ALLOCATE(DummyElemSlabIntervals(1:6,1:nOldElems))
!  IF (.NOT.ALLOCATED(DummyElemSlabIntervals)) CALL abort(&
!    __STAMP__&
!    ,'Could not allocate ElemIndex')
!  DummyElemSlabIntervals=ElemSlabIntervals
!  DEALLOCATE(ElemSlabIntervals)
!  ALLOCATE(ElemSlabIntervals(1:6,1:nTotalElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(&
!    __STAMP__&
!    ,'Could not allocate ElemIndex')
!  ElemSlabIntervals=0
!  ElemSlabIntervals(1:6,1:nOldElems) =DummyElemSlabIntervals(1:6,1:nOldElems)
!  DEALLOCATE(DummyElemSlabIntervals)
END IF

!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!print*,' done elem to side',myrank
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
! HaloToProc
IF(.NOT.ALLOCATED(PartHaloElemToProc))THEN
  !print*,'both here',myrank,nElems,nElems+1,nTotalElems
  !print*,'myrank',myrank,allocstat
  nLower=nElems+1
  ALLOCATE(PartHaloElemToProc(1:3,nLower:nTotalElems),STAT=ALLOCSTAT)                                 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate PartHaloElemToProc')
  PartHaloElemToProc=-1
  !print*,'lower,upper',PP_nElems+1,nTotalElems
ELSE
  nLower=nElems+1
  ALLOCATE(DummyHaloToProc(1:3,nLower:nOldElems))                                 
  IF (.NOT.ALLOCATED(DummyHaloToProc)) CALL abort(&
    __STAMP__&
    ,'Could not allocate DummyPartHaloElemToProc')
  DummyHaloToProc=PartHaloElemToProc
  DEALLOCATE(PartHaloElemToProc)
  ALLOCATE(PartHaloElemToProc(1:3,nLower:nTotalElems),STAT=ALLOCSTAT)                                 
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate PartHaloElemToProc')
  ! copy array to new
  PartHaloElemToProc=-1
  PartHaloElemToProc(1:3,PP_nElems+1:nOldElems)    =DummyHaloToProc(1:3,PP_nElems+1:nOldElems)
  DEALLOCATE(DummyHaloToProc)
END IF
!print*,' done halotoproc',myrank
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!print*,' done halotoproc',myrank
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
! PartSideToElem
ALLOCATE(DummySideToElem(1:5,1:nOldSides))
IF (.NOT.ALLOCATED(DummySideToElem)) CALL abort(&
  __STAMP__&
  ,'Could not allocate ElemIndex')
DummySideToElem=PartSideToElem
DEALLOCATE(PartSideToElem)
ALLOCATE(PartSideToElem(1:5,1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'Could not allocate PartSideToElem')
PartSideToElem=-1
PartSideToElem(1:5,1:nOldSides  )              =DummySideToElem(1:5,1:nOldSides)
DEALLOCATE(DummySideToElem)
!print*,' done side to elem',myrank
! PartElemToElemGlob
ALLOCATE(DummyElemToElem(1:4,1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemToElem)) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
DummyElemToElem=PartElemToElemGlob
DEALLOCATE(PartElemToElemGlob)
ALLOCATE(PartElemToElemGlob(1:4,1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
PartElemToElemGlob=-1
PartElemToElemGlob(1:4,1:6,1:nOldElems)            =DummyElemToElem(1:4,1:6,1:nOldElems)
DEALLOCATE(DummyElemToElem)
IF(DoRefMapping)THEN
  ! BezierControlPoints3D
  ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummyBezierControlPoints3d)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummyBezierControlPoints3d=BezierControlPoints3d
  DEALLOCATE(BezierControlPoints3D)
  ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  BezierControlPoints3d(:,:,:,1:nOldBCSides) =DummyBezierControlPoints3D(:,:,:,1:nOldBCSides)
  DEALLOCATE(DummyBezierControlPoints3D)
  ! SideSlabNormals
  ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySideSlabNormals)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummySideSlabNormals=SideSlabNormals
  DEALLOCATE(SideSlabNormals)
  ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  SideSlabNormals=0
  SideSlabNormals(1:3,1:3,1:nOldBCSides) =DummySideSlabNormals(1:3,1:3,1:nOldBCSides)
  DEALLOCATE(DummySideSlabNormals)
  ! SideSlabIntervals
  ALLOCATE(DummySideSlabIntervals(1:6,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySideSlabIntervals)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummySideSlabIntervals=SideSlabIntervals
  DEALLOCATE(SideSlabIntervals)
  ALLOCATE(SideSlabIntervals(1:6,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  SideSlabIntervals=0
  SideSlabIntervals(1:6,1:nOldBCSides) =DummySideSlabIntervals(1:6,1:nOldBCSides)
  DEALLOCATE(DummySideSlabIntervals)
  ! BoundingBoxIsEmpty
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
  DEALLOCATE(BoundingBoxIsEmpty)
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  BoundingBoxIsEmpty(1:nOldBCSides) =DummyBoundingBoxIsEmpty(1:nOldBCSides)
  DEALLOCATE(DummyBoundingBoxIsEmpty)

ELSE ! no mapping
  ! BezierControlPoints3D
  ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBezierControlPoints3d)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummyBezierControlPoints3d=BezierControlPoints3d
  DEALLOCATE(BezierControlPoints3D)
  ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  BezierControlPoints3d(:,:,:,1:nOldSides) =DummyBezierControlPoints3D(:,:,:,1:nOldSides)
  DEALLOCATE(DummyBezierControlPoints3D)
  ! SideSlabNormals
  ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideSlabNormals)) CALL abort(&
    __STAMP__&
   ,'Could not allocate ElemIndex')
  DummySideSlabNormals=SideSlabNormals
  DEALLOCATE(SideSlabNormals)
  ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  SideSlabNormals=0
  SideSlabNormals(1:3,1:3,1:nOldSides) =DummySideSlabNormals(1:3,1:3,1:nOldSides)
  DEALLOCATE(DummySideSlabNormals)
  ! SideSlabIntervals
  ALLOCATE(DummySideSlabIntervals(1:6,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideSlabIntervals)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummySideSlabIntervals=SideSlabIntervals
  DEALLOCATE(SideSlabIntervals)
  ALLOCATE(SideSlabIntervals(1:6,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  SideSlabIntervals=0
  SideSlabIntervals(1:6,1:nOldSides) =DummySideSlabIntervals(1:6,1:nOldSides)
  DEALLOCATE(DummySideSlabIntervals)
  ! BoundingBoxIsEmpty
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
  DEALLOCATE(BoundingBoxIsEmpty)
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  BoundingBoxIsEmpty(1:nOldSides) =DummyBoundingBoxIsEmpty(1:nOldSides)
  DEALLOCATE(DummyBoundingBoxIsEmpty)
END IF

! ElemBaryNGeo
ALLOCATE(DummyElemBaryNGeo(1:3,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemBaryNGeo)) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemBaryNGeo')
DummyElemBaryNGeo=ElemBaryNGeo
DEALLOCATE(ElemBaryNGeo)
ALLOCATE(ElemBaryNGeo(1:3,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemBaryNGeo')
ElemBaryNGeo=0.
ElemBaryNGeo(1:3,1:nOldElems) =DummyElemBaryNGeo(1:3,1:nOldElems)

! SideBCType
ALLOCATE(DummySideBCType(1:nOldSides))
IF (.NOT.ALLOCATED(DummySideBCType)) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
DummySideBCType(1:nOldSides)=SidePeriodicType(1:nOldSides)
DEALLOCATE(SidePeriodicType)
ALLOCATE(SidePeriodicType(1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
SidePeriodicType=-1
SidePeriodicType(1:nOldSides) =DummySideBCType(1:nOldSides)
DEALLOCATE(DummySideBCType)
! BC
ALLOCATE(DummyBC(1:nOldSides))
IF (.NOT.ALLOCATED(DummyBC)) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
! check
!IF(ALLOCATED(BC)) print*,'yes it is'
!print*,'size',size(BC)
!print*, BC
DummyBC(1:nOldSides)=BC(1:nOldSides)
DEALLOCATE(BC)
ALLOCATE(BC(1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
BC=0
BC(1:nOldSides) =DummyBC(1:nOldSides)
DEALLOCATE(DummyBC)
! finished copying
IF(.NOT.DoRefMapping)THEN
  ! MortarType
  ALLOCATE(DummyMortarType(1:2,1:nOldSides))
  IF (.NOT.ALLOCATED(DummyMortarType)) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  ! check
  DummyMortarType(1:2,1:nOldSides)=MortarType(1:2,1:nOldSides)
  DEALLOCATE(MortarType)
  ALLOCATE(MortarType(1:2,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate ElemIndex')
  MortarType=-1
  MortarType(1:2,1:nOldSides) =DummyMortarType(1:2,1:nOldSides)
  DEALLOCATE(DummyMortarType)
END IF

END SUBROUTINE ResizeParticleMeshData


SUBROUTINE WriteParticlePartitionInformation(nPlanar,nBilinear,nCurved,nPlanarHalo,nBilinearHalo,nCurvedHalo &
                                        ,nBCElems,nLinearElems,nCurvedElems,nBCElemsHalo,nLinearElemsHalo,nCurvedElemsHalo)
!===================================================================================================================================
! write the particle partition information to file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,              ONLY:nSides,nElems
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalSides,nTotalElems
USE MOD_LoadBalance_Vars,       ONLY:DoLoadBalance,nLoadBalance, writePartitionInfo
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              ::nPlanar,nBilinear,nCurved,nPlanarHalo,nBilinearHalo,nCurvedHalo
INTEGER,INTENT(IN)              ::nBCElems,nLinearElems,nCurvedElems,nBCElemsHalo,nLinearElemsHalo,nCurvedElemsHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=10)          :: formatstr
INTEGER,ALLOCATABLE        :: nNBProcs_glob(:), ProcInfo_glob(:,:),NBInfo_glob(:,:), NBInfo(:), tmparray(:,:),ProcInfo(:)
REAL,ALLOCATABLE           :: tmpreal(:,:)
INTEGER                    :: nVars,nNBmax,i,j,ioUnit
CHARACTER(LEN=64)          :: filename
CHARACTER(LEN=4)           :: hilf
!===================================================================================================================================

IF(.NOT.WritePartitionInfo) RETURN

nVars=12
IF(DoRefMapping) nVars=16
Allocate(ProcInfo(nVars))

!output partitioning info
ProcInfo( 1)=nElems
ProcInfo( 2)=nSides
ProcInfo( 3)=nTotalElems-nElems ! number of halo elems
ProcInfo( 4)=nTotalSides-nSides ! number of halo sides
ProcInfo( 5)=nPlanar
ProcInfo( 6)=nBilinear
ProcInfo( 7)=nCurved
ProcInfo( 8)=nPlanarHalo
ProcInfo( 9)=nBilinearHalo
ProcInfo(10)=nCurvedHalo
ProcInfo(11)=nLinearElems
ProcInfo(12)=nCurvedElems
IF(DoRefmapping)THEN
  ProcInfo(13)=nLinearElemsHalo
  ProcInfo(14)=nCurvedElemsHalo
  ProcInfo(15)=nBCElems
  ProcInfo(16)=nBCElemsHalo
END IF

IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:PartMPI%nProcs-1))
  ALLOCATE(ProcInfo_glob(nVars,0:PartMPI%nProcs-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot 
CALL MPI_GATHER(PartMPI%nMPINeighbors,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,PartMPI%COMM,iError)
CALL MPI_GATHER(ProcInfo,nVars,MPI_INTEGER,ProcInfo_glob,nVars,MPI_INTEGER,0,PartMPI%COMM,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(nNBmax,0:PartMPI%nProcs-1))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError) 
ALLOCATE(NBinfo(nNbmax))
NBinfo(1:PartMPI%nMPINeighbors)=PartMPI%MPINeighbor(1:PartMPi%nMPINeighbors)
CALL MPI_GATHER(NBinfo,nNBmax,MPI_INTEGER,NBinfo_glob,nNBmax,MPI_INTEGER,0,PartMPI%COMM,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  IF(DoLoadBalance)THEN
    WRITE( hilf,'(I4.4)') nLoadBalance
    filename='particlepartitionInfo-'//TRIM(hilf)//'.out'
  ELSE
    filename='particlepartitionInfo.out'
  END IF
  OPEN(UNIT=ioUnit,FILE=filename,STATUS='REPLACE')
  WRITE(ioUnit,*)'Particle Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',PartMPI%nProcs
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  IF(DoRefMapping)THEN
    WRITE(ioUnit,'(18(A20))')'Rank','nElems','nSides','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved' &
                                   ,'nPlanarHalo','nBilinearHalo','nCurvedHalo','nLinearElems','nCurvedElems'   &
                                   ,'nLinearElemsHalo','nCurvedElemsHalo','nBCElems','nBCElemsHalo','nNBProcs'
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========================================================================================='
  ELSE
    WRITE(ioUnit,'(14(A20))')'Rank','nElems','nSides','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved' &
                                   ,'nPlanarHalo','nBilinearHalo','nCurvedHalo','nLinearElems','nCurvedElems','nNBProcs'
    WRITE(ioUnit,'(A90,A90,A90,A10)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========='
  END IF

  !statistics
  ALLOCATE(tmparray(nVars+1,0:3),tmpreal(nVars+1,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=0      !max
  tmparray(:,3)=HUGE(1)   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray( 1,0)=Procinfo_glob( 1,i) ! nElems
    tmparray( 2,0)=Procinfo_glob( 2,i) ! nSides
    tmparray( 3,0)=Procinfo_glob( 3,i) ! nHaloElems
    tmparray( 4,0)=Procinfo_glob( 4,i) ! nHaloSides
    tmparray( 5,0)=Procinfo_glob( 5,i) ! nPlanar
    tmparray( 6,0)=Procinfo_glob( 6,i) ! nBilinear
    tmparray( 7,0)=Procinfo_glob( 7,i) ! nCurved
    tmparray( 8,0)=Procinfo_glob( 8,i) ! nPlanarHalo
    tmparray( 9,0)=Procinfo_glob( 9,i) ! nBilinearHalo
    tmparray(10,0)=Procinfo_glob(10,i) ! nCurvedHalo
    tmparray(11,0)=Procinfo_glob(11,i) ! nLinearElems
    tmparray(12,0)=Procinfo_glob(12,i) ! nCurvedElems
    IF(DoRefMapping)THEN
      tmparray(13,0)=Procinfo_glob(13,i) ! nLinearElemsHalo
      tmparray(14,0)=Procinfo_glob(14,i) ! nCurvedElemsHalo
      tmparray(15,0)=Procinfo_glob(15,i) ! nBCElems
      tmparray(16,0)=Procinfo_glob(16,i) ! nBCElemsHalo
    END IF
    tmparray(nVars+1,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,nVars+1
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO !j
  END DO ! i
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(PartMPI%nProcs) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,PartMPI%nProcs-1
    !actual proc
    tmparray( 1,0)=Procinfo_glob( 1,i) ! nElems
    tmparray( 2,0)=Procinfo_glob( 2,i) ! nSides
    tmparray( 3,0)=Procinfo_glob( 3,i) ! nHaloElems
    tmparray( 4,0)=Procinfo_glob( 4,i) ! nHaloSides
    tmparray( 5,0)=Procinfo_glob( 5,i) ! nPlanar
    tmparray( 6,0)=Procinfo_glob( 6,i) ! nBilinear
    tmparray( 7,0)=Procinfo_glob( 7,i) ! nCurved
    tmparray( 8,0)=Procinfo_glob( 8,i) ! nPlanarHalo
    tmparray( 9,0)=Procinfo_glob( 9,i) ! nBilinearHalo
    tmparray(10,0)=Procinfo_glob(10,i) ! nCurvedHalo
    tmparray(11,0)=Procinfo_glob(11,i) ! nLinearElems
    tmparray(12,0)=Procinfo_glob(12,i) ! nCurvedElems
    IF(DoRefMapping)THEN
      tmparray(13,0)=Procinfo_glob(13,i) ! nLinearElemsHalo
      tmparray(14,0)=Procinfo_glob(14,i) ! nCurvedElemsHalo
      tmparray(15,0)=Procinfo_glob(15,i) ! nBCElems
      tmparray(16,0)=Procinfo_glob(16,i) ! nBCElemsHalo
    END IF
    tmparray(nVars+1,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,nVars+1
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2 
    END DO ! j
  END DO ! i
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(PartMPI%nProcs))
  ! output
  IF(DoRefMapping)THEN
    WRITE(ioUnit,'(A15,17(8X,F12.4))')'   MEAN        ',tmpreal(:,1)
    WRITE(ioUnit,'(A90,A90,A90,A90,A20)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '--------------------'
    WRITE(ioUnit,'(A15,17(8X,F12.4))')'   RMS         ',tmpreal(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------'
    WRITE(ioUnit,'(A15,17(8X,I12))')'   MIN         ',tmparray(:,3)
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------'
    WRITE(ioUnit,'(A15,17(8X,I12))')'   MAX         ',tmparray(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========================================================================================='
    DO i=0,PartMPI%nProcs-1
      WRITE(ioUnit,'(17(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
      WRITE(ioUnit,'(A90,A90,A90,A90)')&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------'
    END DO! i
  ELSE
    WRITE(ioUnit,'(A15,13(8X,F12.4))')'   MEAN        ',tmpreal(:,1)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A15,13(8X,F12.4))')'   RMS         ',tmpreal(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A15,13(8X,I12))')'   MIN         ',tmparray(:,3)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A15,13(8X,I12))')'   MAX         ',tmparray(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A90,A90,A90,A10)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========='
    DO i=0,PartMPI%nProcs-1
      WRITE(ioUnit,'(13(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
      WRITE(ioUnit,'(A90,A90,A90,A10)')&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '----------'
    END DO! i
  END IF
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit) 
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)

END SUBROUTINE WriteParticlePartitionInformation


SUBROUTINE WriteParticleMappingPartitionInformation(nPlanar,nBilinear,nCurved,nTotalBCElems)
!===================================================================================================================================
! write the particle partition information to file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,            ONLY:nSides,nElems
USE MOD_Particle_MPI_Vars,    ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,   ONLY:nTotalSides,nTotalElems
USE MOD_LoadBalance_Vars,     ONLY:DoLoadBalance,nLoadBalance, writePartitionInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: nPlanar,nBilinear,nCurved,nTotalBCElems
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=10)          :: formatstr
INTEGER,ALLOCATABLE        :: nNBProcs_glob(:), ProcInfo_glob(:,:),NBInfo_glob(:,:), NBInfo(:), tmparray(:,:)
REAL,ALLOCATABLE           :: tmpreal(:,:)
INTEGER                    :: ProcInfo(8),nNBmax,i,j,ioUnit
CHARACTER(LEN=64)          :: filename
CHARACTER(LEN=4)           :: hilf
!===================================================================================================================================

IF(.NOT.WritePartitionInfo) RETURN

!output partitioning info
ProcInfo(1)=nElems
ProcInfo(2)=nSides
ProcInfo(3)=nTotalBCElems
ProcInfo(4)=nTotalElems-nElems
ProcInfo(5)=nTotalSides-nSides
ProcInfo(6)=nPlanar
ProcInfo(7)=nBilinear
ProcInfo(8)=nCurved
IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:PartMPI%nProcs-1))
  ALLOCATE(ProcInfo_glob(8,0:PartMPI%nProcs-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot 
CALL MPI_GATHER(PartMPI%nMPINeighbors,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,PartMPI%COMM,iError)
CALL MPI_GATHER(ProcInfo,8,MPI_INTEGER,ProcInfo_glob,8,MPI_INTEGER,0,PartMPI%COMM,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(nNBmax,0:PartMPI%nProcs-1))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError) 
ALLOCATE(NBinfo(nNbmax))
NBinfo(1:PartMPI%nMPINeighbors)=PartMPI%MPINeighbor(1:PartMPi%nMPINeighbors)
CALL MPI_GATHER(NBinfo,nNBmax,MPI_INTEGER,NBinfo_glob,nNBmax,MPI_INTEGER,0,PartMPI%COMM,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  IF(DoLoadBalance)THEN
    WRITE( hilf,'(I4.4)') nLoadBalance
    filename='particlepartitionInfo-'//TRIM(hilf)//'.out'
  ELSE
    filename='particlepartitionInfo.out'
  END IF
  OPEN(UNIT=ioUnit,FILE=filename,STATUS='REPLACE')
  WRITE(ioUnit,*)'Particle Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',PartMPI%nProcs
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  WRITE(ioUnit,'(10(A15))')'Rank','nElems','nSides','nBCElems','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved','nNBProcs'
  WRITE(ioUnit,'(A90,A60)')&
      '==========================================================================================',&
      '============================================================'
  !statistics
  ALLOCATE(tmparray(9,0:3),tmpreal(9,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=0      !max
  tmparray(:,3)=HUGE(1)   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i) ! nElems
    tmparray(2,0)=Procinfo_glob(2,i) ! nSides
    tmparray(3,0)=Procinfo_glob(3,i) ! nBCElems
    tmparray(4,0)=Procinfo_glob(4,i) ! nHaloElems
    tmparray(5,0)=Procinfo_glob(5,i) ! nHaloSides
    tmparray(6,0)=Procinfo_glob(6,i) ! nPlanarSides
    tmparray(7,0)=Procinfo_glob(7,i) ! nBilinearSides
    tmparray(8,0)=Procinfo_glob(8,i) ! nCurvedSides
    tmparray(9,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,9
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO !j
  END DO ! i
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(PartMPI%nProcs) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,PartMPI%nProcs-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i)
    tmparray(2,0)=Procinfo_glob(2,i)
    tmparray(3,0)=Procinfo_glob(3,i)
    tmparray(4,0)=Procinfo_glob(4,i)
    tmparray(5,0)=Procinfo_glob(5,i)
    tmparray(6,0)=Procinfo_glob(6,i)
    tmparray(7,0)=Procinfo_glob(7,i)
    tmparray(8,0)=Procinfo_glob(8,i)
    tmparray(9,0)=nNBProcs_glob(i)
    DO j=1,9
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2 
    END DO ! j
  END DO ! i
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(PartMPI%nProcs))
  WRITE(ioUnit,'(A15,9(5X,F10.2))')'   MEAN        ',tmpreal(:,1)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A15,9(5X,F10.2))')'   RMS         ',tmpreal(:,2)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A15,9(5X,I10))')'   MIN         ',tmparray(:,3)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A15,9(5X,I10))')'   MAX         ',tmparray(:,2)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A90,A60)')&
      '==========================================================================================',&
      '============================================================'
  DO i=0,PartMPI%nProcs-1
    WRITE(ioUnit,'(10(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
    WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  END DO! i
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit) 
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)

END SUBROUTINE WriteParticleMappingPartitionInformation
#endif /*MPI*/

END MODULE MOD_Particle_MPI_Halo
