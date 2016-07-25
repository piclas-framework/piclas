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

INTERFACE ExchangeMappedHaloGeometry
  MODULE PROCEDURE ExchangeMappedHaloGeometry
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


PUBLIC :: IdentifyHaloMPINeighborhood,ExchangeHaloGeometry,ExchangeMappedHaloGeometry
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
USE MOD_Mesh_Vars,                  ONLY:NGeo,ElemToSide,nElems,nInnerSides,nBCSides,nSides
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
INTEGER                     :: iBGM,jBGM,kBGM,iNBProc,iElem,ElemID,ilocSide,SideID
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
              IF(SideID.GT.0)THEN
                IF((SideID.LE.nBCSides).OR.(SideID.GT.(nBCSides+nInnerSides)))THEN
                  !IF(SideID.GT.(nInnerSides+nBCSides).AND.(SideIndex(SideID).EQ.0))THEN
                  ! because of implicit, but here I send for checking, other process sends the required halo region
                  IF(SideIndex(SideID).EQ.0)THEN
                    SendMsg%nMPISides=SendMsg%nMPISides+1
                    SideIndex(SideID)=SendMsg%nMPISides
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
DO iElem=1,nElems
  DO iLocSide=1,6
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    !IF(SideID.GT.(nInnerSides+nBCSides)) THEN
    ! only required sides
    IF(SideIndex(SideID).NE.0)THEN
      !IF(DoRefMapping)THEN
      !  SELECT CASE(iLocSide)
      !  CASE(XI_MINUS)
      !    SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,iElem)
      !    SideIndex(SideID)=0
      !  CASE(XI_PLUS)
      !    SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,iElem)
      !    SideIndex(SideID)=0
      !  CASE(ETA_MINUS)
      !    SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,iElem)
      !    SideIndex(SideID)=0
      !  CASE(ETA_PLUS)
      !    SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,iElem)
      !    SideIndex(SideID)=0
      !  CASE(ZETA_MINUS)
      !    SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,iElem)
      !    SideIndex(SideID)=0
      !  CASE(ZETA_PLUS)
      !    SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,iElem)
      !    SideIndex(SideID)=0
      !  END SELECT
      !ELSE ! no ref mapping
        SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      !END IF ! DoRefMapping
    END IF ! SideIndex NE.0
    !END IF ! Side = MPI Side
  END DO ! Side
END DO ! iElem

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
USE MOD_Particle_Mesh_Vars,         ONLY:GEO
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,         ONLY:ElemBaryNGeo,ElemRadiusNGeo
!USE MOD_Particle_Tracking_Vars,     ONLY:DoRefMapping
USE MOD_Mesh_Vars,                  ONLY:NGeo,ElemToSide,nElems,nInnerSides,nBCSides,nSides
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
INTEGER                     :: iBGM,jBGM,kBGM,iNBProc,iElem,ElemID,ilocSide,SideID
INTEGER                     :: DataSize
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
              IF(SideID.GT.0)THEN
                IF((SideID.LE.nBCSides).OR.(SideID.GT.(nBCSides+nInnerSides)))THEN
                  !IF(SideID.GT.(nInnerSides+nBCSides).AND.(SideIndex(SideID).EQ.0))THEN
                  ! because of implicit, but here I send for checking, other process sends the required halo region
                  IF(ElemIndex(SideID).EQ.0)THEN
                    SendMsg%nMPIElems=SendMsg%nMPIElems+1
                    ElemIndex(SideID)=SendMsg%nMPIElems
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
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)   :: ElemIndex(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem, ElemID,iBGMElem,iCase,NbOfElems,iNode
INTEGER                  :: iBGMmin,iBGMmax,jBGMmin,jBGMmax,kBGMmin,kBGMmax,iPBGM,jPBGM,kPBGM
INTEGER                  :: iBGM,jBGM,kBGM
REAL                     :: Vec1(1:3),Vec2(1:3),Vec3(1:3), NodeX(1:3), Radius,Distance,Vec0(1:3)
REAL                     :: xmin, xmax, ymin,ymax,zmin,zmax
!===================================================================================================================================

! For each (NGeo+1)^2 BezierControlPoint of each side, the FIBGM cell(s) in which the side 
! resides is identified and the surrouding nPaddingCells for each neighboring element
! are searched
!--- for idiots: get BezierControlPoints of myProc that are within eps distance to MPI-bound 
!                of iProc
ElemIndex=0
NBOfElems=0

DO iElem=1,nExternalElems
  ! check only the min-max values
  DO iNode=1,8
    ! get elem extension based on barycenter and radius
    NodeX(1:3) = ElemBaryAndRadius(1:3,iElem) 
    Radius     = ElemBaryAndRadius( 4 ,iElem)

    xmin = ElemBaryNGeo(1,iElem) -ElemRadiusNGeo(iElem)
    ymin = ElemBaryNGeo(2,iElem) -ElemRadiusNGeo(iElem)
    zmin = ElemBaryNGeo(3,iElem) -ElemRadiusNGeo(iElem)
    xmax = ElemBaryNGeo(1,iElem) +ElemRadiusNGeo(iElem)
    ymax = ElemBaryNGeo(2,iElem) +ElemRadiusNGeo(iElem)
    zmax = ElemBaryNGeo(3,iElem) +ElemRadiusNGeo(iElem)

    ! BGM mesh cells
    iBGMmin = CEILING((xMin-GEO%xminglob)/GEO%FIBGMdeltas(1))-FIBGMCellPadding(1)
    iBGMmin = MIN(GEO%FIBGMimax,iBGMmin)                             
    iBGMmax = CEILING((xMax-GEO%xminglob)/GEO%FIBGMdeltas(1))+FIBGMCellPadding(1)
    iBGMmax = MIN(GEO%FIBGMimax,iBGMmax)                             

    jBGMmin = CEILING((yMin-GEO%yminglob)/GEO%FIBGMdeltas(2))-FIBGMCellPadding(2)
    jBGMmin = MIN(GEO%FIBGMjmax,jBGMmin)                             
    jBGMmax = CEILING((yMax-GEO%yminglob)/GEO%FIBGMdeltas(2))+FIBGMCellPadding(2)
    jBGMmax = MIN(GEO%FIBGMjmax,jBGMmax)                             

    kBGMmin = CEILING((zMin-GEO%zminglob)/GEO%FIBGMdeltas(3))-FIBGMCellPadding(3)
    kBGMmin = MIN(GEO%FIBGMkmax,kBGMmin)                             
    kBGMmax = CEILING((zMax-GEO%zminglob)/GEO%FIBGMdeltas(3))+FIBGMCellPadding(3)
    kBGMmax = MIN(GEO%FIBGMkmax,kBGMmax)                             

    ! loop over BGM cells    
    DO iPBGM = iBGMmin,iBGMmax
      DO jPBGM = jBGMmin,jBGMmax
        DO kPBGM = kBGMmin,kBGMmax
          DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
            ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
            IF((ElemID.LE.1.).OR.(ElemID.GT.PP_nElems))CYCLE
            IF(ElemIndex(ElemID).EQ.0)THEN
              Vec0=ElemBaryNGeo(1:3,ElemID)-ElemBaryAndRadius(1:3,iElem)
              Distance=SQRT(DOT_PRODUCT(Vec0,Vec0)) &
                      +Radius+ElemRadiusNGeo(ElemID)
              IF(Distance.LE.halo_eps)THEN
                NbOfElems=NbOfElems+1
                ElemIndex(ElemID)=NbofElems
              END IF ! in range
            END IF ! ElemIndex(ElemID).EQ.0
          END DO ! iBGMElem
        END DO ! kPBGM
      END DO ! jPBGM
    END DO ! i PBGM
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
  DO iElem=1,nExternalElems
    DO iCase = 1, NbrOfCases
      IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
          (casematrix(iCase,2).EQ.0) .AND. &
          (casematrix(iCase,3).EQ.0)) CYCLE
      DO iNode=1,8
        NodeX =ElemBaryAndRadius(1:3,iElem)
        Radius=ElemBaryAndRadius( 4 ,iElem)
        SELECT CASE(iNode)
        CASE(1)
          NodeX(1:3)=NodeX(1:3) + (/-Radius,-Radius,-Radius/)
        CASE(2)
          NodeX(1:3)=NodeX(1:3) + (/Radius,-Radius,-Radius/)
        CASE(3)
          NodeX(1:3)=NodeX(1:3) + (/-Radius,Radius,-Radius/)
        CASE(4)
          NodeX(1:3)=NodeX(1:3) + (/Radius,Radius,-Radius/)
        CASE(5)
          NodeX(1:3)=NodeX(1:3) + (/-Radius,-Radius,Radius/)
        CASE(6)
          NodeX(1:3)=NodeX(1:3) + (/Radius,-Radius,Radius/)
        CASE(7)
          NodeX(1:3)=NodeX(1:3) + (/-Radius,Radius,Radius/)
        CASE(8)
          NodeX(1:3)=NodeX(1:3) + (/Radius,Radius,Radius/)
        END SELECT
        NodeX(:) = NodeX(:) + &
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
                IF((ElemID.LE.1.).OR.(ElemID.GT.PP_nElems))CYCLE
                IF(ElemIndex(ElemID).EQ.0)THEN
                  Vec0=ElemBaryNGeo(1:3,ElemID)-ElemBaryAndRadius(1:3,iElem)
                  Distance=SQRT(DOT_PRODUCT(Vec0,Vec0)) &
                          +Radius+ElemRadiusNGeo(ElemID)
                  IF(Distance.LE.halo_eps)THEN
                    NbOfElems=NbOfElems+1
                    ElemIndex(ElemID)=NbofElems
                  END IF ! in range
                END IF ! ElemIndex(ElemID).EQ.0
              END DO ! iBGMElem
            END DO ! kPBGM
          END DO ! jPBGM
        END DO ! iPBGM
      END DO ! Node=1,8
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
USE MOD_Mesh_Vars,                 ONLY:NGeo,ElemToSide,nSides,SideToElem!,XCL_NGeo
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
INTEGER                  :: iSide, NbOfSides,p,q,ElemID,ilocSide,SideID,r,s,iBGMElem,iCase,NbOfElems
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
              DO ilocSide=1,6
                SideID=ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
                !IF(DoRefMapping)THEN
                !  IF(SideIndex(SideID).EQ.0)THEN
                !    leave=.FALSE.
                !    SELECT CASE(ilocSide)
                !    CASE(XI_MINUS)
                !      xNodes=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,ElemID)
                !    CASE(XI_PLUS)
                !      xNodes=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,ElemID)
                !    CASE(ETA_MINUS)
                !      xNodes=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,ElemID)
                !    CASE(ETA_PLUS)
                !      xNodes=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,ElemID)
                !    CASE(ZETA_MINUS)
                !      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,ElemID)
                !    CASE(ZETA_PLUS)
                !      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,ElemID)
                !    END SELECT
                !    leave=.FALSE.
                !    DO s=0,NGeo
                !      DO r=0,NGeo
                !        IF(SQRT(DOT_Product(xNodes(:,r,s)-NodeX &
                !                           ,xNodes(:,r,s)-NodeX )).LE.halo_eps)THEN
                !          NbOfSides=NbOfSides+1
                !          SideIndex(SideID)=NbOfSides
                !          IF(ElemIndex(ElemID).EQ.0)THEN
                !            NbOfElems=NbOfElems+1
                !            ElemIndex(ElemID)=NbofElems
                !          END IF
                !          leave=.TRUE.
                !          EXIT
                !        END IF
                !      END DO ! r
                !      IF(leave) EXIT
                !    END DO ! s
                !  END IF ! SideIndex(SideID).EQ.0
                !ELSE ! norefmapping
                  ! caution, not save if corect
                  IF(SideIndex(SideID).EQ.0)THEN
                    leave=.FALSE.
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
                !END IF
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
                  DO ilocSide=1,6
                    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
                    ! caution, not save if corect
                    !IF(DoRefMapping)THEN
                    !  IF(SideIndex(SideID).EQ.0)THEN
                    !    leave=.FALSE.
                    !    SELECT CASE(ilocSide)
                    !    CASE(XI_MINUS)
                    !      xNodes=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,ElemID)
                    !    CASE(XI_PLUS)
                    !      xNodes=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,ElemID)
                    !    CASE(ETA_MINUS)
                    !      xNodes=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,ElemID)
                    !    CASE(ETA_PLUS)
                    !      xNodes=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,ElemID)
                    !    CASE(ZETA_MINUS)
                    !      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,ElemID)
                    !    CASE(ZETA_PLUS)
                    !      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,ElemID)
                    !    END SELECT
                    !    leave=.FALSE.
                    !    DO s=0,NGeo
                    !      DO r=0,NGeo
                    !        IF(SQRT(DOT_Product(xNodes(:,r,s)-NodeX &
                    !                           ,xNodes(:,r,s)-NodeX )).LE.halo_eps)THEN
                    !          NbOfSides=NbOfSides+1
                    !          SideIndex(SideID)=NbOfSides
                    !          IF(ElemIndex(ElemID).EQ.0)THEN
                    !            NbOfElems=NbOfElems+1
                    !            ElemIndex(ElemID)=NbofElems
                    !          END IF
                    !          leave=.TRUE.
                    !          EXIT
                    !        END IF
                    !      END DO ! r
                    !      IF(leave) EXIT
                    !    END DO ! s
                    !  END IF ! SideIndex(SideID).EQ.0
                    !ELSE ! norefmapping
                      ! caution, not save if corect
                      IF(SideIndex(SideID).EQ.0)THEN
                        leave=.FALSE.
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
                    !END IF ! DoRefMapping
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
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalSides,nTotalElems,SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElem,GEO,nTotalBCSides,ElemBaryNGeo
USE MOD_Particle_Surfaces_Vars, ONLY:ElemSlabNormals,ElemSlabIntervals  
USE MOD_Mesh_Vars,              ONLY:XCL_NGeo,dXCL_NGeo
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicDisplacement
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
  REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ElemSlabNormals    
  REAL,ALLOCATABLE,DIMENSION(:,:)     :: ElemSlabIntervals   
  REAL,ALLOCATABLE,DIMENSION(:,:)     :: ElemBaryNGeo   
  INTEGER,ALLOCATABLE       :: ElemToSide(:,:,:) 
  INTEGER,ALLOCATABLE       :: SideToElem(:,:)
  INTEGER,ALLOCATABLE       :: SideBCType(:)
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
INTEGER                     :: datasize,datasize2,datasize3
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
!  DO ilocSide=1,6
!    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
!    ! CAUTION DEBUG
!    IF(SideList(SideID).NE.0)THEN
!      IF(.NOT.isElem(iElem)) THEN
!        SendMsg%nElems=SendMsg%nElems+1
!        ElemIndex(iElem) = SendMsg%nElems
!        isElem(iElem)=.TRUE.
!      END IF ! NOT isElem
!    END IF
!  END DO ! ilocSide
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

!IPWRITE(*,*) 'sideindex,nsides',SideIndex(:),SendMsg%nSides

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

  ! ElemSlabNormals for exchange
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%ElemSlabNormals(1:3,0:3,1:SendMsg%nElems),STAT=ALLOCSTAT) 
    IF (ALLOCSTAT.NE.0) CALL abort(&
     __STAMP__&
      ,'Could not allocate SendMsg%ElemSlabNormals',SendMsg%nElems)
    SendMsg%ElemSlabNormals(:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ElemSlabNormals(1:3,0:3,1:RecvMsg%nElems),STAT=ALLOCSTAT)  
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemSlabNormals',RecvMsg%nElems)
    RecvMsg%ElemSlabNormals(:,:,:)=0
  END IF
  
  ! ElemSlabIntervals for exchange
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%ElemSlabIntervals(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT) 
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%ElemSlabIntervals',SendMsg%nElems)
    SendMsg%DXCL_NGeo(:,:,:,:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ElemSlabIntervals(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT) 
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemSlabIntervals',RecvMsg%nElems)
    RecvMsg%DXCL_NGeo(:,:,:,:,:,:)=0
  END IF
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
      SendMsg%ElemSlabNormals(:,:,ElemIndex(iElem))=ElemSlabNormals(:,:,iElem)
      SendMsg%ElemSlabIntervals(:,ElemIndex(iElem))=ElemSlabIntervals(:,iElem)
    END IF
    SendMsg%ElemBaryNGeo(:,ElemIndex(iElem))=ElemBaryNGeo(:,iElem)
    DO iLocSide = 1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      !IF(isSide(SideID))THEN
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
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemSlabNormals,SendMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemSlabIntervals,SendMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemBaryNGeo,SendMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,IERROR)

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
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemSlabNormals,RecvMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemSlabIntervals,RecvMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemBaryNGeo,RecvMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,MPISTATUS,IERROR)
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
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemSlabNormals,RecvMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemSlabIntervals,RecvMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemBaryNGeo,RecvMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,MPISTATUS,IERROR)

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
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemSlabNormals,SendMsg%nElems*12,MPI_DOUBLE_PRECISION,iProc,1116,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemSlabIntervals,SendMsg%nElems*6,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemBaryNGeo,SendMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,IERROR)

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
    ! periodic bcs are msising??, only periodic sides can be douoble sides
    ! now, there are now double sides possible
    !DO iSide=nBCSides+nInnerSides+1,nTotalSides
    !  DO iHaloSide=1,RecvMsg%nSides
    !    nDoubleBezier=0
    !    IF(  ALMOSTEQUAL(BezierControlPoints3D(1,0,0,iSide),RecvMsg%BezierControlPoints3D(1,0,0,iHaloSide))   &
    !    .AND.ALMOSTEQUAL(BezierControlPoints3D(2,0,0,iSide),RecvMsg%BezierControlPoints3D(2,0,0,iHaloSide))   &
    !    .AND.ALMOSTEQUAL(BezierControlPoints3D(3,0,0,iSide),RecvMsg%BezierControlPoints3D(3,0,0,iHaloSide)) ) &
    !      nDoubleBezier=nDoubleBezier+1
    !    IF(  ALMOSTEQUAL(BezierControlPoints3D(1,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(1,NGeo,NGeo,iHaloSide))   &
    !    .AND.ALMOSTEQUAL(BezierControlPoints3D(2,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(2,NGeo,NGeo,iHaloSide))   &
    !    .AND.ALMOSTEQUAL(BezierControlPoints3D(3,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(3,NGeo,NGeo,iHaloSide)) ) &
    !      nDoubleBezier=nDoubleBezier+1
    !    IF(nDoubleBezier.EQ.2) THEN
    !      nDoubleSides=nDoubleSides+1
    !      isSide(iHaloSide)=.FALSE.
    !    END IF
    !  END DO ! iHaloSide
    !END DO ! iSide
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
    !print*,'MyRank after resize', PartMPI%MyRank
  
    ! loop over all new elements
    !DO iElem=tmpnElems+1,nTotalElems
    DO iElem=1,RecvMsg%nElems
      !print*,'iElem',iElem
      ! first, new SideID=entry of RecvMsg+tmpnSides
      newElemID=tmpnElems+iElem
      DO ilocSide=1,6
        haloSideID=RecvMsg%ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
        IF(haloSideID.LE.0) CYCLE ! all non BC faces have not to be located
        isDoubleSide=.FALSE.
        IF(isSide(haloSideID)) THEN
          IF(HaloInc(haloSideID).EQ.0) IPWRITE(UNIT_stdOut,*) ' Warning: wrong halo inc'
          newSideID  =tmpnSides+haloinc(haloSideID)
          newBCSideID=tmpBCSides+haloinc(haloSideID)
          IF(newSideID.LT.tmpnSides) IPWRITE(UNIT_stdOut,*) 'Warning: wrong new sideid', newsideid
        ELSE ! find correct side id
          ! check if side is consistent with older side 
          DO iOldSide=nBCSides+1,tmpnSides
            nDoubleBezier=0
            IF(  ALMOSTEQUAL(BezierControlPoints3D(1,0,0,iOldSide),RecvMsg%BezierControlPoints3D(1,0,0,haloSideID))   &
            .AND.ALMOSTEQUAL(BezierControlPoints3D(2,0,0,iOldSide),RecvMsg%BezierControlPoints3D(2,0,0,haloSideID))   &
            .AND.ALMOSTEQUAL(BezierControlPoints3D(3,0,0,iOldSide),RecvMsg%BezierControlPoints3D(3,0,0,haloSideID)) ) & 
              nDoubleBezier=nDoubleBezier+1
            IF(  ALMOSTEQUAL(BezierControlPoints3D(1,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(1,NGeo,NGeo,haloSideID))   &
            .AND.ALMOSTEQUAL(BezierControlPoints3D(2,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(2,NGeo,NGeo,haloSideID))   &
            .AND.ALMOSTEQUAL(BezierControlPoints3D(3,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(3,NGeo,NGeo,haloSideID)) ) &
              nDoubleBezier=nDoubleBezier+1
            IF(nDoubleBezier.EQ.2) THEN
              isDoubleSide=.TRUE.
              EXIT
            END IF
          END DO ! iOldSide
          newSideID=iOldSide
        END IF
        IF(isDoubleSide)THEN ! equals IF(.NOT.isSide(haloSideID))
          ! here something fancy with SideID
          ! check master and slave flip
          oldElemID=PartSideToElem(S2E_NB_ELEM_ID,newSideID)
          ! if neighbor ElemID equals -1, then iElem is unknown in old grid
          IF(oldElemID.EQ.-1)THEN
            ! get elemid
            IF(PartSideToElem(S2E_ELEM_ID,newSideID).EQ.-1) &
              CALL abort(&
__STAMP__&
              ,'Critical error in domain reconstrution.')
              PartSideToElem(S2E_NB_ELEM_ID,newSideID)     = newElemID
              PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = ilocSide
              ! nothing to do, is already filled
              !PartSideToElem(S2E_ELEM_ID       ,newSideID) = 
              !PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = 
              !PartSideToElem(S2E_FLIP          ,newSideID) = 
              ! NeighboreElemID
              PartElemToElem(E2E_NB_ELEM_ID,PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))&
                                                                                                                  =newElemID
              PartElemToElem(E2E_NB_LOC_SIDE_ID,PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))&
                                                                                                                  =ilocSide
              PartElemToElem(E2E_NB_ELEM_ID,ilocSide,newElemID) = PartSideToElem(S2E_ELEM_ID,newSideID)
              PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,newElemID) = PartSideToElem(S2E_LOC_SIDE_ID,newSideID)
          ELSE ! SE2_NB_ELEM_ID=DEFINED
            IF(PartSideToElem(S2E_ELEM_ID,newSideID).NE.-1) &
              CALL abort(&
__STAMP__&
            , 'Critical error in domain reconstrution.')
            PartSideToElem(S2E_ELEM_ID       ,newSideID) = newElemID !root Element
            PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = iLocSide
            PartSideToElem(S2E_FLIP          ,newSideID) = 0
            ! already filled
            !PartSideToElem(S2E_NB_ELEM_ID,newSide)       = 
            !PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = 
          PartElemToElem(E2E_NB_ELEM_ID,PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))&
                                                                                                                =newElemID
          PartElemToElem(E2E_NB_LOC_SIDE_ID,PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))&
                                                                                                                =ilocSide
            PartElemToElem(E2E_NB_ELEM_ID,ilocSide,newElemID) = PartSideToElem(S2E_NB_ELEM_ID,newSideID)
            PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,newElemID) = PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)
          END IF
          isDone(haloSideID)=.TRUE.
        ELSE ! non-double side || new side
          ! cannnot build ElemToElem and PartElemToElem yet
          ! build PartSideToElem, so much as possible
          ! get correct side out of RecvMsg%SideToElem
          IF(iElem.EQ.RecvMsg%SideToElem(S2E_ELEM_ID,haloSideID))THEN
            ! master side
            PartSideToElem(S2E_ELEM_ID,newSideID)    =newElemID
            PartSideToElem(S2E_LOC_SIDE_ID,newSideID)=ilocSide
            PartSideToElem(S2E_FLIP       ,newSideID)=0 !RecvMsg%SideToElem(S2E_FLIP,haloSideID)
          ELSE IF(iElem.EQ.RecvMsg%SideToElem(S2E_NB_ELEM_ID,haloSideID))THEN
            ! slave side
            ! should never happens
            PartSideToElem(S2E_NB_ELEM_ID,newSideID) =newElemID
            PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) =ilocSide
            PartSideToElem(S2E_FLIP      ,newSideID) =RecvMsg%SideToElem(S2E_FLIP,haloSideID)
          ELSE ! should be found, because there should be halo sides without any connection
            CALL abort(&
__STAMP__&
            ,'Non-Critical error in domain reconstrution. IF NOT encountered, something is terrible wrong.')
          END IF
          !BC(1:4,newSideID)=RecvMsg%BC(1:4,haloSideID)
          BC(newSideID)=RecvMsg%BC(haloSideID)
          SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
        END IF ! isDoubleSide
        ! copy Bezier to new side id
        IF(.NOT.isDoubleSide) THEN
          BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newBCSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
          ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
          SideSlabNormals(1:3,1:3,newBCSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
          SideSlabIntervals(1:6,newBCSideID) =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
          BoundingBoxIsEmpty(newBCSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
        ELSE
          IF(RecvMsg%ElemToSide(2,ilocSide,iElem).EQ.0)THEN
            SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
            SideSlabIntervals(1:6,newSideID) =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
            BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
          END IF
        END IF
        ! build entry to PartElemToSide
        BC(newSideID)=RecvMsg%BC(haloSideID)
        SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
        PartBCSideList(newSideID)=newBCSideID !tmpBCSides+haloinc(haloSideID)
        PartElemToSide(E2S_SIDE_ID,iLocSide,newElemId)=newSideID
        PartElemToSide(E2S_FLIP,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
      END DO ! ilocSide
      ! set native elemID
      PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
      PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
      XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
      dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
      ElemSlabNormals(1:3,0:3,newElemID) = RecvMsg%ElemSlabNormals(1:3,0:3,iElem)
      ElemSlabIntervals(1:6,newElemID) = RecvMsg%ElemSlabIntervals(1:6,iElem)
      ElemBaryNGeo(1:3,newElemID) = RecvMsg%ElemBaryNGeo(1:3,iElem)
    END DO ! iElem
    ! missing connection info to MPI-Neighbor-ElemID
    DO iSide=nBCSides+nInnerSides+1,nSides
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      ElemID2=PartSideToElem(S2E_NB_ELEM_ID,iSide)
      IF( ElemID .NE. -1 .AND. ElemID2 .NE. -1 ) CYCLE
      locSideID=-1
      IF( ElemID .GE. 1 .AND. ElemID .LE. PP_nElems) THEN
        HostElemID=ElemID
        locSideID = SideToElem(S2E_LOC_SIDE_ID,iSide)
        flip=0
      END IF
      IF( ElemID2 .GE. 1 .AND. ElemID2 .LE. PP_nElems) THEN
        HostElemID=ElemID2
        locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,iSide)
        flip=1
      END IF
      SELECT CASE(locSideID)
      CASE(XI_MINUS)
        xNodes=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,HostelemID)
      CASE(XI_PLUS)
        xNodes=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,HostelemID)
      CASE(ETA_MINUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,HostelemID)
      CASE(ETA_PLUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,HostelemID)
      CASE(ZETA_MINUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,HostelemID)
      CASE(ZETA_PLUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,HostelemID)
      END SELECT
      IF(GEO%nPeriodicVectors.GT.0)THEN
        iDisplace= SidePeriodicType(iSide)
        IF(iDisplace.EQ.0) CYCLE
        IF(flip.EQ.0)THEN
          DO q=0,NGeo
            DO p=0,NGeo
              xNodes(1:3,q,p)=xNodes(1:3,q,p)+SidePeriodicDisplacement(1:3,iDisplace)
            END DO ! p=0,PP_N
          END DO ! q=0,PP_N
        ELSE
          DO q=0,NGeo
            DO p=0,NGeo
              xNodes(1:3,q,p)=xNodes(1:3,q,p)-SidePeriodicDisplacement(1:3,iDisplace)
            END DO ! p=0,PP_N
          END DO ! q=0,PP_N
        END IF
      END IF
      DO iElem=PP_nElems+1,nTotalElems
        DO ilocSide=1,6
          SELECT CASE(ilocSide)
          CASE(XI_MINUS)
            xNodes2=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,iElem)
          CASE(XI_PLUS)
            xNodes2=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,iElem)
          CASE(ETA_MINUS)
            xNodes2=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,iElem)
          CASE(ETA_PLUS)
            xNodes2=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,iElem)
          CASE(ZETA_MINUS)
            xNodes2=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,iElem)
          CASE(ZETA_PLUS)
            xNodes2=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,iElem)
          END SELECT
  
          ! BUG: works only if both local coord systems are not rotated
          IF(  ALMOSTEQUAL(xNodes(1,0,0),xNodes2(1,0,0))               &
          .AND.ALMOSTEQUAL(xNodes(2,0,0),xNodes2(2,0,0))               &
          .AND.ALMOSTEQUAL(xNodes(3,0,0),xNodes2(3,0,0))               & 
          .AND.ALMOSTEQUAL(xNodes(1,NGeo,NGeo),xNodes2(1,NGeo,NGeo))   &
          .AND.ALMOSTEQUAL(xNodes(2,NGeo,NGeo),xNodes2(2,NGeo,NGeo))   &
          .AND.ALMOSTEQUAL(xNodes(3,NGeo,NGeo),xNodes2(3,NGeo,NGeo)) ) THEN
            IF(flip.EQ.0)THEN
              PartSideToElem(S2E_NB_ELEM_ID,iSide)=iElem
              EXIT
            ELSE
              PartSideToElem(S2E_ELEM_ID,iSide)=iElem
              EXIT
            END IF
          END IF
        END DO ! ilocSide=1,6
      END DO ! iElem=PP_nElems+1,nTotalElems
    END DO ! iSide=nBCSides+nInnerSides+1,nSides
    ! build rest: PartElemToElem, PartLocSideID
    DO iElem=PP_nElems+1,nTotalElems
      DO ilocSide=1,6
        flip   = PartElemToSide(E2S_FLIP,ilocSide,iElem)
        SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        ! check of sideid
        HaloSideID=SideID-tmpnSides
        IF(HaloSideID.LE.tmpnSides) CYCLE 
        IF(isDone(HaloSideID)) CYCLE
        IF(flip.EQ.0)THEN
          ! SideID of slave
          PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)=PartSideToElem(S2E_NB_LOC_SIDE_ID,SideID)
          PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)=PartSideToElem(S2E_NB_ELEM_ID,SideID)
        ELSE
          ! SideID of master
          PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)=PartSideToElem(S2E_LOC_SIDE_ID,SideID)
          PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)=PartSideToElem(S2E_ELEM_ID,SideID)
        END IF
        isDone(HaloSideID)=.TRUE.
      END DO ! ilocSide
    END DO ! Elem
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
          ,' SideSlabNormals is zero!')
      END IF
    END DO 
  END IF ! RecvMsg%nSides>0

ELSE
 IF (RecvMsg%nSides.GT.0) THEN
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
   DO iSide=nBCSides+nInnerSides+1,nTotalSides
     DO iHaloSide=1,RecvMsg%nSides
       nDoubleBezier=0
       IF(  ALMOSTEQUAL(BezierControlPoints3D(1,0,0,iSide),RecvMsg%BezierControlPoints3D(1,0,0,iHaloSide))   &
       .AND.ALMOSTEQUAL(BezierControlPoints3D(2,0,0,iSide),RecvMsg%BezierControlPoints3D(2,0,0,iHaloSide))   &
       .AND.ALMOSTEQUAL(BezierControlPoints3D(3,0,0,iSide),RecvMsg%BezierControlPoints3D(3,0,0,iHaloSide)) ) &
         nDoubleBezier=nDoubleBezier+1
       IF(  ALMOSTEQUAL(BezierControlPoints3D(1,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(1,NGeo,NGeo,iHaloSide))   &
       .AND.ALMOSTEQUAL(BezierControlPoints3D(2,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(2,NGeo,NGeo,iHaloSide))   &
       .AND.ALMOSTEQUAL(BezierControlPoints3D(3,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(3,NGeo,NGeo,iHaloSide)) ) &
         nDoubleBezier=nDoubleBezier+1
       IF(nDoubleBezier.EQ.2) THEN
         nDoubleSides=nDoubleSides+1
         isSide(iHaloSide)=.FALSE.
       END IF
     END DO ! iHaloSide
   END DO ! iSide
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
   !print*,'MyRank,nSides,nnewSides,nDoubleSides', PartMPI%MyRank,nSides,RecvMsg%nSides,nDoubleSides
   tmpnSides =nTotalSides
   tmpnElems=nTotalElems
   nTotalSides=nTotalSides+RecvMsg%nSides-nDoubleSides
   nTotalElems=nTotalElems+RecvMsg%nElems
   CALL ResizeParticleMeshData(tmpnSides,tmpnElems,nTotalSides,nTotalElems)
   !print*,'MyRank after resize', PartMPI%MyRank
 
   ! loop over all new elements
   !DO iElem=tmpnElems+1,nTotalElems
   DO iElem=1,RecvMsg%nElems
     !print*,'iElem',iElem
     ! first, new SideID=entry of RecvMsg+tmpnSides
     newElemID=tmpnElems+iElem
     DO ilocSide=1,6
       haloSideID=RecvMsg%ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
       ! first, set new sideid
     !  print*,'haloSideId',haloSideID
       isDoubleSide=.FALSE.
       IF(isSide(haloSideID)) THEN
         IF(HaloInc(haloSideID).EQ.0) IPWRITE(UNIT_stdOut,*) ' Warning: wrong halo inc'
         newSideID=tmpnSides+haloinc(haloSideID)
         IF(newSideID.LT.tmpnSides) IPWRITE(UNIT_stdOut,*) 'Warning: wrong new sideid', newsideid
       ELSE ! find correct side id
         ! check if side is consistent with older side 
         DO iOldSide=nBCSides+1,tmpnSides
           nDoubleBezier=0
           IF(  ALMOSTEQUAL(BezierControlPoints3D(1,0,0,iOldSide),RecvMsg%BezierControlPoints3D(1,0,0,haloSideID))   &
           .AND.ALMOSTEQUAL(BezierControlPoints3D(2,0,0,iOldSide),RecvMsg%BezierControlPoints3D(2,0,0,haloSideID))   &
           .AND.ALMOSTEQUAL(BezierControlPoints3D(3,0,0,iOldSide),RecvMsg%BezierControlPoints3D(3,0,0,haloSideID)) ) & 
             nDoubleBezier=nDoubleBezier+1
           IF(  ALMOSTEQUAL(BezierControlPoints3D(1,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(1,NGeo,NGeo,haloSideID))   &
           .AND.ALMOSTEQUAL(BezierControlPoints3D(2,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(2,NGeo,NGeo,haloSideID))   &
           .AND.ALMOSTEQUAL(BezierControlPoints3D(3,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(3,NGeo,NGeo,haloSideID)) ) &
             nDoubleBezier=nDoubleBezier+1
           IF(nDoubleBezier.EQ.2) THEN
             isDoubleSide=.TRUE.
             EXIT
           END IF
         END DO ! iOldSide
         newSideID=iOldSide
       END IF
       IF(isDoubleSide)THEN ! equals IF(.NOT.isSide(haloSideID))
         ! here something fancy with SideID
         ! check master and slave flip
         oldElemID=PartSideToElem(S2E_NB_ELEM_ID,newSideID)
         ! if neighbor ElemID equals -1, then iElem is unknown in old grid
         IF(oldElemID.EQ.-1)THEN
           ! get elemid
           IF(PartSideToElem(S2E_ELEM_ID,newSideID).EQ.-1) &
             CALL abort(&
__STAMP__&
            ,'Critical error in domain reconstrution.')
             PartSideToElem(S2E_NB_ELEM_ID,newSideID)     = newElemID
             PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = ilocSide
             ! here? CAUTION// DEBUG
             PartSideToElem(S2E_FLIP          ,newSideID) = RecvMsg%SideToElem(S2E_FLIP,haloSideID)
             ! nothing to do, is already filled
             !PartSideToElem(S2E_ELEM_ID       ,newSideID) = 
             !PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = 
             !PartSideToElem(S2E_FLIP          ,newSideID) = 
             ! NeighboreElemID
             PartElemToElem(E2E_NB_ELEM_ID,PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))&
                                                                                                                     =newElemID
             PartElemToElem(E2E_NB_LOC_SIDE_ID,PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))&
                                                                                                                     =ilocSide
             PartElemToElem(E2E_NB_ELEM_ID,ilocSide,newElemID)    = PartSideToElem(S2E_ELEM_ID,newSideID)
             PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,newElemID) = PartSideToElem(S2E_LOC_SIDE_ID,newSideID)
           !  IF(PartSideToElem(S2E_ELEM_ID,newSideID).EQ.-1) IPWRITE(UNIT_stdOut,*) 'warning'
         ELSE ! SE2_NB_ELEM_ID=DEFINED
           IF(PartSideToElem(S2E_ELEM_ID,newSideID).NE.-1) &
             CALL abort(&
__STAMP__&
           ,'Critical error in domain reconstrution.')
           PartSideToElem(S2E_ELEM_ID       ,newSideID) = newElemID !root Element
           PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = iLocSide
           !PartSideToElem(S2E_FLIP          ,newSideID) = RecvMsg%SideToElem(S2E_FLIP
           ! already filled
           !PartSideToElem(S2E_NB_ELEM_ID,newSide)       = 
           !PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = 
           PartElemToElem(E2E_NB_ELEM_ID,PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID)) &
                                                                                                                       =newElemID
          PartElemToElem(E2E_NB_LOC_SIDE_ID,PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID)) &
                                                                                                                      =ilocSide
           PartElemToElem(E2E_NB_ELEM_ID,ilocSide,newElemID) = PartSideToElem(S2E_NB_ELEM_ID,newSideID)
           PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,newElemID) = PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)
           !IF(PartSideToElem(S2E_NB_ELEM_ID,newSideID).EQ.-1) IPWRITE(UNIT_stdOut,*)'warning'
         END IF
         isDone(haloSideID)=.TRUE.
       ELSE ! non-double side || new side
         ! cannnot build PartElemToElem and PartElemToElem yet
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
         BC(newSideID)=RecvMsg%BC(haloSideID)
         SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
       END IF ! isDoubleSide
       ! copy Bezier to new side id
       IF(.NOT.isDoubleSide) THEN
         BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
         ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
         SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
         SideSlabIntervals(1:6,newSideID) =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
         BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
       ELSE
         IF(RecvMsg%ElemToSide(2,ilocSide,iElem).EQ.0)THEN
           SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
           SideSlabIntervals(1:6,newSideID) =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
           BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
         END IF
       END IF
       ! build entry to PartElemToSide
       PartElemToSide(1,iLocSide,newElemId)=newSideID
       PartElemToSide(2,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
     END DO ! ilocSide
     ! set native elemID
     PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
     PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
     ElemBaryNGeo(1:3,newElemID) = RecvMsg%ElemBaryNGeo(1:3,iElem)
   END DO ! iElem
   ! build rest: PartElemToElem, PartLocSideID
   DO iElem=PP_nElems+1,nTotalElems
     DO ilocSide=1,6
       flip   = PartElemToSide(E2S_FLIP,ilocSide,iElem)
       SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
       ! check of sideid
       !HaloSideID=SideID-tmpnSides
       !print*,'HaloSideID',HaloSideID
       ! do not double sides
       ! debug commented out
       !IF(HaloSideID.LE.tmpnSides) CYCLE 
       !IF(isDone(HaloSideID)) CYCLE
       IF(flip.EQ.0)THEN
         ! SideID of slave
         PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)=PartSideToElem(S2E_NB_LOC_SIDE_ID,SideID)
         PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)=PartSideToElem(S2E_NB_ELEM_ID,SideID)
       ELSE
         ! SideID of master
         PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)=PartSideToElem(S2E_LOC_SIDE_ID,SideID)
         PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)=PartSideToElem(S2E_ELEM_ID,SideID)
       END IF
       !isDone(HaloSideID)=.TRUE.
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
USE MOD_Mesh_Vars,              ONLY:BC,nGeo,nElems,XCL_NGeo,DXCL_NGEO
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElem,ElemBaryNGeo
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars, ONLY:ElemSlabNormals,ElemSlabIntervals  
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
INTEGER,ALLOCATABLE                :: DummySideToElem(:,:)
INTEGER,ALLOCATABLE                :: DummySideBCType(:),DummyPartBCSideList(:)
INTEGER,ALLOCATABLE                :: DummyElemToElem(:,:,:)
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
  ! ElemSlabNormals
  ALLOCATE(DummyElemSlabNormals(1:3,0:3,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyElemSlabNormals)) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  DummyElemSlabNormals=ElemSlabNormals
  DEALLOCATE(ElemSlabNormals)
  ALLOCATE(ElemSlabNormals(1:3,0:3,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  ElemSlabNormals=0
  ElemSlabNormals(1:3,0:3,1:nOldElems) =DummyElemSlabNormals(1:3,0:3,1:nOldElems)
  DEALLOCATE(DummyElemSlabNormals)
  ! ElemSlabIntervals
  ALLOCATE(DummyElemSlabIntervals(1:6,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyElemSlabIntervals)) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  DummyElemSlabIntervals=ElemSlabIntervals
  DEALLOCATE(ElemSlabIntervals)
  ALLOCATE(ElemSlabIntervals(1:6,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate ElemIndex')
  ElemSlabIntervals=0
  ElemSlabIntervals(1:6,1:nOldElems) =DummyElemSlabIntervals(1:6,1:nOldElems)
  DEALLOCATE(DummyElemSlabIntervals)
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
PartSideToElem(:,1:nOldSides  )              =DummySideToElem(:,1:nOldSides)
DEALLOCATE(DummySideToElem)
!print*,' done side to elem',myrank
! PartElemToElem
ALLOCATE(DummyElemToElem(1:2,1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemToElem)) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
DummyElemToElem=PartElemToElem
DEALLOCATE(PartElemToElem)
ALLOCATE(PartElemToElem(1:2,1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not allocate ElemIndex')
PartElemToElem=-1
PartElemToElem(1:2,1:6,1:nOldElems)            =DummyElemToElem(1:2,1:6,1:nOldElems)
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

END SUBROUTINE ResizeParticleMeshData


!SUBROUTINE ExchangeMappedHaloGeometry(iProc,SideList,ElemList)
SUBROUTINE ExchangeMappedHaloGeometry(iProc,ElemList)
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
USE MOD_Mesh_Vars,              ONLY:nElems, nSides, nBCSides, ElemToSide, BC,nGeo,SideToElem,nInnerSides
USE MOD_Mesh_Vars,              ONLY:XCL_NGeo,dXCL_NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalSides,nTotalElems,SidePeriodicType,PartBCSideList,nTotalBCSides,GEO
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElem,SidePeriodicDisplacement
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
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
  INTEGER,ALLOCATABLE       :: ElemToSide(:,:,:) 
  INTEGER,ALLOCATABLE       :: SideToElem(:,:)
  INTEGER,ALLOCATABLE       :: SideBCType(:)
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
REAL                        :: xNodes(1:3,0:NGeo,0:NGeo)
REAL                        :: xNodes2(1:3,0:NGeo,0:NGeo)
INTEGER                     :: ALLOCSTAT,tmpBCSides,newBCSideID,locSideID
INTEGER                     :: newSideID,haloSideID,ioldSide,oldElemID,newElemID,iDisplace,p,q
LOGICAL                     :: isDoubleSide
LOGICAL,ALLOCATABLE         :: isElem(:),isSide(:),isDone(:)
INTEGER, ALLOCATABLE        :: ElemIndex(:), SideIndex(:),HaloInc(:)
INTEGER                     :: iElem, ilocSide,SideID,iSide,iIndex,iHaloSide,flip,ElemID,ElemID2, HostElemID
INTEGER                     :: nDoubleSides,nDoubleBezier,tmpnSides,tmpnElems
INTEGER                     :: datasize,datasize2,datasize3
!===================================================================================================================================


ALLOCATE(isElem(1:nElems))
IF (.NOT.ALLOCATED(isElem)) CALL abort(&
  __STAMP__&
  ,'Could not allocate isElem')
isElem(:) = .FALSE.

ALLOCATE(isSide(1:nSides))
IF (.NOT.ALLOCATED(isElem)) CALL abort(&
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
  END IF ! NOT isElem
  !END IF
!  DO ilocSide=1,6
!    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
!    ! CAUTION DEBUG
!    IF(SideList(SideID).NE.0)THEN
!      IF(.NOT.isElem(iElem)) THEN
!        SendMsg%nElems=SendMsg%nElems+1
!        ElemIndex(iElem) = SendMsg%nElems
!        isElem(iElem)=.TRUE.
!      END IF ! NOT isElem
!    END IF
!  END DO ! ilocSide
END DO ! iElem

! 2) mark all required sides and get number of send sides
DO iElem=1,nElems
  IF(isElem(iElem))THEN
    DO ilocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      ! for reference mapping, always send only bc sides, all other sides are not required
      !IF(.NOT.isSide(SideID).AND.PartBCSideList(SideID).LE.nTotalBCSides) THEN
      IF(.NOT.isSide(SideID).AND.SideID.LE.nBCSides) THEN
        ! missing: what do do with BC sides??"
        SendMsg%nSides=SendMsg%nSides+1
        SideIndex(SideID) = SendMsg%nSides
        isSide(SideID)=.TRUE.
      END IF ! not isSide
    END DO ! ilocSide
  END IF ! Element is marked to send
END DO ! iElem

!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!print*,'MyRank, sides to send, elem to send', PartMPI%MyRank,SendMsg%nSides,SendMsg%nElems
!IF(PartMPI%MPIROOT) print*,'SideIndex',SideIndex(:)
!IF(PartMPI%MPIROOT) print*,'isElem',isElem(:)
!IF(PartMPI%MPIROOT) print*,'sendnelems',SendMsg%nElems
!IF(PartMPI%MPIROOT) print*,'isSide',isSide
!IPWRITE(UNIT_stdOut,*) ' Send number of elems,sides', SendMsg%nElems, SendMsg%nSides

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
    SendMsg%XCL_NGeo(:,:,:,:,ElemIndex(iElem))=XCL_NGeo(:,:,:,:,iElem)
    SendMsg%dXCL_NGeo(:,:,:,:,:,ElemIndex(iElem))=dXCL_NGeo(:,:,:,:,:,iElem)
    DO iLocSide = 1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      IF(isSide(SideID))THEN
        SendMsg%ElemToSide(1,iLocSide,ElemIndex(iElem)) = &
                 SideIndex(ElemToSide(E2S_SIDE_ID,iLocSide,iElem))
             ! CAUTION DEBUG correct sideid????
        SendMsg%ElemToSide(2,iLocSide,ElemIndex(iElem)) = &
                ElemToSide(2,iLocSide,iElem)
      END IF
    END DO
  END IF
END DO
! SideToElem Mapping & BezierControlPoints3D
DO iSide = 1,nSides
  !IF (SideIndex(iSide).NE.0) THEN
  IF(isSide(iSide))THEN
    DO iIndex = 1,2   ! S2E_ELEM_ID, S2E_NB_ELEM_ID
      IF (SideToElem(iIndex,iSide).GT.0) THEN
         SendMsg%SideToElem(iIndex,SideIndex(iSide)) = &
                 ElemIndex(SideToElem(iIndex,iSide))     ! gets 0 if pointing to a non-Halo-cell.
                                                         ! if 0 and no BC is associated this 
                                                         ! is a Halo-bound-cell
         SendMsg%SideBCType(SideIndex(iSide)) = SidePeriodicType(iSide)
      END IF
    END DO ! S2E_LOC_SIDE_ID, S2E_NB_LOC_SIDE_ID, S2E_FLIP
   ! IF(PartMPI%MPIROOT) print*,'iSide,SideIndex',iSide,SideIndex(iSide)
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
  END IF
END DO


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

dataSize =3*(NGeo+1)*(NGeo+1)
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
      CALL MPI_SEND(SendMsg%SideSlabNormals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabIntervals,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)

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
      CALL MPI_RECV(RecvMsg%SideSlabNormals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabIntervals,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)

  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
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
      CALL MPI_RECV(RecvMsg%SideSlabNormals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabIntervals,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)


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
      CALL MPI_SEND(SendMsg%SideSlabNormals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabIntervals,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)

END IF

IF ((RecvMsg%nElems.EQ.0) .AND. (RecvMsg%nSides.GT.0))THEN
    ERRWRITE(*,*)'ERROR: nElems=0 when nSides=',RecvMsg%nSides,' and nSides=',RecvMsg%nSides,'!'
    CALL abort(&
    __STAMP__&
    ,'nElems=0 while nSides=',RecvMsg%nSides)
END IF

DEALLOCATE(isElem,isSide,ElemIndex,SideIndex)


!IF (RecvMsg%nSides.GT.0) THEN
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
  ! periodic bcs are msising??, only periodic sides can be douoble sides
  ! now, there are now double sides possible
  !DO iSide=nBCSides+nInnerSides+1,nTotalSides
  !  DO iHaloSide=1,RecvMsg%nSides
  !    nDoubleBezier=0
  !    IF(  ALMOSTEQUAL(BezierControlPoints3D(1,0,0,iSide),RecvMsg%BezierControlPoints3D(1,0,0,iHaloSide))   &
  !    .AND.ALMOSTEQUAL(BezierControlPoints3D(2,0,0,iSide),RecvMsg%BezierControlPoints3D(2,0,0,iHaloSide))   &
  !    .AND.ALMOSTEQUAL(BezierControlPoints3D(3,0,0,iSide),RecvMsg%BezierControlPoints3D(3,0,0,iHaloSide)) ) &
  !      nDoubleBezier=nDoubleBezier+1
  !    IF(  ALMOSTEQUAL(BezierControlPoints3D(1,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(1,NGeo,NGeo,iHaloSide))   &
  !    .AND.ALMOSTEQUAL(BezierControlPoints3D(2,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(2,NGeo,NGeo,iHaloSide))   &
  !    .AND.ALMOSTEQUAL(BezierControlPoints3D(3,NGeo,NGeo,iSide),RecvMsg%BezierControlPoints3D(3,NGeo,NGeo,iHaloSide)) ) &
  !      nDoubleBezier=nDoubleBezier+1
  !    IF(nDoubleBezier.EQ.2) THEN
  !      nDoubleSides=nDoubleSides+1
  !      isSide(iHaloSide)=.FALSE.
  !    END IF
  !  END DO ! iHaloSide
  !END DO ! iSide
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
  tmpnSides    =nTotalSides
  tmpnElems    =nTotalElems
  tmpBCSides   =nTotalBCSides
  nTotalSides  =nTotalSides+RecvMsg%nSides-nDoubleSides
  nTotalBCSides=nTotalBCSides+RecvMsg%nSides-nDoubleSides
  nTotalElems  =nTotalElems+RecvMsg%nElems
  CALL ResizeParticleMeshData(tmpnSides,tmpnElems,nTotalSides,nTotalElems,tmpBCSides,nTotalBCSides)
  !print*,'MyRank after resize', PartMPI%MyRank

  ! loop over all new elements
  !DO iElem=tmpnElems+1,nTotalElems
  DO iElem=1,RecvMsg%nElems
    !print*,'iElem',iElem
    ! first, new SideID=entry of RecvMsg+tmpnSides
    newElemID=tmpnElems+iElem
    DO ilocSide=1,6
      haloSideID=RecvMsg%ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      IF(haloSideID.EQ.0) CYCLE ! all non BC faces have not to be located
      isDoubleSide=.FALSE.
      IF(isSide(haloSideID)) THEN
        IF(HaloInc(haloSideID).EQ.0) IPWRITE(UNIT_stdOut,*) ' Warning: wrong halo inc'
        newSideID  =tmpnSides+haloinc(haloSideID)
        newBCSideID=tmpBCSides+haloinc(haloSideID)
        IF(newSideID.LT.tmpnSides) IPWRITE(UNIT_stdOut,*) 'Warning: wrong new sideid', newsideid
      ELSE ! find correct side id
        ! check if side is consistent with older side 
        DO iOldSide=nBCSides+1,tmpnSides
          nDoubleBezier=0
          IF(  ALMOSTEQUAL(BezierControlPoints3D(1,0,0,iOldSide),RecvMsg%BezierControlPoints3D(1,0,0,haloSideID))   &
          .AND.ALMOSTEQUAL(BezierControlPoints3D(2,0,0,iOldSide),RecvMsg%BezierControlPoints3D(2,0,0,haloSideID))   &
          .AND.ALMOSTEQUAL(BezierControlPoints3D(3,0,0,iOldSide),RecvMsg%BezierControlPoints3D(3,0,0,haloSideID)) ) & 
            nDoubleBezier=nDoubleBezier+1
          IF(  ALMOSTEQUAL(BezierControlPoints3D(1,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(1,NGeo,NGeo,haloSideID))   &
          .AND.ALMOSTEQUAL(BezierControlPoints3D(2,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(2,NGeo,NGeo,haloSideID))   &
          .AND.ALMOSTEQUAL(BezierControlPoints3D(3,NGeo,NGeo,iOldSide),RecvMsg%BezierControlPoints3D(3,NGeo,NGeo,haloSideID)) ) &
            nDoubleBezier=nDoubleBezier+1
          IF(nDoubleBezier.EQ.2) THEN
            isDoubleSide=.TRUE.
            EXIT
          END IF
        END DO ! iOldSide
        newSideID=iOldSide
      END IF
      IF(isDoubleSide)THEN ! equals IF(.NOT.isSide(haloSideID))
        ! here something fancy with SideID
        ! check master and slave flip
        oldElemID=PartSideToElem(S2E_NB_ELEM_ID,newSideID)
        ! if neighbor ElemID equals -1, then iElem is unknown in old grid
        IF(oldElemID.EQ.-1)THEN
          ! get elemid
          IF(PartSideToElem(S2E_ELEM_ID,newSideID).EQ.-1) &
            CALL abort(&
__STAMP__&
            ,'Critical error in domain reconstrution.')
            PartSideToElem(S2E_NB_ELEM_ID,newSideID)     = newElemID
            PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = ilocSide
            ! nothing to do, is already filled
            !PartSideToElem(S2E_ELEM_ID       ,newSideID) = 
            !PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = 
            !PartSideToElem(S2E_FLIP          ,newSideID) = 
            ! NeighboreElemID
            PartElemToElem(E2E_NB_ELEM_ID,PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))&
                                                                                                                =newElemID
            PartElemToElem(E2E_NB_LOC_SIDE_ID,PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))&
                                                                                                                =ilocSide
            PartElemToElem(E2E_NB_ELEM_ID,ilocSide,newElemID) = PartSideToElem(S2E_ELEM_ID,newSideID)
            PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,newElemID) = PartSideToElem(S2E_LOC_SIDE_ID,newSideID)
        ELSE ! SE2_NB_ELEM_ID=DEFINED
          IF(PartSideToElem(S2E_ELEM_ID,newSideID).NE.-1) &
            CALL abort(&
__STAMP__&
            ,'Critical error in domain reconstrution.')
          PartSideToElem(S2E_ELEM_ID       ,newSideID) = newElemID !root Element
          PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = iLocSide
          PartSideToElem(S2E_FLIP          ,newSideID) = 0
          ! already filled
          !PartSideToElem(S2E_NB_ELEM_ID,newSide)       = 
          !PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = 
          PartElemToElem(E2E_NB_ELEM_ID,PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))&
                                                                                                              =newElemID
          PartElemToElem(E2E_NB_LOC_SIDE_ID,PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))&
                                                                                                              =ilocSide
          PartElemToElem(E2E_NB_ELEM_ID,ilocSide,newElemID) = PartSideToElem(S2E_NB_ELEM_ID,newSideID)
          PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,newElemID) = PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)
        END IF
        isDone(haloSideID)=.TRUE.
      ELSE ! non-double side || new side
        ! cannnot build ElemToElem and PartElemToElem yet
        ! build PartSideToElem, so much as possible
        ! get correct side out of RecvMsg%SideToElem
        IF(iElem.EQ.RecvMsg%SideToElem(S2E_ELEM_ID,haloSideID))THEN
          ! master side
          PartSideToElem(S2E_ELEM_ID,newSideID)    =newElemID
          PartSideToElem(S2E_LOC_SIDE_ID,newSideID)=ilocSide
          PartSideToElem(S2E_FLIP       ,newSideID)=0 !RecvMsg%SideToElem(S2E_FLIP,haloSideID)
        ELSE IF(iElem.EQ.RecvMsg%SideToElem(S2E_NB_ELEM_ID,haloSideID))THEN
          ! slave side
          PartSideToElem(S2E_NB_ELEM_ID,newSideID) =newElemID
          PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) =ilocSide
          PartSideToElem(S2E_FLIP      ,newSideID) =RecvMsg%SideToElem(S2E_FLIP,haloSideID)
        ELSE ! should be found, because there should be halo sides without any connection
          CALL abort(&
__STAMP__&
             , 'Non-Critical error in domain reconstrution. IF NOT encountered, something is terrible wrong.')
        END IF
        !BC(1:4,newSideID)=RecvMsg%BC(1:4,haloSideID)
        BC(newSideID)=RecvMsg%BC(haloSideID)
        SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
      END IF ! isDoubleSide
      ! copy Bezier to new side id
      IF(.NOT.isDoubleSide) THEN
        BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newBCSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
        ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
        SideSlabNormals(1:3,1:3,newBCSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
        SideSlabIntervals(1:6,newBCSideID) =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
        BoundingBoxIsEmpty(newBCSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
      ELSE
        IF(RecvMsg%ElemToSide(2,ilocSide,iElem).EQ.0)THEN
          SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
          SideSlabIntervals(1:6,newSideID) =RecvMsg%SideSlabIntervals(1:6 ,haloSideID) 
          BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
        END IF
      END IF
      ! build entry to PartElemToSide
      BC(newSideID)=RecvMsg%BC(haloSideID)
      SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
      PartBCSideList(newSideID)=newBCSideID !tmpBCSides+haloinc(haloSideID)
      PartElemToSide(E2S_SIDE_ID,iLocSide,newElemId)=newSideID
      PartElemToSide(E2S_FLIP,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
    END DO ! ilocSide
    ! set native elemID
    PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
    PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
    XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
    dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
  END DO ! iElem
  ! missing connection info to MPI-Neighbor-ElemID
  DO iSide=nBCSides+nInnerSides+1,nSides
    ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
    ElemID2=PartSideToElem(S2E_NB_ELEM_ID,iSide)
    IF( ElemID .NE. -1 .AND. ElemID2 .NE. -1 ) CYCLE
    IF( ElemID .GE. 1 .AND. ElemID .LE. PP_nElems) THEN
      HostElemID=ElemID
      locSideID = SideToElem(S2E_LOC_SIDE_ID,iSide)
      flip=0
    END IF
    IF( ElemID2 .GE. 1 .AND. ElemID2 .LE. PP_nElems) THEN
      HostElemID=ElemID2
      locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,iSide)
      flip=1
    END IF
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      xNodes=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,HostelemID)
    CASE(XI_PLUS)
      xNodes=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,HostelemID)
    CASE(ETA_MINUS)
      xNodes=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,HostelemID)
    CASE(ETA_PLUS)
      xNodes=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,HostelemID)
    CASE(ZETA_MINUS)
      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,HostelemID)
    CASE(ZETA_PLUS)
      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,HostelemID)
    END SELECT
    IF(GEO%nPeriodicVectors.GT.0)THEN
      iDisplace= SidePeriodicType(iSide)
      IF(iDisplace.EQ.0) CYCLE
      IF(flip.EQ.0)THEN
        DO q=0,NGeo
          DO p=0,NGeo
            xNodes(1:3,q,p)=xNodes(1:3,q,p)+SidePeriodicDisplacement(1:3,iDisplace)
          END DO ! p=0,PP_N
        END DO ! q=0,PP_N
      ELSE
        DO q=0,NGeo
          DO p=0,NGeo
            xNodes(1:3,q,p)=xNodes(1:3,q,p)-SidePeriodicDisplacement(1:3,iDisplace)
          END DO ! p=0,PP_N
        END DO ! q=0,PP_N
      END IF
    END IF
    DO iElem=PP_nElems+1,nTotalElems
      DO ilocSide=1,6
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          xNodes2=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,iElem)
        CASE(XI_PLUS)
          xNodes2=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,iElem)
        CASE(ETA_MINUS)
          xNodes2=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,iElem)
        CASE(ETA_PLUS)
          xNodes2=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,iElem)
        CASE(ZETA_MINUS)
          xNodes2=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,iElem)
        CASE(ZETA_PLUS)
          xNodes2=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,iElem)
        END SELECT

        IF(  ALMOSTEQUAL(xNodes(1,0,0),xNodes2(1,0,0))               &
        .AND.ALMOSTEQUAL(xNodes(2,0,0),xNodes2(2,0,0))               &
        .AND.ALMOSTEQUAL(xNodes(3,0,0),xNodes2(3,0,0))               & 
        .AND.ALMOSTEQUAL(xNodes(1,NGeo,NGeo),xNodes2(1,NGeo,NGeo))   &
        .AND.ALMOSTEQUAL(xNodes(2,NGeo,NGeo),xNodes2(2,NGeo,NGeo))   &
        .AND.ALMOSTEQUAL(xNodes(3,NGeo,NGeo),xNodes2(3,NGeo,NGeo)) ) THEN
          IF(flip.EQ.0)THEN
            PartSideToElem(S2E_NB_ELEM_ID,iSide)=iElem
            EXIT
          ELSE
            PartSideToElem(S2E_ELEM_ID,iSide)=iElem
            EXIT
          END IF
        END IF
      END DO ! ilocSide=1,6
    END DO ! iElem=PP_nElems+1,nTotalElems
  END DO ! iSide=nBCSides+nInnerSides+1,nSides
  ! build rest: PartElemToElem, PartLocSideID
  DO iElem=PP_nElems+1,nTotalElems
    DO ilocSide=1,6
      flip   = PartElemToSide(E2S_FLIP,ilocSide,iElem)
      SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      ! check of sideid
      HaloSideID=SideID-tmpnSides
      !print*,'HaloSideID',HaloSideID
      ! do not double sides
      IF(HaloSideID.LE.tmpnSides) CYCLE 
      IF(isDone(HaloSideID)) CYCLE
      IF(flip.EQ.0)THEN
        ! SideID of slave
        PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)=PartSideToElem(S2E_NB_LOC_SIDE_ID,SideID)
        PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)=PartSideToElem(S2E_NB_ELEM_ID,SideID)
      ELSE
        ! SideID of master
        PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)=PartSideToElem(S2E_LOC_SIDE_ID,SideID)
        PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)=PartSideToElem(S2E_ELEM_ID,SideID)
      END IF
      isDone(HaloSideID)=.TRUE.
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

END SUBROUTINE ExchangeMappedHaloGeometry

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
