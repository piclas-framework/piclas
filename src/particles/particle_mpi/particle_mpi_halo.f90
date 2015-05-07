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

INTERFACE WriteParticleMappingPartitionInformation
  MODULE PROCEDURE WriteParticleMappingPartitionInformation
END INTERFACE


PUBLIC :: IdentifyHaloMPINeighborhood,ExchangeHaloGeometry,ExchangeMappedHaloGeometry
PUBLIC :: WriteParticlePartitionInformation,WriteParticleMappingPartitionInformation

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
USE MOD_Particle_Surfaces_Vars,     ONLY:BezierControlPoints3D,DoRefMapping
USE MOD_Mesh_Vars,                  ONLY:NGeo,ElemToSide,nElems,nInnerSides,nBCSides,nSides,SideToElem,XCL_NGeo
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
              IF(SideID.GT.(nInnerSides+nBCSides).AND.(SideIndex(SideID).EQ.0))THEN
              ! because of implicit, but here I send for checking, other process sends the required halo region
              !IF(SideIndex(SideID).EQ.0)THEN
                SendMsg%nMPISides=SendMsg%nMPISides+1
                SideIndex(SideID)=SendMsg%nMPISides
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
!IPWRITE(*,*)' nMPISides for iProc=',iProc,':',SendMsg%nMPISides

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
    CALL abort(__STAMP__,&
                         'Could not allocate SendMessage%BezierSides3D ',SendMsg%nMPISides)
  END IF
  SendMsg%BezierSides3D=0.
END IF
IF (RecvMsg%nMPISides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nMPISides), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__,&
                         'Could not allocate RecvMessage%BezierSides3D ',RecvMsg%nMPISides)
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
      IF(DoRefMapping)THEN
        SELECT CASE(iLocSide)
        CASE(XI_MINUS)
          SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,iElem)
          SideIndex(SideID)=0
        CASE(XI_PLUS)
          SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,iElem)
          SideIndex(SideID)=0
        CASE(ETA_MINUS)
          SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,iElem)
          SideIndex(SideID)=0
        CASE(ETA_PLUS)
          SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,iElem)
          SideIndex(SideID)=0
        CASE(ZETA_MINUS)
          SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,iElem)
          SideIndex(SideID)=0
        CASE(ZETA_PLUS)
          SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,iElem)
          SideIndex(SideID)=0
        END SELECT
      ELSE ! no ref mapping
        SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(SideID))=BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      END IF ! DoRefMapping
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
    CALL abort(__STAMP__,&
                         'Could not deallocate SendMessage%BezierSides3D proc ',iProc)
  END IF
END IF
IF (RecvMsg%nMPISides.GT.0) THEN
  DEALLOCATE(RecvMsg%BezierSides3D, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__,&
                         'Could not deallocate RecvMessage%BezierSides3D proc ',iProc)
  END IF
END IF

END SUBROUTINE IdentifyHaloMPINeighborhood


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
USE MOD_Mesh_Vars,                 ONLY:NGeo,ElemToSide,SideToElem,nSides,XCL_NGeo
USE MOD_Particle_Mesh_Vars,        ONLY:GEO, FIBGMCellPadding
USE MOD_Particle_MPI_Vars,         ONLY:halo_eps, NbrOfCases,casematrix, halo_eps2
USE MOD_Particle_Surfaces_Vars,    ONLY:BezierControlPoints3D,DoRefMapping
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
REAL                     :: NodeX(1:3),xNodes(1:3,0:NGeo,0:NGeo)
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
                IF(DoRefMapping)THEN
                  IF(SideIndex(SideID).EQ.0)THEN
                    leave=.FALSE.
                    SELECT CASE(ilocSide)
                    CASE(XI_MINUS)
                      xNodes=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,ElemID)
                    CASE(XI_PLUS)
                      xNodes=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,ElemID)
                    CASE(ETA_MINUS)
                      xNodes=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,ElemID)
                    CASE(ETA_PLUS)
                      xNodes=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,ElemID)
                    CASE(ZETA_MINUS)
                      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,ElemID)
                    CASE(ZETA_PLUS)
                      xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,ElemID)
                    END SELECT
                    leave=.FALSE.
                    DO s=0,NGeo
                      DO r=0,NGeo
                        IF(SQRT(DOT_Product(xNodes(:,r,s)-NodeX &
                                           ,xNodes(:,r,s)-NodeX )).LE.halo_eps)THEN
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
                ELSE ! norefmapping
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
                END IF
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
  WRITE(*,*) 'ups, something missing'
  STOP
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


!IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Number of marked sides:   ', NbOfSides

END SUBROUTINE CheckMPINeighborhoodByFIBGM


SUBROUTINE ExchangeHaloGeometry(iProc,SideList,ElemList)
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
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartHaloToProc
USE MOD_Mesh_Vars,              ONLY:nElems, nSides, nBCSides, nInnerSides, ElemToSide, BC,nGeo,SideToElem
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,nTotalSides,nTotalElems,SidePeriodicType
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartNeighborElemID,PartNeighborLocSideID
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars, ONLY:SlabNormals,SlabIntervalls,BoundingBoxIsEmpty
! should not be needed annymore
!USE MOD_Particle_MPI_Vars,      ONLY:nNbProcs,offsetMPISides_MINE, offsetMPISides_YOUR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iProc       ! MPI proc with which the local proc is to exchange boundary information
INTEGER, INTENT(INOUT)          :: SideList(nSides)
INTEGER, INTENT(INOUT)          :: ElemList(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE tMPISideMessage
  REAL,ALLOCATABLE          :: BezierControlPoints3D(:,:,:,:)
  INTEGER,ALLOCATABLE       :: ElemToSide(:,:,:) 
  INTEGER,ALLOCATABLE       :: SideToElem(:,:)
  INTEGER,ALLOCATABLE       :: SideBCType(:)
  INTEGER,ALLOCATABLE       :: BC(:)
  INTEGER,ALLOCATABLE       :: NativeElemID(:)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: SlabNormals                  ! normal vectors of bounding slab box
  REAL,ALLOCATABLE,DIMENSION(:,:)    :: SlabIntervalls               ! intervalls beta1, beta2, beta3
  LOGICAL,ALLOCATABLE,DIMENSION(:)   :: BoundingBoxIsEmpty
  !INTEGER,ALLOCATABLE       :: PeriodicElemSide(:,:)
  INTEGER                   :: nSides                 ! number of sides to send
  INTEGER                   :: nElems                 ! number of elems to send
END TYPE
TYPE(tMPISideMessage)       :: SendMsg
TYPE(tMPISideMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT
INTEGER                     :: newSideID,haloSideID,ioldSide,haloElemID,oldElemID,newElemID
LOGICAL                     :: isDoubleSide
LOGICAL,ALLOCATABLE         :: isElem(:),isSide(:),isDone(:)
INTEGER, ALLOCATABLE        :: ElemIndex(:), SideIndex(:),HaloInc(:)
INTEGER                     :: iElem, ilocSide,NbOfMarkedSides,SideID,iSide,p,q,iIndex,iHaloSide,flip
INTEGER                     :: nDoubleSides,nDoubleBezier,tmpnSides,tmpnElems
INTEGER                     :: datasize
!===================================================================================================================================

ALLOCATE(isElem(1:nElems))
IF (.NOT.ALLOCATED(isElem)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate isElem')
isElem(:) = .FALSE.

ALLOCATE(isSide(1:nSides))
IF (.NOT.ALLOCATED(isSide)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate isSide')
isSide(:) = .FALSE.

ALLOCATE(ElemIndex(1:nElems))
IF (.NOT.ALLOCATED(ElemIndex)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
ElemIndex(:) = 0

ALLOCATE(SideIndex(1:nSides))
IF (.NOT.ALLOCATED(SideIndex)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate SideIndex')
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
      IF(.NOT.isSide(SideID)) THEN
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
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%ElemToSide',SendMsg%nElems)
  SendMsg%ElemToSide(:,:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemToSide(1:2,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
  RecvMsg%ElemToSide(:,:,:)=0
END IF

! BezierControlPoints3D for exchange
IF (SendMsg%nSides.GT.0) THEN       ! Beziercontrolpoints3d
  ALLOCATE(SendMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%BezierControlPoints3D',SendMsg%nSides)
  SendMsg%BezierControlPoints3D=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%BezierControlPoints3D',RecvMsg%nSides)
  RecvMsg%BezierControlPoints3D=0.
END IF

! SideToElem Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SideToElem(1:2,1:nSides) 
  ALLOCATE(SendMsg%SideToElem(1:5,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SideToElem',SendMsg%nSides)
  SendMsg%SideToElem(:,:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideToElem(1:5,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SideToElem',RecvMsg%nSides)
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
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%BC',SendMsg%nSides,999.)
  SendMsg%BC(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BC(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%BC',RecvMsg%nSides,999.)
  RecvMsg%BC(:)=0
END IF
! SideBCType
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%SideBCType(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SideBCType',SendMsg%nSides,999.)
  SendMsg%SideBCType(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideBCType(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SideBCType',RecvMsg%nSides,999.)
  RecvMsg%SideBCType(:)=0
END IF
! NativeElemID 
IF (SendMsg%nElems.GT.0) THEN 
  ALLOCATE(SendMsg%NativeElemID(1:SendMsg%nElems),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%NativeElemID',SendMsg%nElems,999.)
  SendMsg%NativeElemID(:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%NativeElemID(1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%NativeElemID',RecvMsg%nElems,999.)
  RecvMsg%NativeElemID(:)=0
END IF
! SlabNormals Mapping
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%SlabNormals(1:3,1:3,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SlabNormals',SendMsg%nSides)
  SendMsg%SlabNormals(:,:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SlabNormals(1:3,1:3,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SlabNormals',RecvMsg%nSides)
  RecvMsg%SlabNormals(:,:,:)=0.
END IF
! SlabIntervalls Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SlabIntervalls(1:2,1:nSides) 
  ALLOCATE(SendMsg%SlabIntervalls(1:6,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SlabIntervalls',SendMsg%nSides)
  SendMsg%SlabIntervalls(:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SlabIntervalls(1:6,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SlabIntervalls',RecvMsg%nSides)
  RecvMsg%SlabIntervalls(:,:)=0.
END IF
! BoundingBoxIsEmpty Mapping
IF (SendMsg%nSides.GT.0) THEN       ! BoundingBoxIsEmpty(1:2,1:nSides) 
  ALLOCATE(SendMsg%BoundingBoxIsEmpty(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%BoundingBoxIsEmpty',SendMsg%nSides)
  SendMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BoundingBoxIsEmpty(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%BoundingBoxIsEmpty',RecvMsg%nSides)
  RecvMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
!! PeriodicElemSide Mapping 
!IF (SendMsg%nElems.GT.0) THEN       ! PeriodicElemSide(1:iLocSide,1:nElems)
!  ALLOCATE(SendMsg%PeriodicElemSide(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
!    'Could not allocate SendMsg%PeriodicElemSide',SendMsg%nElems,999.)
!  SendMsg%PeriodicElemSide(:,:)=0
!END IF
!IF (RecvMsg%nElems.GT.0) THEN
!  ALLOCATE(RecvMsg%PeriodicElemSide(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
!    'Could not allocate RecvMsg%PeriodicElemSide',RecvMsg%nElems,999.)
!  RecvMsg%PeriodicElemSide(:,:)=0
!END IF


! fill send buffers with node, side and element data (including connectivity!)
! ElemtoSide 
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
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
                 ElemIndex(SideToElem(iIndex,iSide))     
         SendMsg%SideBCType(SideIndex(iSide)) = SidePeriodicType(iSide)
      END IF
    END DO ! S2E_LOC_SIDE_ID, S2E_NB_LOC_SIDE_ID, S2E_FLIP
   ! IF(PartMPI%MPIROOT) print*,'iSide,SideIndex',iSide,SideIndex(iSide)
    SendMsg%SideToElem(3:5,SideIndex(iSide)) = &
        SideToElem(3:5,iSide)
    SendMsg%BezierControlPoints3D(:,:,:,SideIndex(iSide)) = &
        BezierControlPoints3D(:,:,:,iSide)
    ! slabnormals
    SendMsg%SlabNormals(:,:,SideIndex(iSide)) = &
        SlabNormals(:,:,iSide)
    ! slabintervalls
    SendMsg%SlabIntervalls(:,SideIndex(iSide)) = &
        SlabIntervalls(:,iSide)
    ! BoundingBoxIsEmpty
    SendMsg%BoundingBoxIsEmpty(SideIndex(iSide)) = &
        BoundingBoxIsEmpty(iSide)
  END IF
END DO
!--- BC Mapping ------------------------------------------------------!
DO iSide = 1,nBCSides  ! no need to go through all side since BC(1:nBCSides)
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
!IPWRITE(*,*) " Now MPI exchange"

dataSize=3*(NGeo+1)*(NGeo+1)
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
      CALL MPI_SEND(SendMsg%SlabNormals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SlabIntervalls,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)

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
      CALL MPI_RECV(RecvMsg%SlabNormals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SlabIntervalls,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
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
      CALL MPI_RECV(RecvMsg%SlabNormals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SlabIntervalls,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)

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
      CALL MPI_SEND(SendMsg%SlabNormals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SlabIntervalls,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
END IF

IF ((RecvMsg%nElems.EQ.0) .AND. (RecvMsg%nSides.GT.0))THEN
    ERRWRITE(*,*)'ERROR: nElems=0 when nSides=',RecvMsg%nSides,' and nSides=',RecvMsg%nSides,'!'
    CALL abort(&
   __STAMP__,&
      'nElems=0 while nSides=',RecvMsg%nSides)
END IF

DEALLOCATE(isElem,isSide,ElemIndex,SideIndex)

!print*,'Rank,sendsides',PartMPI%MyRank,SendMsg%nSides
!print*,'Rank,recvsides',PartMPI%MyRank,RecvMsg%nSides
!print*,'Rank,recvelemtoside',PartMPI%MyRank,RecvMsg%ElemToSide
!print*,'iproc',iproc

!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!IPWRITE(*,*) " Recive stuff"

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
        IF(HaloInc(haloSideID).EQ.0) IPWRITE(*,*) ' Warning: wrong halo inc'
        newSideID=tmpnSides+haloinc(haloSideID)
        IF(newSideID.LT.tmpnSides) IPWRITE(*,*) 'Warning: wrong new sideid', newsideid
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
            CALL abort(__STAMP__,&
            'Critical error in domain reconstrution.')
            PartSideToElem(S2E_NB_ELEM_ID,newSideID)     = newElemID
            PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = ilocSide
            ! here? CAUTION// DEBUG
            PartSideToElem(S2E_FLIP          ,newSideID) = RecvMsg%SideToElem(S2E_FLIP,haloSideID)
            ! nothing to do, is already filled
            !PartSideToElem(S2E_ELEM_ID       ,newSideID) = 
            !PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = 
            !PartSideToElem(S2E_FLIP          ,newSideID) = 
            ! NeighboreElemID
            PartneighborElemID   (PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))=newElemID
            PartNeighborlocSideID(PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))=ilocSide
            PartNeighborElemID(ilocSide,newElemID)    = PartSideToElem(S2E_ELEM_ID,newSideID)
            PartNeighborlocSideID(ilocSide,newElemID) = PartSideToElem(S2E_LOC_SIDE_ID,newSideID)
          !  IF(PartSideToElem(S2E_ELEM_ID,newSideID).EQ.-1) IPWRITE(*,*) 'warning'
        ELSE ! SE2_NB_ELEM_ID=DEFINED
          IF(PartSideToElem(S2E_ELEM_ID,newSideID).NE.-1) &
            CALL abort(__STAMP__,&
            'Critical error in domain reconstrution.')
          PartSideToElem(S2E_ELEM_ID       ,newSideID) = newElemID !root Element
          PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = iLocSide
          !PartSideToElem(S2E_FLIP          ,newSideID) = RecvMsg%SideToElem(S2E_FLIP
          ! already filled
          !PartSideToElem(S2E_NB_ELEM_ID,newSide)       = 
          !PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = 
          PartNeighborElemID(   PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))=newElemID
          PartNeighborlocSideID(PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))=ilocSide
          PartNeighborElemID(ilocSide,newElemID)    = PartSideToElem(S2E_NB_ELEM_ID,newSideID)
          PartNeighborlocSideID(ilocSide,newElemID) = PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)
          !IF(PartSideToElem(S2E_NB_ELEM_ID,newSideID).EQ.-1) IPWRITE(*,*)'warning'
        END IF
        isDone(haloSideID)=.TRUE.
      ELSE ! non-double side || new side
        ! cannnot build PartNeighborElemID and PartNeighborlocSideID yet
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
          CALL abort(__STAMP__,&
              'Non-Critical error in domain reconstrution. IF NOT encountered, something is terrible wrong.')
        END IF
        !BC(1:4,newSideID)=RecvMsg%BC(1:4,haloSideID)
        BC(newSideID)=RecvMsg%BC(haloSideID)
        SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
      END IF ! isDoubleSide
      ! copy Bezier to new side id
      IF(.NOT.isDoubleSide) THEN
        BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
        ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
        SlabNormals(1:3,1:3,newSideID)=RecvMsg%SlabNormals(1:3,1:3,haloSideID)
        SlabIntervalls(1:6,newSideID) =RecvMsg%SlabIntervalls(1:6 ,haloSideID) 
        BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
      ELSE
        IF(RecvMsg%ElemToSide(2,ilocSide,iElem).EQ.0)THEN
          SlabNormals(1:3,1:3,newSideID)=RecvMsg%SlabNormals(1:3,1:3,haloSideID)
          SlabIntervalls(1:6,newSideID) =RecvMsg%SlabIntervalls(1:6 ,haloSideID) 
          BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
        END IF
      END IF
      ! build entry to PartElemToSide
      PartElemToSide(1,iLocSide,newElemId)=newSideID
      PartElemToSide(2,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
    END DO ! ilocSide
    ! set native elemID
    PartHaloToProc(1,newElemId)=RecvMsg%NativeElemID(iElem)
    PartHaloToProc(2,newElemId)=iProc
  END DO ! iElem
  ! build rest: PartNeighborElemID, PartLocSideID
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
        PartneighborlocSideID(ilocSide,iElem)=PartSideToElem(S2E_NB_LOC_SIDE_ID,SideID)
        PartneighborElemID   (ilocSide,iElem)=PartSideToElem(S2E_NB_ELEM_ID,SideID)
      ELSE
        ! SideID of master
        PartneighborlocSideID(ilocSide,iElem)=PartSideToElem(S2E_LOC_SIDE_ID,SideID)
        PartneighborElemID   (ilocSide,iElem)=PartSideToElem(S2E_ELEM_ID,SideID)
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

END SUBROUTINE ExchangeHaloGeometry


SUBROUTINE ResizeParticleMeshData(nOldSides,nOldElems,nTotalSides,nTotalElems,nOldBCSides,nTotalBCSides)
!===================================================================================================================================
! resize the partilce mesh data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartHaloToProc
USE MOD_Mesh_Vars,              ONLY:BC,nGeo,nElems,XCL_NGeo,DXCL_NGEO
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartNeighborElemID,PartNeighborLocSideID
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,DoRefMapping
USE MOD_Particle_Surfaces_Vars, ONLY:SlabNormals,SlabIntervalls,BoundingBoxIsEmpty
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
INTEGER,ALLOCATABLE                :: DummyHaloToProc(:,:)                                 
INTEGER,ALLOCATABLE                :: DummySideToElem(:,:)
INTEGER,ALLOCATABLE                :: DummySideBCType(:),DummyPartBCSideList(:)
INTEGER,ALLOCATABLE                :: DummyNeighborElemID(:,:)
INTEGER,ALLOCATABLE                :: DummyNeighborlocSideID(:,:)
REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummySlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummySlabIntervalls               ! intervalls beta1, beta2, beta3
LOGICAL,ALLOCATABLE,DIMENSION(:)   :: DummyBoundingBoxIsEmpty
!===================================================================================================================================

! reallocate shapes
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!print*,'rank, in reshape',myrank
!print*,'rank, in reshape',myrank, nOldSides,nOldElems,nTotalSides,nTotalElems
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
! PartElemToSide
ALLOCATE(DummyElemToSide(1:2,1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemToSide)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummyElemToSide=PartElemToSide
!IPWRITE(*,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
DEALLOCATE(PartElemToSide)
ALLOCATE(PartElemToSide(1:2,1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate PartElemToSide')
PartElemToSide=-1
PartElemToSide(:,:,1:nOldElems) =DummyElemToSide(:,:,1:nOldElems)
DEALLOCATE(DummyElemToSide)
IF(DoRefMapping)THEN
  ! XCL_NGeo
  ALLOCATE(DummyXCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyXCL_NGeo)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummyXCL_NGeo=XCL_NGeo
  !IPWRITE(*,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
  DEALLOCATE(XCL_NGeo)
  ALLOCATE(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate XCL_NGeo')
  XCL_NGeo=0
  XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems) =DummyXCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems)
  DEALLOCATE(DummyXCL_NGeo)
  ! dXCL_NGeo
  ALLOCATE(DummydXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems))
  IF (.NOT.ALLOCATED(DummydXCL_NGeo)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummydXCL_NGeo=dXCL_NGeo
  !IPWRITE(*,*)"not allocated partelemtoside",ALLOCATED(PartElemToSide)
  DEALLOCATE(dXCL_NGeo)
  ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate dXCL_NGeo')
  dXCL_NGeo=0
  dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems) =DummydXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems)
  DEALLOCATE(DummydXCL_NGeo)
  ! PartBCSideList
  ALLOCATE(DummyPartBCSideList(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyPartBCSideList)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate PartBCSideList')
  DummyPartBCSideList(1:nOldSides)=PartBCSideList(1:nOldSides)
  DEALLOCATE(PartBCSideList)
  ALLOCATE(PartBCSideList(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  PartBCSideList=HUGE(1)
  PartBCSideList(1:nOldSides) =DummyPartBCSideList(1:nOldSides)
  DEALLOCATE(DummyPartBCSideList)
END IF

!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!print*,' done elem to side',myrank
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
! HaloToProc
IF(.NOT.ALLOCATED(PartHaloToProc))THEN
  !print*,'both here',myrank,nElems,nElems+1,nTotalElems
  !print*,'myrank',myrank,allocstat
  nLower=nElems+1
  ALLOCATE(PartHaloToProc(1:3,nLower:nTotalElems),STAT=ALLOCSTAT)                                 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate PartHaloToProc')
  PartHaloToProc=-1
  !print*,'lower,upper',PP_nElems+1,nTotalElems
ELSE
  nLower=nElems+1
  ALLOCATE(DummyHaloToProc(1:3,nLower:nOldElems))                                 
  IF (.NOT.ALLOCATED(DummyHaloToProc)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate DummyPartHaloToProc')
  DummyHaloToProc=PartHaloToProc
  DEALLOCATE(PartHaloToProc)
  ALLOCATE(PartHaloToProc(1:3,nLower:nTotalElems),STAT=ALLOCSTAT)                                 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate PartHaloToProc')
  ! copy array to new
  PartHaloToProc=-1
  PartHaloToProc(1:2,PP_nElems+1:nOldElems)    =DummyHaloToProc(1:2,PP_nElems+1:nOldElems)
  DEALLOCATE(DummyHaloToProc)
END IF
!print*,' done halotoproc',myrank
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!print*,' done halotoproc',myrank
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
! PartSideToElem
ALLOCATE(DummySideToElem(1:5,1:nOldSides))
IF (.NOT.ALLOCATED(DummySideToElem)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummySideToElem=PartSideToElem
DEALLOCATE(PartSideToElem)
ALLOCATE(PartSideToElem(1:5,1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate PartSideToElem')
PartSideToElem=-1
PartSideToElem(:,1:nOldSides  )              =DummySideToElem(:,1:nOldSides)
DEALLOCATE(DummySideToElem)
!print*,' done side to elem',myrank
! PartNeighborElemID
ALLOCATE(DummyNeighborElemID(1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyNeighborElemID)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummyNeighborElemID=PartNeighborElemID
DEALLOCATE(PartNeighborElemID)
ALLOCATE(PartNeighborElemID(1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
PartNeighborElemID=-1
PartNeighborElemID(:,1:nOldElems)            =DummyNeighborElemID(:,1:nOldElems)
DEALLOCATE(DummyNeighborElemID)
! PartNeighborlocSideID
ALLOCATE(DummyNeighborlocSideID(1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyNeighborLocSideID)) CALL abort(&
    __STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummyNeighborlocSideID=PartNeighborlocSideID
DEALLOCATE(PartNeighborlocSideID)
ALLOCATE(PartNeighborlocSideID(1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
PartNeighborlocSideID=-1
PartNeighborlocSideID(:,1:nOldElems)         =DummyNeighborlocSideID(:,1:nOldElems)
DEALLOCATE(DummyNeighborlocSideID)
IF(DoRefMapping)THEN
  ! BezierControlPoints3D
  ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummyBezierControlPoints3d)) CALL abort(&
      __STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummyBezierControlPoints3d=BezierControlPoints3d
  DEALLOCATE(BezierControlPoints3D)
  ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  BezierControlPoints3d(:,:,:,1:nOldBCSides) =DummyBezierControlPoints3D(:,:,:,1:nOldBCSides)
  DEALLOCATE(DummyBezierControlPoints3D)
  ! SlabNormals
  ALLOCATE(DummySlabNormals(1:3,1:3,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySlabNormals)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummySlabNormals=SlabNormals
  DEALLOCATE(SlabNormals)
  ALLOCATE(SlabNormals(1:3,1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  SlabNormals=0
  SlabNormals(1:3,1:3,1:nOldBCSides) =DummySlabNormals(1:3,1:3,1:nOldBCSides)
  DEALLOCATE(DummySlabNormals)
  ! SlabIntervalls
  ALLOCATE(DummySlabIntervalls(1:6,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySlabIntervalls)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummySlabIntervalls=SlabIntervalls
  DEALLOCATE(SlabIntervalls)
  ALLOCATE(SlabIntervalls(1:6,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  SlabIntervalls=0
  SlabIntervalls(1:6,1:nOldBCSides) =DummySlabIntervalls(1:6,1:nOldBCSides)
  DEALLOCATE(DummySlabIntervalls)
  ! BoundingBoxIsEmpty
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
      __STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
  DEALLOCATE(BoundingBoxIsEmpty)
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  BoundingBoxIsEmpty(1:nOldBCSides) =DummyBoundingBoxIsEmpty(1:nOldBCSides)
  DEALLOCATE(DummyBoundingBoxIsEmpty)

ELSE ! no mapping
  ! BezierControlPoints3D
  ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBezierControlPoints3d)) CALL abort(&
      __STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummyBezierControlPoints3d=BezierControlPoints3d
  DEALLOCATE(BezierControlPoints3D)
  ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  BezierControlPoints3d(:,:,:,1:nOldSides) =DummyBezierControlPoints3D(:,:,:,1:nOldSides)
  DEALLOCATE(DummyBezierControlPoints3D)
  ! SlabNormals
  ALLOCATE(DummySlabNormals(1:3,1:3,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySlabNormals)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummySlabNormals=SlabNormals
  DEALLOCATE(SlabNormals)
  ALLOCATE(SlabNormals(1:3,1:3,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  SlabNormals=0
  SlabNormals(1:3,1:3,1:nOldSides) =DummySlabNormals(1:3,1:3,1:nOldSides)
  DEALLOCATE(DummySlabNormals)
  ! SlabIntervalls
  ALLOCATE(DummySlabIntervalls(1:6,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySlabIntervalls)) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummySlabIntervalls=SlabIntervalls
  DEALLOCATE(SlabIntervalls)
  ALLOCATE(SlabIntervalls(1:6,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  SlabIntervalls=0
  SlabIntervalls(1:6,1:nOldSides) =DummySlabIntervalls(1:6,1:nOldSides)
  DEALLOCATE(DummySlabIntervalls)
  ! BoundingBoxIsEmpty
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
      __STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
  DEALLOCATE(BoundingBoxIsEmpty)
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
    'Could not allocate ElemIndex')
  BoundingBoxIsEmpty(1:nOldSides) =DummyBoundingBoxIsEmpty(1:nOldSides)
  DEALLOCATE(DummyBoundingBoxIsEmpty)
END IF
! SideBCType
ALLOCATE(DummySideBCType(1:nOldSides))
IF (.NOT.ALLOCATED(DummySideBCType)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummySideBCType(1:nOldSides)=SidePeriodicType(1:nOldSides)
DEALLOCATE(SidePeriodicType)
ALLOCATE(SidePeriodicType(1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
SidePeriodicType=-1
SidePeriodicType(1:nOldSides) =DummySideBCType(1:nOldSides)
DEALLOCATE(DummySideBCType)
! BC
ALLOCATE(DummyBC(1:nOldSides))
IF (.NOT.ALLOCATED(DummyBC)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
! check
!IF(ALLOCATED(BC)) print*,'yes it is'
!print*,'size',size(BC)
!print*, BC
DummyBC(1:nOldSides)=BC(1:nOldSides)
DEALLOCATE(BC)
ALLOCATE(BC(1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
BC=0
BC(1:nOldSides) =DummyBC(1:nOldSides)
DEALLOCATE(DummyBC)
! finished copying

END SUBROUTINE ResizeParticleMeshData


SUBROUTINE ExchangeMappedHaloGeometry(iProc,SideList,ElemList)
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
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartHaloToProc
USE MOD_Mesh_Vars,              ONLY:nElems, nSides, nBCSides, nInnerSides, ElemToSide, BC,nGeo,SideToElem
USE MOD_Mesh_Vars,              ONLY:XCL_NGeo,dXCL_NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,nTotalSides,nTotalElems,SidePeriodicType,PartBCSideList,nTotalBCSides
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartNeighborElemID,PartNeighborLocSideID
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars, ONLY:SlabNormals,SlabIntervalls,BoundingBoxIsEmpty
! should not be needed annymore
!USE MOD_Particle_MPI_Vars,      ONLY:nNbProcs,offsetMPISides_MINE, offsetMPISides_YOUR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iProc       ! MPI proc with which the local proc is to exchange boundary information
INTEGER, INTENT(INOUT)          :: SideList(nSides)
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
  REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: SlabNormals                  ! normal vectors of bounding slab box
  REAL,ALLOCATABLE,DIMENSION(:,:)    :: SlabIntervalls               ! intervalls beta1, beta2, beta3
  LOGICAL,ALLOCATABLE,DIMENSION(:)   :: BoundingBoxIsEmpty
  !INTEGER,ALLOCATABLE       :: PeriodicElemSide(:,:)
  INTEGER                   :: nSides                 ! number of sides to send
  INTEGER                   :: nElems                 ! number of elems to send
END TYPE
TYPE(tMPISideMessage)       :: SendMsg
TYPE(tMPISideMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT,tmpBCSides,newBCSideID
INTEGER                     :: newSideID,haloSideID,ioldSide,haloElemID,oldElemID,newElemID
LOGICAL                     :: isDoubleSide
LOGICAL,ALLOCATABLE         :: isElem(:),isSide(:),isDone(:)
INTEGER, ALLOCATABLE        :: ElemIndex(:), SideIndex(:),HaloInc(:)
INTEGER                     :: iElem, ilocSide,NbOfMarkedSides,SideID,iSide,p,q,iIndex,iHaloSide,flip
INTEGER                     :: nDoubleSides,nDoubleBezier,tmpnSides,tmpnElems
INTEGER                     :: datasize,datasize2,datasize3
INTEGER                     :: markElem
!===================================================================================================================================

ALLOCATE(isElem(1:nElems))
IF (.NOT.ALLOCATED(isElem)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate isElem')
isElem(:) = .FALSE.

ALLOCATE(isSide(1:nSides))
IF (.NOT.ALLOCATED(isElem)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate isSide')
isSide(:) = .FALSE.

ALLOCATE(ElemIndex(1:nElems))
IF (.NOT.ALLOCATED(ElemIndex)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
ElemIndex(:) = 0

ALLOCATE(SideIndex(1:nSides))
IF (.NOT.ALLOCATED(SideIndex)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate SideIndex')
SideIndex(:) = 0

!--- First, count marker node indices (nNeighborhoodNodes are within eps distance of at least one MPI-neighbor's node)
!--- For each MPI neighbor, identify the number of sides and elements to be sent
SendMsg%nElems=0
SendMsg%nSides=0
!LOGWRITE(*,*)'nNeighborhoodNodes=',nNeighborhoodNodes

! 1) get number of elements and sides
DO iElem=1,nElems
  SWRITE(*,*) ElemList(iElem)
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
      !IF(.NOT.isSide(SideID).AND.SideID.LE.nBCSides) THEN
      IF(.NOT.isSide(SideID).AND.PartBCSideList(SideID).LE.nTotalBCSides) THEN
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
!IPWRITE(*,*) ' Send number of elems,sides', SendMsg%nElems, SendMsg%nSides

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
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%ElemToSide',SendMsg%nElems)
  SendMsg%ElemToSide(:,:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemToSide(1:2,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
  RecvMsg%ElemToSide(:,:,:)=0
END IF

! BezierControlPoints3D for exchange
IF (SendMsg%nSides.GT.0) THEN       ! Beziercontrolpoints3d
  ALLOCATE(SendMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%BezierControlPoints3D',SendMsg%nSides)
  SendMsg%BezierControlPoints3D=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%BezierControlPoints3D',RecvMsg%nSides)
  RecvMsg%BezierControlPoints3D=0.
END IF

! XCL_NGeo for exchange
IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%XCL_NGeo',SendMsg%nElems)
  SendMsg%XCL_NGeo(:,:,:,:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%XCL_NGeo',RecvMsg%nElems)
  RecvMsg%XCL_NGeo(:,:,:,:,:)=0
END IF

! DXCL_NGeo for exchange
IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%DXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%DXCL_NGeo',SendMsg%nElems)
  SendMsg%DXCL_NGeo(:,:,:,:,:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%DXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
  RecvMsg%DXCL_NGeo(:,:,:,:,:,:)=0
END IF

! SideToElem Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SideToElem(1:2,1:nSides) 
  ALLOCATE(SendMsg%SideToElem(1:5,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SideToElem',SendMsg%nSides)
  SendMsg%SideToElem(:,:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideToElem(1:5,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SideToElem',RecvMsg%nSides)
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
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%BC',SendMsg%nSides,999.)
  SendMsg%BC(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BC(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%BC',RecvMsg%nSides,999.)
  RecvMsg%BC(:)=0
END IF
! SideBCType
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%SideBCType(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SideBCType',SendMsg%nSides,999.)
  SendMsg%SideBCType(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideBCType(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SideBCType',RecvMsg%nSides,999.)
  RecvMsg%SideBCType(:)=0
END IF
! NativeElemID 
IF (SendMsg%nElems.GT.0) THEN 
  ALLOCATE(SendMsg%NativeElemID(1:SendMsg%nElems),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%NativeElemID',SendMsg%nElems,999.)
  SendMsg%NativeElemID(:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%NativeElemID(1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%NativeElemID',RecvMsg%nElems,999.)
  RecvMsg%NativeElemID(:)=0
END IF
! SlabNormals Mapping
IF (SendMsg%nSides.GT.0) THEN       
  ALLOCATE(SendMsg%SlabNormals(1:3,1:3,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SlabNormals',SendMsg%nSides)
  SendMsg%SlabNormals(:,:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SlabNormals(1:3,1:3,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SlabNormals',RecvMsg%nSides)
  RecvMsg%SlabNormals(:,:,:)=0.
END IF
! SlabIntervalls Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SlabIntervalls(1:2,1:nSides) 
  ALLOCATE(SendMsg%SlabIntervalls(1:6,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%SlabIntervalls',SendMsg%nSides)
  SendMsg%SlabIntervalls(:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SlabIntervalls(1:6,1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%SlabIntervalls',RecvMsg%nSides)
  RecvMsg%SlabIntervalls(:,:)=0.
END IF
! BoundingBoxIsEmpty Mapping
IF (SendMsg%nSides.GT.0) THEN       ! BoundingBoxIsEmpty(1:2,1:nSides) 
  ALLOCATE(SendMsg%BoundingBoxIsEmpty(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h 
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate SendMsg%BoundingBoxIsEmpty',SendMsg%nSides)
  SendMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BoundingBoxIsEmpty(1:RecvMsg%nSides),STAT=ALLOCSTAT)  
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
    'Could not allocate RecvMsg%BoundingBoxIsEmpty',RecvMsg%nSides)
  RecvMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
!! PeriodicElemSide Mapping 
!IF (SendMsg%nElems.GT.0) THEN       ! PeriodicElemSide(1:iLocSide,1:nElems)
!  ALLOCATE(SendMsg%PeriodicElemSide(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
!    'Could not allocate SendMsg%PeriodicElemSide',SendMsg%nElems,999.)
!  SendMsg%PeriodicElemSide(:,:)=0
!END IF
!IF (RecvMsg%nElems.GT.0) THEN
!  ALLOCATE(RecvMsg%PeriodicElemSide(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
!  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,&
!    'Could not allocate RecvMsg%PeriodicElemSide',RecvMsg%nElems,999.)
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
    SendMsg%SlabNormals(:,:,SideIndex(iSide)) = &
        SlabNormals(:,:,iSide)
    ! slabintervalls
    SendMsg%SlabIntervalls(:,SideIndex(iSide)) = &
        SlabIntervalls(:,iSide)
    ! BoundingBoxIsEmpty
    SendMsg%BoundingBoxIsEmpty(SideIndex(iSide)) = &
        BoundingBoxIsEmpty(iSide)
  END IF
END DO
!--- BC Mapping ------------------------------------------------------!
DO iSide = 1,nBCSides  ! no need to go through all side since BC(1:nBCSides)
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
!IPWRITE(*,*) " Now MPI exchange"

dataSize =3*(NGeo+1)*(NGeo+1)
dataSize2=3*(NGeo+1)*(NGeo+1)*(NGeo+1)
dataSize3=9*(NGeo+1)*(NGeo+1)*(NGeo+1)
IF (PartMPI%MyRank.LT.iProc) THEN
  ! Send:
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BezierControlPoints3D,SendMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
!  IF(GEO%nPeriodicVectors.GT.0)THEN
!    IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%PeriodicElemSide,SendMsg%nElems*6,MPI_INTEGER ,iProc,1109,PartMPI%COMM,IERROR)
!  END IF
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideBCType,SendMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SlabNormals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SlabIntervalls,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror

  ! Receive:
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
  CALL MPI_RECV(RecvMsg%BezierControlPoints3D,RecvMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror

  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,MPISTATUS,IERROR)
  !!IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  !IF(GEO%nPeriodicVectors.GT.0)THEN
  !  IF (RecvMsg%nElems.GT.0) &
  !       CALL MPI_RECV(RecvMsg%PeriodicElemSide,RecvMsg%nElems*6,MPI_INTEGER,iProc,1109,PartMPI%COMM,MPISTATUS,IERROR)
  !END IF
  IF (RecvMsg%nSides.GT.0) CALL MPI_RECV(RecvMsg%SideBCType,RecvMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SlabNormals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SlabIntervalls,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror

  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  ! Receive:
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
  CALL MPI_RECV(RecvMsg%BezierControlPoints3D,RecvMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  !IF(GEO%nPeriodicVectors.GT.0)THEN
  !  IF (RecvMsg%nElems.GT.0) &
  !       CALL MPI_RECV(RecvMsg%PeriodicElemSide,RecvMsg%nElems*6,MPI_INTEGER,iProc,1109,PartMPI%COMM,MPISTATUS,IERROR)
  !END IF
  IF (RecvMsg%nSides.GT.0) CALL MPI_RECV(RecvMsg%SideBCType,RecvMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SlabNormals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SlabIntervalls,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
  !IPWRITE(*,*) 'ierror',ierror


  ! Send:
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BezierControlPoints3D,SendMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  !IF(GEO%nPeriodicVectors.GT.0)THEN
  !  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%PeriodicElemSide,SendMsg%nElems*6,MPI_INTEGER ,iProc,1109,PartMPI%COMM,IERROR)
  !END IF
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideBCType,SendMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SlabNormals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SlabIntervalls,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)
  !IPWRITE(*,*) 'ierror',ierror

END IF

IF ((RecvMsg%nElems.EQ.0) .AND. (RecvMsg%nSides.GT.0))THEN
    ERRWRITE(*,*)'ERROR: nElems=0 when nSides=',RecvMsg%nSides,' and nSides=',RecvMsg%nSides,'!'
    CALL abort(&
   __STAMP__,&
      'nElems=0 while nSides=',RecvMsg%nSides)
END IF

DEALLOCATE(isElem,isSide,ElemIndex,SideIndex)

!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!IPWRITE(*,*) " Recive stuff"

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
  !IPWRITE(*,*) 'recnsides',RecvMsg%nSides
  !IPWRITE(*,*) 'nSides,nnewSides,nDoubleSides', nSides,RecvMsg%nSides,nDoubleSides
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
      !IPWRITE(*,*) 'halosideid', haloSideID
      IF(haloSideID.EQ.0) CYCLE ! all non BC faces have not to be located
      isDoubleSide=.FALSE.
      IF(isSide(haloSideID)) THEN
        IF(HaloInc(haloSideID).EQ.0) IPWRITE(*,*) ' Warning: wrong halo inc'
        newSideID  =tmpnSides+haloinc(haloSideID)
        newBCSideID=tmpBCSides+haloinc(haloSideID)
        IF(newSideID.LT.tmpnSides) IPWRITE(*,*) 'Warning: wrong new sideid', newsideid
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
            CALL abort(__STAMP__,&
            'Critical error in domain reconstrution.')
            PartSideToElem(S2E_NB_ELEM_ID,newSideID)     = newElemID
            PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = ilocSide
            ! nothing to do, is already filled
            !PartSideToElem(S2E_ELEM_ID       ,newSideID) = 
            !PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = 
            !PartSideToElem(S2E_FLIP          ,newSideID) = 
            ! NeighboreElemID
            PartneighborElemID(PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))=newElemID
            PartNeighborlocSideID(PartSideToElem(S2E_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_ELEM_ID,newSideID))=ilocSide
            PartNeighborElemID(ilocSide,newElemID)    = PartSideToElem(S2E_ELEM_ID,newSideID)
            PartNeighborlocSideID(ilocSide,newElemID) = PartSideToElem(S2E_LOC_SIDE_ID,newSideID)
        ELSE ! SE2_NB_ELEM_ID=DEFINED
          IF(PartSideToElem(S2E_ELEM_ID,newSideID).NE.-1) &
            CALL abort(__STAMP__,&
            'Critical error in domain reconstrution.')
          PartSideToElem(S2E_ELEM_ID       ,newSideID) = newElemID !root Element
          PartSideToElem(S2E_LOC_SIDE_ID   ,newSideID) = iLocSide
          PartSideToElem(S2E_FLIP          ,newSideID) = 0
          ! already filled
          !PartSideToElem(S2E_NB_ELEM_ID,newSide)       = 
          !PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) = 
          PartNeighborElemID(PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))=newElemID
          PartNeighborlocSideID(PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID),PartSideToElem(S2E_NB_ELEM_ID,newSideID))=ilocSide
          PartNeighborElemID(ilocSide,newElemID)    = PartSideToElem(S2E_NB_ELEM_ID,newSideID)
          PartNeighborlocSideID(ilocSide,newElemID) = PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)
        END IF
        isDone(haloSideID)=.TRUE.
      ELSE ! non-double side || new side
        ! cannnot build PartNeighborElemID and PartNeighborlocSideID yet
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
          CALL abort(__STAMP__,&
              'Non-Critical error in domain reconstrution. IF NOT encountered, something is terrible wrong.')
        END IF
        !BC(1:4,newSideID)=RecvMsg%BC(1:4,haloSideID)
        BC(newSideID)=RecvMsg%BC(haloSideID)
        SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)
      END IF ! isDoubleSide
      ! copy Bezier to new side id
      IF(.NOT.isDoubleSide) THEN
        BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newBCSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
        ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
        SlabNormals(1:3,1:3,newBCSideID)=RecvMsg%SlabNormals(1:3,1:3,haloSideID)
        SlabIntervalls(1:6,newBCSideID) =RecvMsg%SlabIntervalls(1:6 ,haloSideID) 
        BoundingBoxIsEmpty(newBCSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
      ELSE
        IF(RecvMsg%ElemToSide(2,ilocSide,iElem).EQ.0)THEN
          SlabNormals(1:3,1:3,newSideID)=RecvMsg%SlabNormals(1:3,1:3,haloSideID)
          SlabIntervalls(1:6,newSideID) =RecvMsg%SlabIntervalls(1:6 ,haloSideID) 
          BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( haloSideID) 
        END IF
      END IF
      ! build entry to PartElemToSide
      PartBCSideList(newSideID)=newBCSideID !tmpBCSides+haloinc(haloSideID)
      PartElemToSide(E2S_SIDE_ID,iLocSide,newElemId)=newSideID
      PartElemToSide(E2S_FLIP,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
    END DO ! ilocSide
    ! set native elemID
    PartHaloToProc(1,newElemId)=RecvMsg%NativeElemID(iElem)
    PartHaloToProc(2,newElemId)=iProc
    XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
    dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
  END DO ! iElem
  ! build rest: PartNeighborElemID, PartLocSideID
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
        PartneighborlocSideID(ilocSide,iElem)=PartSideToElem(S2E_NB_LOC_SIDE_ID,SideID)
        PartneighborElemID   (ilocSide,iElem)=PartSideToElem(S2E_NB_ELEM_ID,SideID)
      ELSE
        ! SideID of master
        PartneighborlocSideID(ilocSide,iElem)=PartSideToElem(S2E_LOC_SIDE_ID,SideID)
        PartneighborElemID   (ilocSide,iElem)=PartSideToElem(S2E_ELEM_ID,SideID)
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


SUBROUTINE WriteParticlePartitionInformation(nPlanar,nBilinear,nCurved)
!===================================================================================================================================
! write the particle partition information to file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,            ONLY:nSides,nElems
USE MOD_Particle_MPI_Vars,    ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,   ONLY:nTotalSides,nTotalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: nPlanar,nBilinear,nCurved
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=10)          :: formatstr
INTEGER,ALLOCATABLE        :: nNBProcs_glob(:), ProcInfo_glob(:,:),NBInfo_glob(:,:), NBInfo(:), tmparray(:,:)
REAL,ALLOCATABLE           :: tmpreal(:,:)
INTEGER                    :: ProcInfo(7),nNBmax,iProc,i,j,ioUnit
!===================================================================================================================================


!output partitioning info
ProcInfo(1)=nElems
ProcInfo(2)=nSides
ProcInfo(3)=nTotalElems-nElems
ProcInfo(4)=nTotalSides-nSides
ProcInfo(5)=nPlanar
ProcInfo(6)=nBilinear
ProcInfo(7)=nCurved
IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:PartMPI%nProcs-1))
  ALLOCATE(ProcInfo_glob(7,0:PartMPI%nProcs-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot 
CALL MPI_GATHER(PartMPI%nMPINeighbors,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,PartMPI%COMM,iError)
CALL MPI_GATHER(ProcInfo,7,MPI_INTEGER,ProcInfo_glob,7,MPI_INTEGER,0,PartMPI%COMM,iError)
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
  OPEN(UNIT=ioUnit,FILE='particlepartitionInfo.out',STATUS='REPLACE')
  WRITE(ioUnit,*)'Particle Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',PartMPI%nProcs
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  WRITE(ioUnit,'(9(A15))')'Rank','nElems','nSides','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved','nNBProcs'
  WRITE(ioUnit,'(A90,A45)')&
      '==========================================================================================',&
      '============================================='
  !statistics
  ALLOCATE(tmparray(8,0:3),tmpreal(8,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=0      !max
  tmparray(:,3)=HUGE(1)   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i) ! nElems
    tmparray(2,0)=Procinfo_glob(2,i) ! nSides
    tmparray(3,0)=Procinfo_glob(3,i) ! nHaloElems
    tmparray(4,0)=Procinfo_glob(4,i) ! nHaloSides
    tmparray(5,0)=Procinfo_glob(5,i) ! nHaloSides
    tmparray(6,0)=Procinfo_glob(6,i) ! nHaloSides
    tmparray(7,0)=Procinfo_glob(7,i) ! nHaloSides
    tmparray(8,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,8
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
    tmparray(8,0)=nNBProcs_glob(i)
    DO j=1,5
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2 
    END DO ! j
  END DO ! i
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(PartMPI%nProcs))
  WRITE(ioUnit,'(A15,8(5X,F10.2))')'   MEAN        ',tmpreal(:,1)
  WRITE(ioUnit,'(A90,A45)')&
      '------------------------------------------------------------------------------------------',&
      '---------------------------------------------'
  WRITE(ioUnit,'(A15,8(5X,F10.2))')'   RMS         ',tmpreal(:,2)
  WRITE(ioUnit,'(A90,A45)')&
      '------------------------------------------------------------------------------------------',&
      '---------------------------------------------'
  WRITE(ioUnit,'(A15,8(5X,I10))')'   MIN         ',tmparray(:,3)
  WRITE(ioUnit,'(A90,A45)')&
      '------------------------------------------------------------------------------------------',&
      '---------------------------------------------'
  WRITE(ioUnit,'(A15,8(5X,I10))')'   MAX         ',tmparray(:,2)
  WRITE(ioUnit,'(A90,A45)')&
      '------------------------------------------------------------------------------------------',&
      '---------------------------------------------'
  WRITE(ioUnit,'(A90,A45)')&
      '==========================================================================================',&
      '============================================='
  DO i=0,PartMPI%nProcs-1
    WRITE(ioUnit,'(9(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
    WRITE(ioUnit,'(A90,A45)')&
      '------------------------------------------------------------------------------------------',&
      '---------------------------------------------'
  END DO! i
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
CHARACTER(LEN=10)          :: formatstr
INTEGER,ALLOCATABLE        :: nNBProcs_glob(:), ProcInfo_glob(:,:),NBInfo_glob(:,:), NBInfo(:), tmparray(:,:)
REAL,ALLOCATABLE           :: tmpreal(:,:)
INTEGER                    :: ProcInfo(8),nNBmax,iProc,i,j,ioUnit
!===================================================================================================================================


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
  OPEN(UNIT=ioUnit,FILE='particlepartitionInfo.out',STATUS='REPLACE')
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
