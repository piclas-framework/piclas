!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Photon_Tracking
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE PhotonTriaTracking
  MODULE PROCEDURE PhotonTriaTracking
END INTERFACE

PUBLIC::PhotonTriaTracking, Photon2DSymTracking
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE PhotonTriaTracking()
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
! 1) Loop over all particles that are still inside
!    2) Perform tracking until the particle is considered "done" (either localized or deleted)
!       2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
!           triangle (ParticleInsideQuad3D)
!       2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
!           each side of the element (ParticleThroughSideCheck3DFast)
!       2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
!           treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
!           the determinants
!    3) In case of a boundary, determine the intersection and perform the appropriate boundary interaction (GetBoundaryInteraction)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars,      ONLY:PartBound
USE MOD_RadiationTrans_Vars,         ONLY:PhotonProps,RadObservation_Emission, CalcRadObservationPoint, RadObservation_EmissionPart
USE MOD_Photon_TrackingTools,        ONLY:PhotonThroughSideCheck3DFast, PhotonIntersectionWithSide,CalcAbsoprtion
USE MOD_Photon_TrackingTools,        ONLY:PerfectPhotonReflection, DiffusePhotonReflection, CalcWallAbsoprtion, PointInObsCone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NblocSideID, NbElemID, ind, nbSideID, nMortarElems, BCType, localSideID
INTEGER                          :: ElemID,OldElemID,nlocSides
INTEGER                          :: LocalSide
INTEGER                          :: NrOfThroughSides, ind2
INTEGER                          :: SideID,TempSideID,iLocSide
INTEGER                          :: TriNum, LocSidesTemp(1:6),TriNumTemp(1:6), GlobSideTemp(1:6)
INTEGER                          :: SecondNrOfThroughSides, indSide
INTEGER                          :: DoneLastElem(1:4,1:6) ! 1:3: 1=Element,2=LocalSide,3=TriNum 1:2: 1=last 2=beforelast
LOGICAL                          :: ThroughSide, Done
LOGICAL                          :: oldElemIsMortar, isMortarSideTemp(1:6), doCheckSide
REAL                             :: minRatio, intersecDist, intersecDistVec(3)
REAL                             :: IntersectionPos(1:3), IntersectionPosTemp(1:3)

!===================================================================================================================================
Done = .FALSE.
ElemID = PhotonProps%ElemID
SideID = 0
DoneLastElem(:,:) = 0
! 1) Loop tracking until Photon is considered "done" (either absorbed or deleted)
DO WHILE (.NOT.Done)
  oldElemIsMortar = .FALSE.
  NrOfThroughSides = 0
  LocSidesTemp(:) = 0
  TriNumTemp(:) = 0
  GlobSideTemp = 0
  isMortarSideTemp = .FALSE.
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  DO iLocSide=1,nlocSides
    TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
    ! Side is not one of the 6 local sides
    IF (localSideID.LE.0) CYCLE
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)
    IF (NbElemID.LT.0) THEN ! Mortar side
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.-1)    
      DO ind = 1, nMortarElems
        nbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
        NbElemID = SideInfo_Shared(SIDE_NBELEMID,nbSideID)
        ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
        ! no way to recover it during runtime
        IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)
        ! For small mortar sides, SIDE_NBSIDEID contains the SideID of the corresponding big mortar side
        nbSideID = SideInfo_Shared(SIDE_NBSIDEID,nbSideID)
        NblocSideID =  SideInfo_Shared(SIDE_LOCALID,nbSideID)
        DO TriNum = 1,2
          ThroughSide = .FALSE.
          CALL PhotonThroughSideCheck3DFast(NblocSideID,NbElemID,ThroughSide,TriNum, .TRUE.)
          IF (ThroughSide) THEN
            ! Store the information for this side for future checks, if this side was already treated
            oldElemIsMortar = .TRUE.
            NrOfThroughSides = NrOfThroughSides + 1
            LocSidesTemp(NrOfThroughSides) = NblocSideID
            TriNumTemp(NrOfThroughSides) = TriNum
            GlobSideTemp(NrOfThroughSides) = nbSideID
            isMortarSideTemp(NrOfThroughSides)  = .TRUE.
            SideID = nbSideID
            LocalSide = NblocSideID
          END IF
        END DO
      END DO
    ELSE  ! Regular side
      DO TriNum = 1,2
        ThroughSide = .FALSE.
        CALL PhotonThroughSideCheck3DFast(localSideID,ElemID,ThroughSide,TriNum)
        IF (ThroughSide) THEN
          NrOfThroughSides = NrOfThroughSides + 1
          LocSidesTemp(NrOfThroughSides) = localSideID
          TriNumTemp(NrOfThroughSides) = TriNum
          GlobSideTemp(NrOfThroughSides) = TempSideID
          SideID = TempSideID
          LocalSide = localSideID
        END IF
      END DO
    END IF  ! Mortar or regular side
  END DO  ! iLocSide=1,6
  TriNum = TriNumTemp(1)
  ! ----------------------------------------------------------------------------
  ! Addition treatment if particle did not cross any sides or it crossed multiple sides
  IF (NrOfThroughSides.NE.1) THEN
    ! 2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
    ! treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
    ! the determinants
    IF (NrOfThroughSides.EQ.0) THEN
      ! Particle appears to have not crossed any of the checked sides. Deleted!
      IPWRITE(*,*) 'Error in Photon TriaTracking! Photon lost. Element:', ElemID
      IPWRITE(*,*) 'LastPos: ', PhotonProps%PhotonLastPos(1:3)
      IPWRITE(*,*) 'Pos:     ', PhotonProps%PhotonPos(1:3)
      IPWRITE(*,*) 'Direction:', PhotonProps%PhotonDirection(1:3)
      IPWRITE(*,*) 'Photon deleted!'
      Done = .TRUE.
      EXIT
    ELSE IF (NrOfThroughSides.GT.1) THEN
      ! Use the slower search method if particle appears to have crossed more than one side (possible for irregular hexagons
      ! and in the case of mortar elements)
      SecondNrOfThroughSides = 0
      minRatio = 1E90
      oldElemIsMortar = .FALSE.
      DO ind2 = 1, NrOfThroughSides
        doCheckSide = .TRUE.
        ! Check if this side was already treated
        DO indSide = 2, 6
          IF((DoneLastElem(1,indSide).EQ.ElemID).AND. &
             (DoneLastElem(4,indSide).EQ.GlobSideTemp(ind2)).AND. &
             (DoneLastElem(3,indSide).EQ.TriNumTemp(ind2))) THEN
            doCheckSide = .FALSE.
          END IF
        END DO
        IF (doCheckSide) THEN
          IF (isMortarSideTemp(ind2)) THEN  ! Mortar side
            NbElemID = SideInfo_Shared(SIDE_ELEMID,GlobSideTemp(ind2))
            ! Get the determinant between the old and new particle position and the nodes of the triangle which was crossed
            CALL PhotonIntersectionWithSide(LocSidesTemp(ind2),NbElemID,TriNumTemp(ind2), IntersectionPosTemp, .TRUE.)
            intersecDistVec(1:3) = IntersectionPosTemp(1:3) - PhotonProps%PhotonLastPos(1:3)
            intersecDist = DOT_PRODUCT(intersecDistVec, intersecDistVec)
            ! If the particle is inside the neighboring mortar element, it moved through this side            
            ! Ratio is always negative since detM(=detLastPartPos) is negative or zero, i.e. maximum abs is wanted
            ! The closer the intersected side is to the last particle position the greater the absolute ratio will be
            IF (intersecDist.LT.minRatio) THEN
              IntersectionPos = IntersectionPosTemp
              minRatio = intersecDist
              SecondNrOfThroughSides = SecondNrOfThroughSides + 1
              SideID = GlobSideTemp(ind2)
              LocalSide = LocSidesTemp(ind2)
              TriNum = TriNumTemp(ind2)
              oldElemIsMortar = .TRUE.
            END IF
          ELSE  ! Regular side
            CALL PhotonIntersectionWithSide(LocSidesTemp(ind2),NbElemID,TriNumTemp(ind2), IntersectionPosTemp)
            intersecDistVec(1:3) = IntersectionPosTemp(1:3) - PhotonProps%PhotonLastPos(1:3)
            intersecDist = DOT_PRODUCT(intersecDistVec, intersecDistVec)
            IF (intersecDist.LT.minRatio) THEN
              IntersectionPos = IntersectionPosTemp
              minRatio = intersecDist
              SecondNrOfThroughSides = SecondNrOfThroughSides + 1
              SideID = GlobSideTemp(ind2)
              LocalSide = LocSidesTemp(ind2)
              TriNum = TriNumTemp(ind2)
              oldElemIsMortar = .FALSE.
            END IF
          END IF  ! isMortarSideTemp = T/F
        END IF  ! doCheckSide
      END DO  ! ind2 = 1, NrOfThroughSides
      ! Particle that went through multiple sides first, but did not cross any sides during the second check -> Deleted!
      IF (SecondNrOfThroughSides.EQ.0) THEN
        IPWRITE(*,*) 'Error in Photon TriaTracking! Photon lost on second check. Element:', ElemID
        IPWRITE(*,*) 'LastPos: ', PhotonProps%PhotonLastPos(1:3)
        IPWRITE(*,*) 'Pos:     ', PhotonProps%PhotonPos(1:3)
        IPWRITE(*,*) 'Direction:', PhotonProps%PhotonDirection(1:3)
        IPWRITE(*,*) 'Photon deleted!'
        Done = .TRUE.
        EXIT
      END IF
    END IF  ! NrOfThroughSides.EQ.0/.GT.1
  END IF  ! NrOfThroughSides.NE.1
  ! ----------------------------------------------------------------------------
  ! 3) In case of a boundary, perform the appropriate boundary interaction
  IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
    OldElemID=ElemID
    BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
    SELECT CASE(BCType)
    CASE(1) !PartBound%OpenBC)
      IF (NrOfThroughSides.LT.2) CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos)
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      IF (CalcRadObservationPoint) THEN
        IF (PointInObsCone(IntersectionPos(1:3))) THEN
          RadObservation_Emission(PhotonProps%WaveLength) = RadObservation_Emission(PhotonProps%WaveLength) + PhotonProps%PhotonEnergy
          RadObservation_EmissionPart(PhotonProps%WaveLength) = RadObservation_EmissionPart(PhotonProps%WaveLength) + 1
        END IF
      END IF
      DONE = .TRUE.
    CASE(2)
      IF (PartBound%PhotonSpecularReflection(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) THEN
        IF (NrOfThroughSides.LT.2) THEN
          CALL PerfectPhotonReflection(LocalSide,ElemID,TriNum, IntersectionPos, .FALSE.)
        ELSE 
          CALL PerfectPhotonReflection(LocalSide,ElemID,TriNum, IntersectionPos, .TRUE.)
        END IF
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      ELSE
        IF (NrOfThroughSides.LT.2) THEN
          CALL DiffusePhotonReflection(LocalSide,ElemID,TriNum, IntersectionPos, .FALSE.)
        ELSE 
          CALL DiffusePhotonReflection(LocalSide,ElemID,TriNum, IntersectionPos, .TRUE.)
        END IF      
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
        IF (.NOT.DONE) CALL CalcWallAbsoprtion(SideID, DONE)
      END IF
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,' ERROR: PartBound not associated!. (unknown case)',BCType,999.)
    END SELECT !PartBound%MapToPartBC(BC(SideID)


    IF ((BCType.EQ.2).OR.(BCType.EQ.10)) THEN
      DoneLastElem(:,:) = 0
    ELSE
      DO ind2= 5, 1, -1
        DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
      END DO
      DoneLastElem(1,1) = OldElemID
      DoneLastElem(2,1) = LocalSide
      DoneLastElem(3,1) = TriNum
      DoneLastElem(4,1) = SideID
    END IF
  ELSE  ! BC(SideID).LE.0
    DO ind2= 5, 1, -1
      DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
    END DO
    DoneLastElem(1,1) = ElemID
    DoneLastElem(2,1) = LocalSide
    DoneLastElem(3,1) = TriNum
    DoneLastElem(4,1) = SideID
    IF (oldElemIsMortar) THEN
      ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos,.TRUE.)
      CALL CalcAbsoprtion(IntersectionPos(1:3),DoneLastElem(1,1), DONE)
    ELSE
      CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos)
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
    END IF
  END IF  ! BC(SideID).GT./.LE. 0
  IF (ElemID.LT.1) THEN
    CALL abort(&
     __STAMP__ &
     ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
  END IF
END DO  ! .NOT.PartisDone


END SUBROUTINE PhotonTriaTracking



SUBROUTINE Photon2DSymTracking()
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
! 1) Loop over all particles that are still inside
!    2) Perform tracking until the particle is considered "done" (either localized or deleted)
!       2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
!           triangle (ParticleInsideQuad3D)
!       2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
!           each side of the element (ParticleThroughSideCheck3DFast)
!       2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
!           treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
!           the determinants
!    3) In case of a boundary, determine the intersection and perform the appropriate boundary interaction (GetBoundaryInteraction)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars,      ONLY:PartBound
USE MOD_RadiationTrans_Vars,         ONLY:PhotonProps, RadObservation_Emission, CalcRadObservationPoint,RadObservation_EmissionPart
USE MOD_Photon_TrackingTools,        ONLY:CalcAbsoprtion, CalcWallAbsoprtion, DiffusePhotonReflection2D, PointInObsCone
USE MOD_Photon_TrackingTools,        ONLY:PhotonIntersectionWithSide2D, RotatePhotonIn2DPlane, PerfectPhotonReflection2D
USE MOD_Photon_TrackingTools,        ONLY:PhotonIntersectSensor
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NbElemID, ind, nbSideID, nMortarElems, BCType, NblocSideID, OrigElem
INTEGER                          :: ElemID, OldElemID, LocalSide, NrOfThroughSides, localSideID, nlocSides
INTEGER                          :: SideID, TempSideID, iLocSide, correctSide, LastSide
INTEGER                          :: LocSidesTemp(1:6), GlobSideTemp(1:6)
LOGICAL                          :: oldElemIsMortar, isMortarSideTemp(1:6), isLastSide, ThroughSide, Done
REAL                             :: IntersectionPos(1:6), IntersectionPosTemp(1:6,3), DistanceTemp(1:6), Distance
!===================================================================================================================================
Done = .FALSE.
ElemID = PhotonProps%ElemID
OrigElem = ElemID
SideID = 0
LastSide = 0
! 1) Loop tracking until Photon is considered "done" (either absorbed or deleted)
DO WHILE (.NOT.Done)
  oldElemIsMortar = .FALSE.
  NrOfThroughSides = 0
  LocSidesTemp = 0
  GlobSideTemp = 0
  DistanceTemp = 0.0
  IntersectionPosTemp = 0.0
  isMortarSideTemp = .FALSE.
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  DO iLocSide=1,nlocSides
    isLastSide = .FALSE. 
    TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
    ! Side is not one of the 6 local sides
    IF (localSideID.LE.0) CYCLE
    IF (SideIsSymSide(TempSideID)) CYCLE
    IF (SideInfo_Shared(SIDE_BCID,TempSideID).GT.0) THEN
      IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,TempSideID))).EQ.11) CYCLE
    END IF
    IF (LastSide.EQ.TempSideID) isLastSide = .TRUE.
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)
    IF (NbElemID.LT.0) THEN ! Mortar side
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.-1)   
      DO ind = 1, nMortarElems
        isLastSide = .FALSE. 
        nbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
        NbElemID = SideInfo_Shared(SIDE_NBELEMID,nbSideID)
        ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
        ! no way to recover it during runtime
        IF (LastSide.EQ.nbSideID) isLastSide = .TRUE.
        IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)
        ! For small mortar sides, SIDE_NBSIDEID contains the SideID of the corresponding big mortar side
        nbSideID = SideInfo_Shared(SIDE_NBSIDEID,nbSideID)
        NblocSideID =  SideInfo_Shared(SIDE_LOCALID,nbSideID)
        ThroughSide = .FALSE.
        CALL PhotonIntersectionWithSide2D(NblocSideID,NbElemID,ThroughSide,IntersectionPos, isLastSide,Distance)          
        IF (ThroughSide) THEN               
          ! Store the information for this side for future checks, if this side was already treated
          oldElemIsMortar = .TRUE.
          NrOfThroughSides = NrOfThroughSides + 1
          LocSidesTemp(NrOfThroughSides) = NblocSideID
          GlobSideTemp(NrOfThroughSides) = nbSideID
          DistanceTemp(NrOfThroughSides) = Distance
          IntersectionPosTemp(1:3,NrOfThroughSides) = IntersectionPos(1:3)
          isMortarSideTemp(NrOfThroughSides)  = .TRUE.
          SideID = nbSideID
          LocalSide = NblocSideID
        END IF
      END DO
    ELSE  ! Regular side
      ThroughSide = .FALSE.
      CALL PhotonIntersectionWithSide2D(localSideID,ElemID,ThroughSide,IntersectionPos, isLastSide, Distance)
      IF (ThroughSide) THEN
        NrOfThroughSides = NrOfThroughSides + 1
        SideID = TempSideID
        LocalSide = localSideID
        LocSidesTemp(NrOfThroughSides) = localSideID
        GlobSideTemp(NrOfThroughSides) = TempSideID
        DistanceTemp(NrOfThroughSides) = Distance
        IntersectionPosTemp(1:3,NrOfThroughSides) = IntersectionPos(1:3)
      END IF
    END IF  ! Mortar or regular side
  END DO  ! iLocSide=1,4
  ! ----------------------------------------------------------------------------
  ! Addition treatment if particle did not cross any sides or it crossed multiple sides
  IF (NrOfThroughSides.NE.1) THEN
    ! 2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
    ! treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
    ! the determinants
    IF (NrOfThroughSides.EQ.0) THEN
      ! Particle appears to have not crossed any of the checked sides. Deleted!
      IPWRITE(*,*) 'Error in Photon 2DAxisTracking! Photon lost. Element:', ElemID
      IPWRITE(*,*) 'LastPos: ', PhotonProps%PhotonLastPos(1:3)
      IPWRITE(*,*) 'Pos:     ', PhotonProps%PhotonPos(1:3)
      IPWRITE(*,*) 'Direction:', PhotonProps%PhotonDirection(1:3)
      IPWRITE(*,*) 'OrigElem:', OrigElem
      IPWRITE(*,*) 'Photon deleted!'
      Done = .TRUE.
      EXIT
    ELSE IF (NrOfThroughSides.GT.1) THEN
      correctSide = MINLOC(DistanceTemp(1:NrOfThroughSides),DIM=1)
      oldElemIsMortar = isMortarSideTemp(correctSide)
      SideID = GlobSideTemp(correctSide)
      LocalSide = LocSidesTemp(correctSide)
      IntersectionPos(1:3) = IntersectionPosTemp(1:3,correctSide)
    END IF
  END IF  ! NrOfThroughSides.NE.1
  ! ----------------------------------------------------------------------------
  ! 3) In case of a boundary, perform the appropriate boundary interaction
  IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
    OldElemID=ElemID
    BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
    SELECT CASE(BCType)
    CASE(1) !PartBound%OpenBC)
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      DONE = .TRUE.
      IF (CalcRadObservationPoint) THEN        
        IF (PointInObsCone(IntersectionPos(1:3))) THEN
          IF (PhotonIntersectSensor(IntersectionPos(1:3), PhotonProps%PhotonDirection(1:3))) THEN
            RadObservation_Emission(PhotonProps%WaveLength) = RadObservation_Emission(PhotonProps%WaveLength) + PhotonProps%PhotonEnergy
            RadObservation_EmissionPart(PhotonProps%WaveLength) = RadObservation_EmissionPart(PhotonProps%WaveLength) + 1
          END IF
        END IF
      END IF
      CYCLE
    CASE(2)
      IF (PartBound%PhotonSpecularReflection(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) THEN
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
        IF (.NOT.DONE) CALL PerfectPhotonReflection2D(LocalSide,ElemID, IntersectionPos) 
      ELSE
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)      
        IF (.NOT.DONE) CALL CalcWallAbsoprtion(SideID, DONE)
        IF (.NOT.DONE) CALL DiffusePhotonReflection2D(LocalSide,ElemID, IntersectionPos) 
      END IF
      LastSide = SideID
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,' ERROR: PartBound not associated!. (unknown case)',999,999.)
    END SELECT !PartBound%MapToPartBC(BC(SideID) 
  ELSE  ! BC(SideID).LE.0
    IF (oldElemIsMortar) THEN
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      LastSide = SideID
    ELSE
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
      LastSide = SideInfo_Shared(SIDE_NBSIDEID,SideID)
    END IF
  END IF  ! BC(SideID).GT./.LE. 0
  
  IF (.NOT.DONE) CALL RotatePhotonIn2DPlane(IntersectionPos(1:3))
  IF (ElemID.LT.1) THEN
    CALL abort(&
     __STAMP__ &
     ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
  END IF
END DO  ! .NOT.PartisDone
END SUBROUTINE Photon2DSymTracking

END MODULE MOD_Photon_Tracking
