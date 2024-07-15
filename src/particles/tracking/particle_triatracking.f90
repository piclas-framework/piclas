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

MODULE MOD_Particle_TriaTracking
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

! Define an interface for the function pointer
ABSTRACT INTERFACE
  SUBROUTINE SingleParticleTriaTrackingInterface(i,IsInterPlanePart)
    INTEGER,INTENT(IN)                :: i
    LOGICAL,INTENT(IN),OPTIONAL       :: IsInterPlanePart
  END SUBROUTINE
END INTERFACE

!> Pointer defining the particle tracking routine based on the symmetry order
PROCEDURE(SingleParticleTriaTrackingInterface),POINTER :: SingleParticleTriaTracking => NULL()
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC::InitSingleParticleTriaTracking,ParticleTriaTracking,SingleParticleTriaTracking
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

!==================================================================================================================================!
!> Initialize SingleParticleTriaTracking depending on symmetry dimension using a function pointer
!==================================================================================================================================!
SUBROUTINE InitSingleParticleTriaTracking()
! MODULES
USE MOD_Globals
USE MOD_Symmetry_Vars            ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================

IF (Symmetry%Order.EQ.3) THEN
  SingleParticleTriaTracking => SingleParticleTriaTracking3D
ELSE IF (Symmetry%Order.EQ.2.OR.Symmetry%Order.EQ.1) THEN
  SingleParticleTriaTracking => SingleParticleTriaTracking1D2D
ELSE
  CALL abort(__STAMP__,'ERROR in InitSingleParticleTriaTracking: Function pointer could not be properly defined!')
END IF

END SUBROUTINE InitSingleParticleTriaTracking


SUBROUTINE ParticleTriaTracking()
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
! 1) Loop over all particles that are still inside and call SingleParticleTriaTracking
! 2) Loop again over all particles that are created in case of using rotational periodic inter planes
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY: PEM,PDM,InterPlanePartNumber, InterPlanePartIndx, UseRotRefSubCycling,nSubCyclingSteps
USE MOD_Particle_Vars               ,ONLY: NewPosSubCycling, GlobalElemIDSubCycling, LastPartPosSubCycling
USE MOD_Particle_Vars               ,ONLY: InRotRefFrame,InRotRefFrameSubCycling, PartVeloRotRefSubCycling, LastVeloRotRefSubCycling
USE MOD_DSMC_Vars                   ,ONLY: RadialWeighting
USE MOD_DSMC_Symmetry               ,ONLY: DSMC_2D_RadialWeighting, DSMC_2D_SetInClones
USE MOD_part_tools                  ,ONLY: ParticleOnProc
!----- Used for RotRef Subcycling
USE MOD_part_RHS                    ,ONLY: CalcPartPosInRotRef
USE MOD_TimeDisc_Vars               ,ONLY: dt
USE MOD_Particle_Vars               ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep, Species, PartState, LastPartPos, PartSpecies
USE MOD_Particle_Vars               ,ONLY: PartVeloRotRef, LastPartVeloRotRef
!-----
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i, InterPartID, iStep
REAL                             :: dtVar
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

IF(RadialWeighting%PerformCloning) CALL DSMC_2D_SetInClones()
InterPlanePartNumber = 0
! 1) Loop over all particles that are still inside
DO i = 1,PDM%ParticleVecLength
  ! 1.) PDM%ParticleInside(i) indicates that the particle is active and not deleted
  ! 2.) PEM%LastGlobalElemID(i) is not set in ParticleInserting() to the actual element index but initialized with zero.
  !     ParticleInserting() has been moved in front of PerformTracking() for some time disc methods, hence, skip these
  !     particles here as they have not yet been pushed. GlobalElemID is copied to LastGlobalElemID at the beginning of the time
  !     stepping method. Note that PEM%LastGlobalElemID(i)<0 is used for Symmetry%Order.LE.2 where the information
  !     PEM%LastGlobalElemID(i) = -SideID is stored (surface flux).
  IF (PDM%ParticleInside(i).AND.PEM%LastGlobalElemID(i).NE.0) THEN
    IF(UseRotRefSubCycling) THEN
      CALL SingleParticleTriaTrackingSubcycling(iPart=i)
    ELSE
      CALL SingleParticleTriaTracking(i=i)
    END IF
  END IF
  ! Particle treatment for an axisymmetric simulation (cloning/deleting particles)
  IF(RadialWeighting%PerformCloning) THEN
    IF(PDM%ParticleInside(i).AND.(ParticleOnProc(i))) THEN
      CALL DSMC_2D_RadialWeighting(i,PEM%GlobalElemID(i))
    END IF
  END IF
END DO ! i = 1,PDM%ParticleVecLength
! 2) Loop again over all inter plane particles
IF(InterPlanePartNumber.GT.0) THEN
  DO i = 1,InterPlanePartNumber
    InterPartID = InterPlanePartIndx(i)
    PDM%ParticleInside(InterPartID) = .TRUE.
    IF(UseRotRefSubCycling) THEN
      CALL SingleParticleTriaTrackingSubcycling(iPart=InterPartID,IsInterPlanePart=.TRUE.)
    ELSE
      CALL SingleParticleTriaTracking(i=InterPartID,IsInterPlanePart=.TRUE.)
    END IF
    ! Particle treatment for an axisymmetric simulation (cloning/deleting particles)
    IF(RadialWeighting%PerformCloning) THEN
      IF(PDM%ParticleInside(InterPartID).AND.(ParticleOnProc(InterPartID))) THEN
        CALL DSMC_2D_RadialWeighting(InterPartID,PEM%GlobalElemID(InterPartID))
      END IF
    END IF
  END DO ! i = 1,InterPlanePartNumber
END IF

END SUBROUTINE ParticleTriaTracking


SUBROUTINE SingleParticleTriaTracking3D(i,IsInterPlanePart)
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
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
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID
USE MOD_Particle_Vars               ,ONLY: PEM,PDM,PartSpecies
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Symmetry_Vars               ,ONLY: Symmetry
USE MOD_Particle_Vars               ,ONLY: UseRotRefFrame, InRotRefFrame, RotRefFrameOmega, PartVeloRotRef
USE MOD_Particle_Mesh_Tools         ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Tracking_vars      ,ONLY: ntracks,MeasureTrackTime,CountNbrOfLostParts, NbrOfLostParticles, TrackInfo
USE MOD_Particle_Tracking_vars      ,ONLY: DisplayLostParticles
USE MOD_Part_Tools                  ,ONLY: StoreLostParticleProperties
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound
USE MOD_Particle_Intersection       ,ONLY: ParticleThroughSideCheck3DFast, ParticleThroughSideLastPosCheck, IntersectionWithWall
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_part_tools                  ,ONLY: ParticleOnProc, InRotRefFrameCheck
USE MOD_part_operations             ,ONLY: RemoveParticle
#if USE_LOADBALANCE
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime, LBElemSplitTime, LBElemPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: i
LOGICAL,INTENT(IN),OPTIONAL       :: IsInterPlanePart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NblocSideID, NbElemID, CNElemID, ind, nbSideID, nMortarElems,BCType
INTEGER                          :: ElemID,flip,OldElemID,nlocSides
INTEGER                          :: LocalSide
INTEGER                          :: NrOfThroughSides, ind2
INTEGER                          :: SideID,TempSideID,iLocSide, localSideID
INTEGER                          :: TriNum, LocSidesTemp(1:6),TriNumTemp(1:6), GlobSideTemp(1:6)
INTEGER                          :: SecondNrOfThroughSides, indSide
INTEGER                          :: DoneLastElem(1:4,1:6) ! 1:3: 1=Element,2=LocalSide,3=TriNum 1:2: 1=last 2=beforelast
LOGICAL                          :: ThroughSide, InElementCheck,PartisDone
LOGICAL                          :: crossedBC, oldElemIsMortar, isMortarSideTemp(1:6), doCheckSide
REAL                             :: det(6,2),detM,ratio,minRatio, detPartPos
REAL, PARAMETER                  :: eps = 0
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
IF (MeasureTrackTime) nTracks=nTracks+1
PartisDone = .FALSE.
ElemID = PEM%LastGlobalElemID(i)
TrackInfo%CurrElem=ElemID
SideID = 0
DoneLastElem(:,:) = 0
! 2) Loop tracking until particle is considered "done" (either localized or deleted)
DO WHILE (.NOT.PartisDone)
  ! 2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
  !     triangle (ParticleInsideQuad3D)
  oldElemIsMortar = .FALSE.
  CALL ParticleInsideQuad3D(PartState(1:3,i),ElemID,InElementCheck,det)
  IF (InElementCheck) THEN
    ! If particle is inside the given ElemID, set new PEM%GlobalElemID and stop tracking this particle ->PartisDone
    PEM%GlobalElemID(i) = ElemID
    PartisDone = .TRUE.
  ELSE
    ! 2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
    !     each side of the element (ParticleThroughSideCheck3DFast)
    NrOfThroughSides = 0
    LocSidesTemp(:) = 0
    TriNumTemp(:) = 0
    GlobSideTemp = 0
    isMortarSideTemp = .FALSE.
    TrackInfo%PartTrajectory(1:3)=PartState(1:3,i) - LastPartPos(1:3,i)
    TrackInfo%lengthPartTrajectory=VECNORM(TrackInfo%PartTrajectory(1:3))
    IF(ABS(TrackInfo%lengthPartTrajectory).GT.0.) TrackInfo%PartTrajectory = TrackInfo%PartTrajectory &
                                                                              / TrackInfo%lengthPartTrajectory
    nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
    DO iLocSide=1,nlocSides
      TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
      ! Skip symmetry side
      IF(Symmetry%Order.EQ.2) THEN
        IF(SideIsSymSide(TempSideID)) CYCLE
      END IF
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
            CALL ParticleThroughSideCheck3DFast(i,NblocSideID,NbElemID,ThroughSide,TriNum,.TRUE.)
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
          IF (det(localSideID,TriNum).LE.-eps) THEN
            ThroughSide = .FALSE.
            CALL ParticleThroughSideCheck3DFast(i,localSideID,ElemID,ThroughSide,TriNum)
            IF (ThroughSide) THEN
              NrOfThroughSides = NrOfThroughSides + 1
              LocSidesTemp(NrOfThroughSides) = localSideID
              TriNumTemp(NrOfThroughSides) = TriNum
              GlobSideTemp(NrOfThroughSides) = TempSideID
              SideID = TempSideID
              LocalSide = localSideID
            END IF
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
        IF(DisplayLostParticles)THEN
          IPWRITE(*,*) 'Error in Particle TriaTracking! Particle Number',i,'lost. Element:', ElemID,'(species:',PartSpecies(i),')'
          IPWRITE(*,*) 'LastPos: ', LastPartPos(1:3,i)
          IPWRITE(*,*) 'Pos:     ', PartState(1:3,i)
          IPWRITE(*,*) 'Velo:    ', PartState(4:6,i)
          IPWRITE(*,*) 'Particle deleted!'
        END IF ! DisplayLostParticles
        IF(CountNbrOfLostParts) THEN
          CALL StoreLostParticleProperties(i, ElemID)
          NbrOfLostParticles=NbrOfLostParticles+1
        END IF
        CALL RemoveParticle(i)
        PartisDone = .TRUE.
        EXIT
      ELSE IF (NrOfThroughSides.GT.1) THEN
        ! Use the slower search method if particle appears to have crossed more than one side (possible for irregular hexagons
        ! and in the case of mortar elements)
        SecondNrOfThroughSides = 0
        minRatio = 0
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
              ! Get the element number of the smaller neighboring element
              NbElemID = SideInfo_Shared(SIDE_ELEMID,GlobSideTemp(ind2))
              ! Get the determinant between the old and new particle position and the nodes of the triangle which was crossed
              CALL ParticleThroughSideLastPosCheck(i,LocSidesTemp(ind2),NbElemID,InElementCheck,TriNumTemp(ind2),detM, &
                                                    isMortarSide=.TRUE.,detPartPos=detPartPos)
              ! If the particle is inside the neighboring mortar element, it moved through this side
              IF (InElementCheck) THEN
                IF((detM.EQ.0).AND.(detPartPos.EQ.0)) CYCLE ! Particle moves within side
                ! Determining which side was crossed first by comparing the ratio of the spatial product of PartPos->Tri-Nodes
                ! and LastPartPos->Tri-Nodes
                IF((detM.EQ.0).AND.(minRatio.EQ.0))THEN
                  SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                  SideID = GlobSideTemp(ind2)
                  LocalSide = LocSidesTemp(ind2)
                  TriNum = TriNumTemp(ind2)
                  oldElemIsMortar = .TRUE.
                ELSE
                  IF(detM.EQ.0.0) CYCLE   ! For the extremely unlikely case that the particle landed exactly on the side
                  ratio = detPartPos/detM
                  ! Ratio is always negative since detM(=detLastPartPos) is negative or zero, i.e. maximum abs is wanted
                  ! The closer the intersected side is to the last particle position the greater the absolute ratio will be
                  IF (ratio.LT.minRatio) THEN
                    minRatio = ratio
                    SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                    SideID = GlobSideTemp(ind2)
                    LocalSide = LocSidesTemp(ind2)
                    TriNum = TriNumTemp(ind2)
                    oldElemIsMortar = .TRUE.
                  END IF
                END IF
              END IF  ! InElementCheck
            ELSE  ! Regular side
              CALL ParticleThroughSideLastPosCheck(i,LocSidesTemp(ind2),ElemID,InElementCheck,TriNumTemp(ind2),detM)
              IF (InElementCheck) THEN
                IF((detM.EQ.0).AND.(det(LocSidesTemp(ind2),TriNumTemp(ind2)).EQ.0)) CYCLE ! particle moves within side
                ! Determining which side was crossed first by comparing the ratio of the spatial product of PartPos->Tri-Nodes
                ! and LastPartPos->Tri-Nodes
                IF((detM.EQ.0).AND.(minRatio.EQ.0))THEN
                  SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                  SideID = GlobSideTemp(ind2)
                  LocalSide = LocSidesTemp(ind2)
                  TriNum = TriNumTemp(ind2)
                  oldElemIsMortar = .FALSE.
                ELSE
                  IF(detM.EQ.0) CYCLE   ! For the extremely unlikely case that the particle landed exactly on the side
                  ratio = det(LocSidesTemp(ind2),TriNumTemp(ind2))/detM
                  ! Ratio is always negative since detM(=detLastPartPos) is negative or zero, i.e. maximum abs is wanted
                  ! The closer the intersected side is to the last particle position the greater the absolute ratio will be
                  IF (ratio.LT.minRatio) THEN
                    minRatio = ratio
                    SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                    SideID = GlobSideTemp(ind2)
                    LocalSide = LocSidesTemp(ind2)
                    TriNum = TriNumTemp(ind2)
                    oldElemIsMortar = .FALSE.
                  END IF
                END IF
              END IF  ! InElementCheck
            END IF  ! isMortarSideTemp = T/F
          END IF  ! doCheckSide
        END DO  ! ind2 = 1, NrOfThroughSides
        ! Particle that went through multiple sides first, but did not cross any sides during the second check -> Deleted!
        IF (SecondNrOfThroughSides.EQ.0) THEN
          IF(DisplayLostParticles)THEN
            IPWRITE(*,*) 'Error in Particle TriaTracking! Particle Number',i,'lost. Element:', ElemID,'(species:',PartSpecies(i),')'
            IPWRITE(*,*) 'LastPos: ', LastPartPos(1:3,i)
            IPWRITE(*,*) 'Pos:     ', PartState(1:3,i)
            IPWRITE(*,*) 'Velo:    ', PartState(4:6,i)
            IPWRITE(*,*) 'Particle deleted!'
          END IF ! DisplayLostParticles
          IF(CountNbrOfLostParts) THEN
            CALL StoreLostParticleProperties(i, ElemID)
            NbrOfLostParticles=NbrOfLostParticles+1
          END IF
          CALL RemoveParticle(i)
          PartisDone = .TRUE.
          EXIT
        END IF
      END IF  ! NrOfThroughSides.EQ.0/.GT.1
    END IF  ! NrOfThroughSides.NE.1
    ! ----------------------------------------------------------------------------
    ! 3) In case of a boundary, perform the appropriate boundary interaction
    crossedBC=.FALSE.
    flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)
    IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
      OldElemID=ElemID
      BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
      ! Calculate the intersection with the wall and determine alpha (= fraction of trajectory to the intersection)
      IF(BCType.NE.1) CALL IntersectionWithWall(i,LocalSide,ElemID,TriNum)
      IF(PRESENT(IsInterPlanePart)) THEN
        CALL GetBoundaryInteraction(i,SideID,flip,ElemID,crossedBC,TriNum=TriNum,IsInterPlanePart=IsInterPlanePart)
      ELSE
        CALL GetBoundaryInteraction(i,SideID,flip,ElemID,crossedBC,TriNum=TriNum)
      END IF
      IF(.NOT.PDM%ParticleInside(i)) PartisDone = .TRUE.
#if USE_LOADBALANCE
      IF (OldElemID.GE.offsetElem+1.AND.OldElemID.LE.offsetElem+PP_nElems) CALL LBElemSplitTime(OldElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
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
    ELSE  ! SideInfo_Shared(SIDE_BCID,SideID).LE.0
      DO ind2= 5, 1, -1
        DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
      END DO
      DoneLastElem(1,1) = ElemID
      DoneLastElem(2,1) = LocalSide
      DoneLastElem(3,1) = TriNum
      DoneLastElem(4,1) = SideID
      IF (oldElemIsMortar) THEN
        ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      ELSE
        ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
      END IF
    END IF  ! SideInfo_Shared(SIDE_BCID,SideID).GT./.LE. 0
    CNElemID = GetCNElemID(ElemID)
    IF (CNElemID.LT.1) THEN
      IPWRITE(UNIT_StdOut,*) "VECNORM(PartState(1:3,i)-LastPartPos(1:3,i)): ", VECNORM(PartState(1:3,i)-LastPartPos(1:3,i))
      IPWRITE(UNIT_StdOut,*) " PartState(1:3,i)  : ", PartState(1:3,i)
      IPWRITE(UNIT_StdOut,*) " LastPartPos(1:3,i): ", LastPartPos(1:3,i)
      IPWRITE(UNIT_StdOut,*) " PartState(4:6,i)  : ", PartState(4:6,i)
      IPWRITE(UNIT_StdOut,*) "           ElemID: ", ElemID
      IPWRITE(UNIT_StdOut,*) "         CNElemID: ", CNElemID
      IPWRITE(UNIT_stdout,*) 'Particle Velocity: ',SQRT(DOTPRODUCT(PartState(4:6,i)))
      CALL abort(&
       __STAMP__ &
       ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
    END IF
    IF(UseRotRefFrame) THEN
      IF(PDM%ParticleInside(i)) THEN
        IF(InRotRefFrameCheck(i)) THEN
          ! Particle moved into the rotational frame of reference, initialize velocity
          IF(.NOT.InRotRefFrame(i)) THEN
            PartVeloRotRef(1:3,i) = PartState(4:6,i) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,i))
          END IF
        ELSE
          ! Particle left (or never was in) the rotational frame of reference
          PartVeloRotRef(1:3,i) = 0.
        END IF
        InRotRefFrame(i) = InRotRefFrameCheck(i)
      END IF
    END IF
    TrackInfo%CurrElem = ElemID
  END IF  ! InElementCheck = T/F
END DO  ! .NOT.PartisDone
#if USE_LOADBALANCE
IF(ParticleOnProc(i)) CALL LBElemPauseTime(PEM%LocalElemID(i),tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE SingleParticleTriaTracking3D


SUBROUTINE SingleParticleTriaTracking1D2D(i,IsInterPlanePart)
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
!    1) Perform tracking until the particle is considered "done" (either localized or deleted)
!       1a) The side through which the particle went is determined by checking each side of the element (ParticleThroughSideCheck1D2D)
!    2) In case of a boundary, determine the intersection and perform the appropriate boundary interaction (GetBoundaryInteraction)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Tracking_vars      ,ONLY: ntracks,MeasureTrackTime, TrackInfo
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound
USE MOD_Particle_Intersection       ,ONLY: ParticleThroughSideCheck1D2D
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
#if USE_LOADBALANCE
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime, LBElemSplitTime, LBElemPauseTime
USE MOD_part_tools                  ,ONLY: ParticleOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: i
LOGICAL,INTENT(IN),OPTIONAL       :: IsInterPlanePart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NblocSideID, NbElemID, CNElemID, ind, nbSideID, nMortarElems,BCType
INTEGER                          :: ElemID,flip,OldElemID,nlocSides
INTEGER                          :: LocalSide, NrOfThroughSides
INTEGER                          :: SideID,TempSideID,iLocSide, localSideID
LOGICAL                          :: ThroughSide, PartisDone
LOGICAL                          :: crossedBC, oldElemIsMortar
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
IF (MeasureTrackTime) nTracks=nTracks+1
PartisDone = .FALSE.
IF (PEM%LastGlobalElemID(i).LT.0) THEN
  TrackInfo%LastSide = -PEM%LastGlobalElemID(i)
  ElemID = PEM%GlobalElemID(i)
ELSE
  ElemID = PEM%LastGlobalElemID(i)
  TrackInfo%LastSide = 0
END IF
TrackInfo%CurrElem=ElemID
SideID = 0

! 1) Loop tracking until particle is considered "done" (either localized or deleted)
DO WHILE (.NOT.PartisDone)
  oldElemIsMortar = .FALSE.
  ! 1a) The side through which the particle went is determined by checking each side of the element (ParticleThroughSideCheck1D2D)
  NrOfThroughSides = 0
  TrackInfo%PartTrajectory(1:3)=PartState(1:3,i) - LastPartPos(1:3,i)
  TrackInfo%lengthPartTrajectory=VECNORM(TrackInfo%PartTrajectory(1:3))
  IF(ABS(TrackInfo%lengthPartTrajectory).GT.0.) TrackInfo%PartTrajectory = TrackInfo%PartTrajectory &
                                                                            / TrackInfo%lengthPartTrajectory
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  SideLoop: DO iLocSide=1,nlocSides
    TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    ! Skip symmetry side
    IF(SideIsSymSide(TempSideID)) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
    ! Side is not one of the 6 local sides
    IF (localSideID.LE.0) CYCLE
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)
    IF (NbElemID.LT.0) THEN ! Mortar side
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.-1)
      DO ind = 1, nMortarElems
        nbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
        IF (nbSideID.EQ.TrackInfo%LastSide) CYCLE
        NbElemID = SideInfo_Shared(SIDE_NBELEMID,nbSideID)
        ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
        ! no way to recover it during runtime
        IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)
        ! For small mortar sides, SIDE_NBSIDEID contains the SideID of the corresponding big mortar side
        nbSideID = SideInfo_Shared(SIDE_NBSIDEID,nbSideID)
        NblocSideID =  SideInfo_Shared(SIDE_LOCALID,nbSideID)
        ThroughSide = .FALSE.
        CALL ParticleThroughSideCheck1D2D(i,NblocSideID,NbElemID,ThroughSide)
        IF (ThroughSide) THEN
          ! Store the information for this side for future checks, if this side was already treated
          oldElemIsMortar = .TRUE.
          NrOfThroughSides = NrOfThroughSides + 1
          SideID = nbSideID
          TrackInfo%LastSide = nbSideID
          LocalSide = NblocSideID
          EXIT SideLoop
        END IF
      END DO
    ELSE  ! Regular side
      IF (TempSideID.EQ.TrackInfo%LastSide) CYCLE
      ThroughSide = .FALSE.
      CALL ParticleThroughSideCheck1D2D(i,localSideID,ElemID,ThroughSide)
      IF (ThroughSide) THEN
        NrOfThroughSides = NrOfThroughSides + 1
        SideID = TempSideID
        LocalSide = localSideID
        TrackInfo%LastSide = SideInfo_Shared(SIDE_NBSIDEID,TempSideID)
        EXIT SideLoop
      END IF
    END IF  ! Mortar or regular side
  END DO  SideLoop

  ! ----------------------------------------------------------------------------
  ! Addition treatment if particle did not cross any sides or it crossed multiple sides
  IF (NrOfThroughSides.EQ.0) THEN
    PEM%GlobalElemID(i) = ElemID
    PartisDone = .TRUE.
  ELSE
    ! ----------------------------------------------------------------------------
    ! 2) In case of a boundary, perform the appropriate boundary interaction
    crossedBC=.FALSE.
    flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)
    IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
      TrackInfo%LastSide=SideID
      OldElemID=ElemID
      BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
      ! Calculate the intersection with the wall and determine alpha (= fraction of trajectory to the intersection)
      CALL GetBoundaryInteraction(i,SideID,flip,ElemID,crossedBC,TriNum=1)
      IF(.NOT.PDM%ParticleInside(i)) PartisDone = .TRUE.
#if USE_LOADBALANCE
      IF (OldElemID.GE.offsetElem+1.AND.OldElemID.LE.offsetElem+PP_nElems) CALL LBElemSplitTime(OldElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
    ELSE  ! SideInfo_Shared(SIDE_BCID,SideID).LE.0
      IF (oldElemIsMortar) THEN
        ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      ELSE
        ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
      END IF
    END IF  ! SideInfo_Shared(SIDE_BCID,SideID).GT./.LE. 0
    CNElemID = GetCNElemID(ElemID)

    IF (CNElemID.LT.1) THEN
      IPWRITE(UNIT_StdOut,*) "VECNORM(PartState(1:3,i)-LastPartPos(1:3,i)): ", VECNORM(PartState(1:3,i)-LastPartPos(1:3,i))
      IPWRITE(UNIT_StdOut,*) " PartState(1:3,i)  : ", PartState(1:3,i)
      IPWRITE(UNIT_StdOut,*) " LastPartPos(1:3,i): ", LastPartPos(1:3,i)
      IPWRITE(UNIT_StdOut,*) " PartState(4:6,i)  : ", PartState(4:6,i)
      IPWRITE(UNIT_StdOut,*) "           ElemID: ", ElemID
      IPWRITE(UNIT_StdOut,*) "         CNElemID: ", CNElemID
      IPWRITE(UNIT_stdout,*) 'Particle Velocity: ',SQRT(DOTPRODUCT(PartState(4:6,i)))
      CALL abort(__STAMP__ ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
    END IF
    TrackInfo%CurrElem = ElemID
  END IF
END DO  ! .NOT.PartisDone
#if USE_LOADBALANCE
IF(ParticleOnProc(i)) CALL LBElemPauseTime(PEM%LocalElemID(i),tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE SingleParticleTriaTracking1D2D


SUBROUTINE SingleParticleTriaTrackingSubcycling(iPart,IsInterPlanePart)
!===================================================================================================================================
!> Store Particle information before tracking for sub-cycling.
!> it must be stored before first call of "SingleParticleTriaTracking" because particle information like LastPartPos & PartState
!> can be changed within "SingleParticleTriaTracking" e.g. in RotPeriodicInterPlaneBoundary
!> Loop over nSubCyclingSteps and perform SingleParticleTriaTracking
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars       ,ONLY: PEM,PDM,nSubCyclingSteps
USE MOD_Particle_Vars       ,ONLY: NewPosSubCycling, GlobalElemIDSubCycling, LastPartPosSubCycling
USE MOD_Particle_Vars       ,ONLY: InRotRefFrame,InRotRefFrameSubCycling, PartVeloRotRefSubCycling, LastVeloRotRefSubCycling
USE MOD_part_RHS            ,ONLY: CalcPartPosInRotRef
USE MOD_TimeDisc_Vars       ,ONLY: dt
USE MOD_Particle_Vars       ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep, Species, PartState, LastPartPos, PartSpecies
USE MOD_Particle_Vars       ,ONLY: PartVeloRotRef, LastPartVeloRotRef
!-----
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: iPart
LOGICAL,INTENT(IN),OPTIONAL :: IsInterPlanePart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iStep
REAL                        :: dtVar
LOGICAL                     :: IsInterPlanePart_loc
!===================================================================================================================================

IF(PRESENT(IsInterPlanePart)) THEN
  IsInterPlanePart_loc = IsInterPlanePart
ELSE
  IsInterPlanePart_loc = .FALSE.
END IF

LastPartPosSubCycling(1:3)    = LastPartPos(1:3,iPart)
NewPosSubCycling(1:3)         = PartState(1:3,iPart)
PartVeloRotRefSubCycling(1:3) = PartVeloRotRef(1:3,iPart)
LastVeloRotRefSubCycling(1:3) = LastPartVeloRotRef(1:3,iPart)
GlobalElemIDSubCycling        = PEM%LastGlobalElemID(iPart)
InRotRefFrameSubCycling       = InRotRefFrame(iPart)
!--- split time step in 10 sub-steps
IF (UseVarTimeStep) THEN
  dtVar = dt * PartTimeStep(iPart)
ELSE
  dtVar = dt
END IF
IF(VarTimeStep%UseSpeciesSpecific) dtVar = dtVar * Species(PartSpecies(iPart))%TimeStepFactor
dtVar = dtVar / REAL(nSubCyclingSteps)
!--- Reset Particle Push from timedisc
PartState(1:3,iPart) = LastPartPosSubCycling(1:3)
LastPartPos(1:3,iPart) = LastPartPosSubCycling(1:3)
PartVeloRotRef(1:3,iPart)=LastVeloRotRefSubCycling(1:3)
PEM%GlobalElemID(iPart) = GlobalElemIDSubCycling
PEM%LastGlobalElemID(iPart)= GlobalElemIDSubCycling
LastPartVeloRotRef(1:3,iPart) = LastVeloRotRefSubCycling(1:3)
InRotRefFrame(iPart) = InRotRefFrameSubCycling
!--- Loop over sub time steps
DO iStep = 1, nSubCyclingSteps
  IF (PDM%ParticleInside(iPart)) THEN
    IF(iStep.GT.1) THEN
      LastPartPos(1:3,iPart)=PartState(1:3,iPart)
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
      LastPartVeloRotRef(1:3,iPart)=PartVeloRotRef(1:3,iPart)
    END IF
    CALL CalcPartPosInRotRef(iPart, dtVar)
    CALL SingleParticleTriaTracking(i=iPart,IsInterPlanePart=IsInterPlanePart_loc)
  END IF
END DO
!--- Reset stored particle information
LastPartPosSubCycling    = 0.0
NewPosSubCycling         = 0.0
PartVeloRotRefSubCycling = 0.0
LastVeloRotRefSubCycling = 0.0
GlobalElemIDSubCycling   = 0
InRotRefFrameSubCycling  = .FALSE.

END SUBROUTINE SingleParticleTriaTrackingSubcycling

END MODULE MOD_Particle_TriaTracking
