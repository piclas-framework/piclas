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

MODULE MOD_Particle_Tracking
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE ParticleTriaTracking
  MODULE PROCEDURE ParticleTriaTracking
END INTERFACE

INTERFACE ParticleTracing
  MODULE PROCEDURE ParticleTracing
END INTERFACE

INTERFACE ParticleRefTracking
  MODULE PROCEDURE ParticleRefTracking
END INTERFACE

INTERFACE ParticleCollectCharges
  MODULE PROCEDURE ParticleCollectCharges
END INTERFACE

INTERFACE ParticleSanityCheck
  MODULE PROCEDURE ParticleSanityCheck
END INTERFACE

PUBLIC::ParticleTriaTracking
PUBLIC::ParticleTracing
PUBLIC::ParticleRefTracking
PUBLIC::ParticleCollectCharges
PUBLIC::ParticleSanityCheck
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


#ifdef IMPA
SUBROUTINE ParticleTriaTracking(doParticle_In)
#else
SUBROUTINE ParticleTriaTracking()
#endif /*IMPA*/
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
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Particle_Mesh               ,ONLY: SingleParticleToExactElement
USE MOD_Particle_Mesh_Tools         ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars          ,ONLY: PartElemToSide, PartSideToElem, PartSideToElem, PartElemToElemAndSide
USE MOD_Particle_Tracking_vars      ,ONLY: ntracks,MeasureTrackTime,CountNbOfLostParts,nLostParts, TrackInfo
USE MOD_Mesh_Vars                   ,ONLY: MortarType
USE MOD_Mesh_Vars                   ,ONLY: BC
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound
USE MOD_Particle_Intersection       ,ONLY: IntersectionWithWall
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_DSMC_Vars                   ,ONLY: RadialWeighting
USE MOD_DSMC_Symmetry2D             ,ONLY: DSMC_2D_RadialWeighting, DSMC_2D_SetInClones
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime, LBElemSplitTime, LBElemPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef IMPA
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
#endif /*IMPA*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef IMPA
LOGICAL                          :: doParticle
LOGICAL                          :: doPartInExists
#endif
INTEGER                          :: i, NblocSideID, NbElemID, ind, nbSideID, nMortarElems, SideIDMortar, BCType
INTEGER                          :: ElemID,flip,OldElemID
INTEGER                          :: LocalSide
INTEGER                          :: NrOfThroughSides, ind2
INTEGER                          :: SideID,TempSideID,iLocSide
INTEGER                          :: TriNum, LocSidesTemp(1:6),TriNumTemp(1:6), GlobSideTemp(1:6)
INTEGER                          :: SecondNrOfThroughSides, indSide
INTEGER                          :: DoneLastElem(1:4,1:6) ! 1:3: 1=Element,2=LocalSide,3=TriNum 1:2: 1=last 2=beforelast
LOGICAL                          :: ThroughSide, InElementCheck,PartisDone
LOGICAL                          :: crossedBC, oldElemIsMortar, isMortarSideTemp(1:6), doCheckSide
REAL                             :: det(6,2),detM,ratio,minRatio, detPartPos
REAL                             :: PartTrajectory(1:3),lengthPartTrajectory
REAL                             :: xi = -1. , eta = -1. , alpha = -1.
REAL, PARAMETER                  :: eps = 0
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#ifdef IMPA
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.
#endif /*IMPA*/

IF(RadialWeighting%DoRadialWeighting) CALL DSMC_2D_SetInClones()

! 1) Loop over all particles that are still inside
DO i = 1,PDM%ParticleVecLength
#ifdef IMPA
  IF(doPartInExists)THEN
    DoParticle=PDM%ParticleInside(i).AND.DoParticle_In(i)
  ELSE
    DoParticle=PDM%ParticleInside(i)
  END IF
  IF(DoParticle)THEN
#else
  IF (PDM%ParticleInside(i)) THEN
#endif /*IMPA*/
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    IF (MeasureTrackTime) nTracks=nTracks+1
    PartisDone = .FALSE.
    ElemID = PEM%lastElement(i)
    TrackInfo%CurrElem=ElemID
    SideID = 0
    DoneLastElem(:,:) = 0
    ! 2) Loop tracking until particle is considered "done" (either localized or deleted)
    DO WHILE (.NOT.PartisDone)
      ! 2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
      !     triangle (ParticleInsideQuad3D)
      oldElemIsMortar = .FALSE.
      CALL ParticleInsideQuad3D(PartState(i,1:3),ElemID,InElementCheck,det)
      IF (InElementCheck) THEN
        ! If particle is inside the given ElemID, set new PEM%Element and stop tracking this particle ->PartisDone
        PEM%Element(i) = ElemID
        PartisDone = .TRUE.
      ELSE
        ! 2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
        !     each side of the element (ParticleThroughSideCheck3DFast)
        NrOfThroughSides = 0
        LocSidesTemp(:) = 0
        TriNumTemp(:) = 0
        GlobSideTemp = 0
        isMortarSideTemp = .FALSE.
        PartTrajectory=PartState(i,1:3) - LastPartPos(i,1:3)
        lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                                 +PartTrajectory(2)*PartTrajectory(2) &
                                 +PartTrajectory(3)*PartTrajectory(3) )
        PartTrajectory=PartTrajectory/lengthPartTrajectory
        DO iLocSide=1,6
          TempSideID=PartElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
          SideIDMortar=MortarType(2,TempSideID)
          IF (SideIDMortar.GT.0) THEN ! Mortar side
            IF (MortarType(1,TempSideID).EQ.1) THEN
              nMortarElems = 4
            ELSE
              nMortarElems = 2
            END IF
            DO ind = 1, nMortarElems
              NbElemID = PartElemToElemAndSide(ind,iLocSide,ElemID)
              ! If small mortar element not defined, skip it for now, likely not inside the halo region (additional check is
              ! performed after the MPI communication: ParticleInsideQuad3D_MortarMPI)
              IF (NbElemID.LT.1) CYCLE
              NblocSideID = PartElemToElemAndSide(ind+4,iLocSide,ElemID)
              nbSideID = PartElemToSide(E2S_SIDE_ID,NblocSideID,NbElemID)
              DO TriNum = 1,2
                ThroughSide = .FALSE.
                CALL ParticleThroughSideCheck3DFast(i,PartTrajectory,NblocSideID,NbElemID,ThroughSide,TriNum, .TRUE.)
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
              IF (det(iLocSide,TriNum).le.-eps) THEN
                ThroughSide = .FALSE.
                CALL ParticleThroughSideCheck3DFast(i,PartTrajectory,iLocSide,ElemID,ThroughSide,TriNum)
                IF (ThroughSide) THEN
                  NrOfThroughSides = NrOfThroughSides + 1
                  LocSidesTemp(NrOfThroughSides) = iLocSide
                  TriNumTemp(NrOfThroughSides) = TriNum
                  GlobSideTemp(NrOfThroughSides) = TempSideID
                  SideID = TempSideID
                  LocalSide = iLocSide
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
            IPWRITE(*,*) 'Error in Particle TriaTracking! Particle Number',i,'lost. Element:', ElemID
            IPWRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
            IPWRITE(*,*) 'Pos:     ', PartState(i,1:3)
            IPWRITE(*,*) 'Velo:    ', PartState(i,4:6)
            IPWRITE(*,*) 'Particle deleted!'
            PDM%ParticleInside(i) = .FALSE.
            IF(CountNbOfLostParts) nLostParts=nLostParts+1
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
                  IF (PartSideToElem(S2E_ELEM_ID,GlobSideTemp(ind2)).GT.0) THEN
                    NbElemID = PartSideToElem(S2E_ELEM_ID,GlobSideTemp(ind2))
                  ELSE
                    NbElemID = PartSideToElem(S2E_NB_ELEM_ID,GlobSideTemp(ind2))
                  END IF
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
                      SideID = PartElemToSide(E2S_SIDE_ID,LocSidesTemp(ind2),ElemID)
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
                        SideID = PartElemToSide(E2S_SIDE_ID,LocSidesTemp(ind2),ElemID)
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
              IPWRITE(*,*) 'Error in Particle TriaTracking! Particle Number',i,'lost on second check. Element:', ElemID
              IPWRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
              IPWRITE(*,*) 'Pos:     ', PartState(i,1:3)
              IPWRITE(*,*) 'Velo:    ', PartState(i,4:6)
              IPWRITE(*,*) 'Particle deleted!'
              PDM%ParticleInside(i) = .FALSE.
              IF(CountNbOfLostParts) nLostParts=nLostParts+1
              PartisDone = .TRUE.
              EXIT
            END IF
          END IF  ! NrOfThroughSides.EQ.0/.GT.1
        END IF  ! NrOfThroughSides.NE.1
        ! ----------------------------------------------------------------------------
        ! 3) In case of a boundary, perform the appropriate boundary interaction
        crossedBC=.FALSE.
        flip =PartElemToSide(E2S_FLIP,LocalSide,ElemID)
        IF(BC(SideID).GT.0) THEN
          OldElemID=ElemID
          BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID)))
          IF(BCType.NE.1) CALL IntersectionWithWall(PartTrajectory,alpha,i,LocalSide,ElemID,TriNum)
          CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha &
                                                                       ,xi    &
                                                                       ,eta   ,i,SideID,flip,LocalSide,ElemID,crossedBC&
                                                                       ,TriNum)
          IF(.NOT.PDM%ParticleInside(i)) PartisDone = .TRUE.
#if USE_LOADBALANCE
          IF (OldElemID.LE.PP_nElems) CALL LBElemSplitTime(OldElemID,tLBStart)
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
        ELSE  ! BC(SideID).LE.0
          DO ind2= 5, 1, -1
            DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
          END DO
          DoneLastElem(1,1) = ElemID
          DoneLastElem(2,1) = LocalSide
          DoneLastElem(3,1) = TriNum
          DoneLastElem(4,1) = SideID
          IF (oldElemIsMortar) THEN
            IF (PartSideToElem(S2E_NB_ELEM_ID,SideID).EQ.-1) THEN
             ElemID = PartSideToElem(S2E_ELEM_ID,SideID)
            ELSE
             ElemID = PartSideToElem(S2E_NB_ELEM_ID,SideID)
            END IF
          ELSE
            ElemID = PartElemToElemAndSide(1  ,LocalSide,ElemID)
          END IF
        END IF  ! BC(SideID).GT./.LE. 0
        IF (ElemID.LT.1) THEN
          IPWRITE(UNIT_stdout,*) 'Particle Velocity: ',SQRT(DOT_PRODUCT(PartState(i,4:6),PartState(i,4:6)))
          CALL abort(&
           __STAMP__ &
           ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
        END IF
        TrackInfo%CurrElem = ElemID
      END IF  ! InElementCheck = T/F
    END DO  ! .NOT.PartisDone
#if USE_LOADBALANCE
    IF (PEM%Element(i).LE.PP_nElems) CALL LBElemPauseTime(PEM%Element(i),tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
  ! Particle treatment for an axisymmetric simulation (cloning/deleting particles)
  IF(RadialWeighting%DoRadialWeighting) THEN
    IF(PDM%ParticleInside(i)) THEN
      IF (PEM%Element(i).LE.PP_nElems) CALL DSMC_2D_RadialWeighting(i,PEM%Element(i))
    END IF
  END IF
END DO ! i = 1,PDM%ParticleVecLength

END SUBROUTINE ParticleTriaTracking


#ifdef IMPA
SUBROUTINE ParticleTracing(doParticle_In,nInnerNewton_In)
#else
SUBROUTINE ParticleTracing()
#endif /*NOT IMPA*/
!===================================================================================================================================
! Routine for tracing moving particles, calculate intersection and boundary interaction
! in case of no reference tracking (dorefmapping = false)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Mesh_Vars          ,ONLY: PartElemToSide,ElemType,ElemRadiusNGeo,ElemHasAuxBCs
USE MOD_Particle_Boundary_Vars      ,ONLY: nAuxBCs,UseAuxBCs
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteractionAuxBC
USE MOD_Utils                       ,ONLY: InsertionSort
USE MOD_Particle_Tracking_vars      ,ONLY: ntracks, MeasureTrackTime, CountNbOfLostParts , nLostParts
USE MOD_Particle_Mesh               ,ONLY: SingleParticleToExactElementNoMap,PartInElemCheck
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectInterSection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeAuxBCIntersection
USE MOD_Mesh_Vars                   ,ONLY: OffSetElem
USE MOD_Eval_xyz                    ,ONLY: GetPositionInRefElem
#if USE_MPI
USE MOD_Particle_MPI_Vars           ,ONLY: PartHaloElemToProc
USE MOD_MPI_Vars                    ,ONLY: offsetElemMPI
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
#ifdef IMPA
USE MOD_Particle_Vars               ,ONLY: PartIsImplicit,PartDtFrac
USE MOD_Particle_Vars               ,ONLY: PartStateN
#endif /*IMPA*/
USE MOD_Particle_Intersection       ,ONLY: OutputTrajectory
USE MOD_Particle_Tracking_Vars      ,ONLY: PartOut,MPIRankOut
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO
USE MOD_TimeDisc_Vars               ,ONLY: iStage
USE MOD_Globals_Vars                ,ONLY: epsMach
USE MOD_Mesh_Vars                   ,ONLY: ElemBaryNGeo
#endif /*CODE_ANALYZE*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime,LBElemPauseTime,LBElemSplitTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef IMPA
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
INTEGER,INTENT(IN),OPTIONAL   :: nInnerNewton_In
#endif /*IMPA*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef IMPA
LOGICAL                       :: doParticle
LOGICAL                       :: doPartInExists
#endif /*IMPA*/
INTEGER                       :: iPart,ElemID,flip,OldElemID,firstElem,iAuxBC,AuxBCsToCheck
INTEGER                       :: ilocSide,SideID, locSideList(1:6), hitlocSide,nInterSections
LOGICAL                       :: PartisDone,dolocSide(1:6),isHit,markTol,crossedBC,SwitchedElement,isCriticalParallelInFace
LOGICAL                       :: HasAuxBC,OnlyAuxBC,IsIntersec,IsAuxBC
REAL                          :: localpha(1:6),xi(1:6),eta(1:6),refpos(1:3)
REAL,ALLOCATABLE              :: locAlphaAll(:)
INTEGER,ALLOCATABLE           :: locListAll(:)
!INTEGER                       :: lastlocSide
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
#if USE_LOADBALANCE
REAL                          :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
INTEGER                       :: PartDoubleCheck
REAL                          :: alphaOld
#if USE_MPI
INTEGER                       :: inElem
#endif
#ifdef CODE_ANALYZE
REAL                          :: IntersectionPoint(1:3)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

IF (UseAuxBCs) THEN
  ALLOCATE(locAlphaAll(1:6+nAuxBCs) &
    ,locListAll(1:6+nAuxBCs))
END IF

#ifdef IMPA
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.
#endif /*IMPA*/

! IF(RadialWeighting%DoRadialWeighting) CALL DSMC_2D_SetInClones()

DO iPart=1,PDM%ParticleVecLength
  PartDoubleCheck=0
  alphaOld = -1.0
#ifdef IMPA
  IF(doPartInExists)THEN
    DoParticle=PDM%ParticleInside(iPart).AND.DoParticle_In(iPart)
  ELSE
    DoParticle=PDM%ParticleInside(iPart)
  END IF
  IF(DoParticle)THEN
#else
  IF (PDM%ParticleInside(iPart)) THEN
#endif /*IMPA*/
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    IF (MeasureTrackTime) nTracks=nTracks+1
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(iPart)
#ifdef CODE_ANALYZE
    IF(GEO%nPeriodicVectors.EQ.0)THEN
      IF(   (LastPartPos(iPart,1).GT.GEO%xmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(iPart,1),GEO%xmaxglob) &
        .OR.(LastPartPos(iPart,1).LT.GEO%xminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(iPart,1),GEO%xminglob) &
        .OR.(LastPartPos(iPart,2).GT.GEO%ymaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(iPart,2),GEO%ymaxglob) &
        .OR.(LastPartPos(iPart,2).LT.GEO%yminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(iPart,2),GEO%yminglob) &
        .OR.(LastPartPos(iPart,3).GT.GEO%zmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(iPart,3),GEO%zmaxglob) &
        .OR.(LastPartPos(iPart,3).LT.GEO%zminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(iPart,3),GEO%zminglob) ) THEN
        IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(iPart)
#ifdef IMPA
        IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(iPart)
        IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(iPart)
        IF(PRESENT(nInnerNewton_In))THEN
          IPWRITE(UNIt_stdOut,'(I0,A18,I0)')                           ' nInnerNewton_In ', nInnerNewton_In
        END IF
#endif /*IMPA*/
        IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
        IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
        IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, LastPartPos(iPart,1), GEO%xmaxglob
        IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, LastPartPos(iPart,2), GEO%ymaxglob
        IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, LastPartPos(iPart,3), GEO%zmaxglob
        CALL abort(&
           __STAMP__ &
           ,' LastPartPos outside of mesh. iPart=, iStage',iPart,REAL(iStage))
      END IF
    END IF
#endif /*CODE_ANALYZE*/
    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
    lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             +PartTrajectory(2)*PartTrajectory(2) &
                             +PartTrajectory(3)*PartTrajectory(3) )

    IF(.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(ElemID)))THEN
      PEM%Element(iPart)=ElemID
      PartisDone=.TRUE.
      CYCLE
    END IF
    PartTrajectory=PartTrajectory/lengthPartTrajectory
#ifdef CODE_ANALYZE
    IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
      IF(iPart.EQ.PARTOUT)THEN
        WRITE(UNIT_stdout,'(A32)')         ' ---------------------------------------------------------------'
        WRITE(UNIT_stdout,'(A)')         '     | Output of Particle information '
        CALL OutputTrajectory(iPart,PartState(iPart,1:3),PartTrajectory,lengthPartTrajectory)
#if USE_MPI
        InElem=PEM%LastElement(iPart)
        IF(InElem.LE.PP_nElems)THEN
          WRITE(UNIT_stdOut,'(A,I0)') '     | global ElemID       ', InElem+offSetElem
        ELSE
          WRITE(UNIT_stdOut,'(A,I0)') '     | global ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
        END IF
#else
        WRITE(UNIT_stdOut,'(A,I0)') '     | global ElemID         ', PEM%LastElement(iPart)+offSetElem
#endif /*USE_MPI*/
      END IF
    END IF
    ! caution: reuse of variable, isHit=TRUE == inside
    CALL PartInElemCheck(LastPartPos(iPart,1:3),iPart,ElemID,isHit,IntersectionPoint,CodeAnalyze_Opt=.TRUE.)
    IF(.NOT.isHit)THEN  ! particle not inside
     IPWRITE(UNIT_stdOut,'(I0,A)') ' LastPartPos not inside of element! '
#ifdef IMPA
     IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
     IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%ParticleInside ', PDM%ParticleInside(iPart)
#endif /*IMPA*/
     IF(ElemID.LE.PP_nElems)THEN
       IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', ElemID+offSetElem
     ELSE
#if USE_MPI
       IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,ElemID)) &
                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,ElemID)
#endif /*USE_MPI*/
     END IF
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,ElemID)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' IntersectionPoint: ', IntersectionPoint
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos:       ', LastPartPos(iPart,1:3)
#ifdef IMPA
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartStateN:       ', PartStateN(iPart,1:3)
#endif /*IMPA*/
     IPWRITE(UNIT_stdOut,'(I0,A)')            ' PartPos:           '
     IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
     IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(iPart,1), GEO%xmaxglob
     IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(iPart,2), GEO%ymaxglob
     IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(iPart,3), GEO%zmaxglob
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartTrajectory:    ', PartTrajectory
     IPWRITE(UNIT_stdOut,'(I0,A,E15.8)')      ' lengthPT:          ', lengthPartTrajectory
     CALL abort(&
     __STAMP__ &
     ,'iPart=. ',iPart)
    END IF
#endif /*CODE_ANALYZE*/
    ! track particle vector until the final particle position is achieved
    dolocSide=.TRUE.
    firstElem=ElemID
    OnlyAuxBC=.FALSE.
    HasAuxBC=.FALSE.
    IF (UseAuxBCs) THEN
      IF (ANY(ElemHasAuxBCs(ElemID,:))) THEN
        HasAuxBC=.TRUE.
      END IF
    END IF
    IF (ElemType(ElemID).EQ.1) THEN
      !removed CheckPlanarInside since it can be inconsistent for planar-assumed sides:
      !they can still be planar-nonrect for which the bilin-algorithm will be used which might give a different result
      !(anyway, this was a speed-up for completely planar meshes only, but those should be now calculated with triatracking)
      !CALL CheckPlanarInside(iPart,ElemID,lengthPartTrajectory,PartisDone)
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(110("="))')
          WRITE(UNIT_stdout,'(A,L)')    '     | Elem has AuxBC: ',HasAuxBC
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      IF (PartisDone) THEN
        IF (HasAuxBC) THEN
          OnlyAuxBC=.TRUE.
          PartisDone=.FALSE.
        ELSE
          PEM%Element(iPart) = ElemID
          CYCLE
        END IF !HasAuxBC
      END IF !inside
    END IF !planar elem
    markTol =.FALSE.
    DO WHILE (.NOT.PartisDone)
      locAlpha=-1.
      nInterSections=0
      markTol =.FALSE.
      DO ilocSide=1,6
        IF (HasAuxBC) THEN
          locListAll(ilocSide)=ilocSide
          IF (OnlyAuxBC) CYCLE
        END IF
        locSideList(ilocSide)=ilocSide
        IF(.NOT.dolocSide(ilocSide)) CYCLE
        !SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID)
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID)
        flip  =PartElemToSide(E2S_FLIP,ilocSide,ElemID)
        isCriticalParallelInFace=.FALSE.
        IF (PartDoubleCheck.EQ.1) THEN
#ifdef CODE_ANALYZE
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
            IF(iPart.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(110("="))')
              WRITE(UNIT_stdout,'(A)')    '     | Particle is double checked: '
            END IF
          END IF
#endif /*CODE_ANALYZE*/
          SELECT CASE(SideType(SideID))
          CASE(PLANAR_RECT)
            CALL ComputePlanarRectInterSection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide)   &
                                                                                          ,xi (ilocSide)      &
                                                                                          ,eta(ilocSide)      &
                                                                                          ,iPart,flip,SideID  &
                                                                                          ,isCriticalParallelInFace)
          CASE(BILINEAR,PLANAR_NONRECT)
            CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                          ,xi (ilocSide)      &
                                                                                          ,eta(ilocSide)      &
                                                                                          ,iPart,SideID,alpha2=alphaOld)
          CASE(PLANAR_CURVED)
            CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                          ,xi (ilocSide)      &
                                                                                          ,eta(ilocSide)   ,iPart,flip,SideID &
                                                                                          ,isCriticalParallelInFace)

          CASE(CURVED)
            CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                    ,xi (ilocSide)      &
                                                                                    ,eta(ilocSide)      ,iPart,SideID &
                                                                                    ,isCriticalParallelInFace)
          CASE DEFAULT
            CALL abort(&
            __STAMP__ &
            ,' Missing required side-data. Please increase halo region. ',SideID)
          END SELECT
        ELSE
          SELECT CASE(SideType(SideID))
          CASE(PLANAR_RECT)
            CALL ComputePlanarRectInterSection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide)   &
                                                                                          ,xi (ilocSide)      &
                                                                                          ,eta(ilocSide)      &
                                                                                          ,iPart,flip,SideID  &
                                                                                          ,isCriticalParallelInFace)
          CASE(BILINEAR,PLANAR_NONRECT)
            CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                          ,xi (ilocSide)      &
                                                                                          ,eta(ilocSide)      &
                                                                                          ,iPart,SideID)
          CASE(PLANAR_CURVED)
            CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                          ,xi (ilocSide)      &
                                                                                          ,eta(ilocSide)   ,iPart,flip,SideID &
                                                                                          ,isCriticalParallelInFace)

          CASE(CURVED)
            CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                    ,xi (ilocSide)      &
                                                                                    ,eta(ilocSide)      ,iPart,SideID &
                                                                                    ,isCriticalParallelInFace)
          CASE DEFAULT
            CALL abort(&
            __STAMP__ &
            ,' Missing required side-data. Please increase halo region. ',SideID)
          END SELECT
        END IF
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(iPart.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(30("-"))')
            WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (particle tracing): '
            WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,' | Hit: ',isHit
            WRITE(UNIT_stdout,'(2(A,G0))') '     | Alpha: ',locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
            WRITE(UNIT_stdout,'(A,2(X,G0))') '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
          END IF
        END IF
#endif /*CODE_ANALYZE*/
        IF(isCriticalParallelInFace)THEN
          IPWRITE(UNIT_stdOut,'(I0,A)') ' Warning: Particle located inside of face and moves parallel to side. Undefined position. '
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Removing particle with id: ',iPart
          PartIsDone=.TRUE.
          PDM%ParticleInside(iPart)=.FALSE.
#ifdef IMPA
          DoParticle=.FALSE.
#endif /*IMPA*/
          IF(CountNbOfLostParts) nLostParts=nLostParts+1
          EXIT
        END IF
        IF(isHit) THEN
          nInterSections=nInterSections+1
          IF((ABS(xi(ilocSide)).GE.0.99).OR.(ABS(eta(ilocSide)).GE.0.99)) markTol=.TRUE.
          IF(ALMOSTZERO(locAlpha(ilocSide))) markTol=.TRUE.
          IF(locAlpha(ilocSide)/lengthPartTrajectory.GE.0.99) markTol=.TRUE.
        END IF
      END DO ! ilocSide
      IF (HasAuxBC) THEN
        locAlphaAll=-1.
        DO iAuxBC=1,nAuxBCs
          locListAll(6+iAuxBC)=6+iAuxBC
          isCriticalParallelInFace=.FALSE.
          IF (ElemHasAuxBCs(ElemID,iAuxBC)) THEN
            CALL ComputeAuxBCIntersection(isHit,PartTrajectory,lengthPartTrajectory &
              ,iAuxBC,locAlphaAll(6+iAuxBC) &
              ,iPart &
              ,isCriticalParallelInFace)
          ELSE
            isHit=.FALSE.
          END IF
#ifdef CODE_ANALYZE
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
            IF(iPart.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(30("-"))')
              WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (particle tracing): '
              WRITE(UNIT_stdout,'(A,I0,A,L)') '     | AuxBC: ',iAuxBC,' | Hit: ',isHit
              WRITE(UNIT_stdout,'(2(A,G0))') '     | Alpha: ',locAlphaAll(6+iAuxBC),' | LengthPartTrajectory: ',lengthPartTrajectory
            END IF
          END IF
#endif /*CODE_ANALYZE*/
          IF(isCriticalParallelInFace)THEN
            IPWRITE(UNIT_stdOut,'(I0,A)') ' Warning: Particle located inside of BC and moves parallel to side. Undefined position. '
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Removing particle with id: ',iPart
            PartIsDone=.TRUE.
            PDM%ParticleInside(iPart)=.FALSE.
#ifdef IMPA
            DoParticle=.FALSE.
#endif /*IMPA*/
            IF(CountNbOfLostParts) nLostParts=nLostParts+1
            EXIT
          END IF
          IF(isHit) THEN
            nInterSections=nInterSections+1
            IF(ALMOSTZERO(locAlphaAll(6+iAuxBC))) markTol=.TRUE.
          END IF
        END DO !iAuxBC
      END IF !HasAuxBC

#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A,I0)') '     > Number of found intersections: ',nIntersections
          IF(markTol)THEN
            WRITE(UNIT_stdout,'(A)') '     | Tolerance marked ... '
          END IF
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      SELECT CASE(nInterSections)
      CASE(0) ! no intersection
        PEM%Element(iPart)=ElemID
        PartisDone=.TRUE.
        SwitchedElement=.FALSE.
        crossedBC=.FALSE.
      CASE(1) ! one intersection
        ! get intersection side
        SwitchedElement=.FALSE.
        crossedBC=.FALSE.
        DO ilocSide=1,6
          IF(locAlpha(ilocSide).GT.-1.0) THEN
            IF (PartDoubleCheck.EQ.0) THEN
              alphaOld = locAlpha(ilocSide)
            END IF
            hitlocSide=ilocSide
            SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
            flip  =PartElemToSide(E2S_FLIP,hitlocSide,ElemID)
            OldElemID=ElemID
            CALL SelectInterSectionType(PartIsDone,crossedBC,doLocSide,flip,hitlocSide,ilocSide,PartTrajectory &
              ,lengthPartTrajectory,xi(hitlocSide),eta(hitlocSide),localpha(ilocSide),iPart,SideID,SideType(SideID),ElemID)
            IF(ElemID.NE.OldElemID)THEN
              ! particle moves in new element, do not check yet, because particle may encounter a boundary condition
              SwitchedElement=.TRUE.
              IF(ALMOSTZERO(lengthPartTrajectory))THEN
                PartisDone=.TRUE.
              END IF
#if USE_LOADBALANCE
              IF (OldElemID.LE.PP_nElems) CALL LBElemSplitTime(OldElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
              EXIT
            END IF
            IF(crossedBC) THEN
              firstElem=ElemID
              EXIT
            END IF
          END IF
        END DO ! ilocSide
        !-- check for AuxBC interactions (if one exists, ilocSide-loop could not have found one, since nInterSections=1)
        IF (HasAuxBC) THEN
          DO iAuxBC=1,nAuxBCs
            IF(locAlphaAll(6+iAuxBC).GT.-1.0) THEN
              CALL GetBoundaryInteractionAuxBC(PartTrajectory,lengthPartTrajectory,locAlphaAll(6+iAuxBC),iPart,iAuxBC,crossedBC)
              IF(.NOT.PDM%ParticleInside(iPart)) PartisDone = .TRUE.
              dolocSide=.TRUE. !important when before there was an elemchange !
              OnlyAuxBC=.FALSE. !important, since a new elem could have been reached now !
              IF(crossedBC) THEN
                firstElem=ElemID
                EXIT
              END IF
            END IF
          END DO !iAuxBC
        END IF !HasAuxBC

        IF((.NOT.crossedBC).AND.(.NOT.SwitchedElement)) THEN
          IF (PartDoubleCheck.EQ.0) THEN
            PartDoubleCheck = 1
            PartIsDone= .FALSE.
          ELSE
            PartIsDone= .TRUE.
            PEM%Element(iPart)=ElemID !periodic BC always exits with one hit from outside
            EXIT
          END IF
        ELSE !IF(CrossedBC.OR.SwitchedElem)
          IF (PartDoubleCheck.EQ.1) THEN
            PartDoubleCheck=0
            alphaOld = -1.0
          END IF
        END IF
#ifdef IMPA
        IF(CrossedBC)THEN
          IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
        END IF
#endif /*IMPA*/
      CASE DEFAULT ! two or more hits
        ! more careful witEh bc elems
          IF (HasAuxBC) THEN
            locAlphaAll(1:6)=locAlpha
            CALL InsertionSort(locAlphaAll,locListAll,6+nAuxBCs)
            AuxBCsToCheck=nAuxBCs
          ELSE
            CALL InsertionSort(locAlpha,locSideList,6)
            AuxBCsToCheck=0
          END IF
          IF (PartDoubleCheck.EQ.0) THEN
            DO iLocSide=1,6+nAuxBCs
              IF (HasAuxBC) THEN
                IF (locAlphaAll(ilocSide).GT.-1.0) THEN
                  alphaOld = locAlphaAll(ilocSide)
                  EXIT
                END IF
              ELSE
                IF (locAlpha(ilocSide).GT.-1.0) THEN
                  alphaOld = locAlpha(ilocSide)
                  EXIT
                END IF
              END IF
            END DO
          END IF
          SwitchedElement=.FALSE.
          crossedBC=.FALSE.
          DO ilocSide=1,6+AuxBCsToCheck
            IsIntersec=.FALSE.
            IsAuxBC=.FALSE.
            IF (HasAuxBC) THEN
              IF (locAlphaAll(ilocSide).GT.-1.0) IsIntersec=.TRUE.
              IF (locListAll(ilocSide).GT.6) IsAuxBC=.TRUE.
            ELSE
              IF (locAlpha(ilocSide).GT.-1.0) IsIntersec=.TRUE.
            END IF
            IF(IsIntersec)THEN
              IF (.NOT.IsAuxBC) THEN
                IF (HasAuxBC) THEN
                  hitlocSide=locListAll(ilocSide)
                ELSE
                  hitlocSide=locSideList(ilocSide)
                END IF
                SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
                flip  =PartElemToSide(E2S_FLIP,hitlocSide,ElemID)
                OldElemID=ElemID
                IF (HasAuxBC) THEN
                  CALL SelectInterSectionType(PartIsDone,crossedBC,doLocSide,flip,hitlocSide,ilocSide,PartTrajectory &
                    ,lengthPartTrajectory,xi(hitlocSide),eta(hitlocSide),locAlphaAll(ilocSide),iPart,SideID,SideType(SideID),ElemID)
                ELSE
                  CALL SelectInterSectionType(PartIsDone,crossedBC,doLocSide,flip,hitlocSide,ilocSide,PartTrajectory &
                    ,lengthPartTrajectory,xi(hitlocSide),eta(hitlocSide),localpha(ilocSide),iPart,SideID,SideType(SideID),ElemID)
                END IF

                IF(ElemID.NE.OldElemID)THEN
                  IF(.NOT.CrossedBC) SwitchedElement=.TRUE.
                  IF(ALMOSTZERO(lengthPartTrajectory))THEN
                    PartisDone=.TRUE.
                  END IF
                  !PartTrajectory=PartTrajectory/lengthPartTrajectory
#if USE_LOADBALANCE
                  IF (OldElemID.LE.PP_nElems) CALL LBElemSplitTime(OldElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
                  !EXIT
                END IF
                IF(SwitchedElement) EXIT
                IF(crossedBC) THEN
                  firstElem=ElemID
                  EXIT
                END IF
              ELSE !IsAuxBC=.TRUE.
                CALL GetBoundaryInteractionAuxBC(&
                  PartTrajectory,lengthPartTrajectory,locAlphaAll(ilocSide),iPart,locListAll(ilocSide)-6,crossedBC)
                IF(.NOT.PDM%ParticleInside(iPart)) PartisDone = .TRUE.
                dolocSide=.TRUE. !important when before there was an elemchange !
                OnlyAuxBC=.FALSE. !important, since a new elem could have been reached now !
                IF(crossedBC) THEN
                  firstElem=ElemID
                  EXIT
                END IF
              END IF !IsAuxBC
            END IF !IsIntersec
          END DO ! ilocSide
          IF((.NOT.crossedBC).AND.(.NOT.SwitchedElement)) THEN
            IF (PartDoubleCheck.EQ.0) THEN
              PartDoubleCheck = 1
              PartIsDone= .FALSE.
            ELSE
              PartIsDone= .TRUE.
              EXIT
            END IF
          ELSE !IF(CrossedBC.OR.SwitchedElem)
            IF (PartDoubleCheck.EQ.1) THEN
              PartDoubleCheck=0
              alphaOld = -1.0
            END IF
          END IF
#ifdef IMPA
          IF(CrossedBC)THEN
            IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
          END IF
#endif /*IMPA*/
         !! particle moves close to an edge or corner. this is a critical movement because of possible tolerance issues
!        END IF
      END SELECT
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)') '     | Output of new Element after intersections number check: '
          WRITE(UNIT_stdout,'(A,L,A,L,A,L)') '     | crossed Side: ',crossedBC,' switched Element: ',SwitchedElement,&
                  ' Particle tracking done: ',PartisDone
          IF(SwitchedElement) THEN
            WRITE(UNIT_stdout,'(A,I0,A,I0)') '     | First_ElemID: ',PEM%LastElement(iPart),' | new Element: ',ElemID
#if USE_MPI
            InElem=PEM%LastElement(iPart)
            IF(InElem.LE.PP_nElems)THEN
              WRITE(UNIT_stdOut,'(A,I0)') '     | first global ElemID       ', InElem+offSetElem
            ELSE
              WRITE(UNIT_stdOut,'(A,I0)') '     | first global ElemID       ' &
                , offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
            END IF
#else
            WRITE(UNIT_stdOut,'(A,I0)') '     | first global ElemID         ', PEM%LastElement(iPart)+offSetElem
#endif
#if USE_MPI
            InElem=ElemID
            IF(InElem.LE.PP_nElems)THEN
              WRITE(UNIT_stdOut,'(A,I0)') '     | new global ElemID       ', InElem+offSetElem
            ELSE
              WRITE(UNIT_stdOut,'(A,I0)') '     | new global ElemID       ' &
                , offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
            END IF
#else
            WRITE(UNIT_stdOut,'(A,I0)') '     | new global ElemID         ', ElemID+offSetElem
#endif
          END IF
          IF( crossedBC) THEN
            WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Last    PartPos:       ',lastPartPos(iPart,1:3)
            WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Current PartPos:       ',PartState(iPart,1:3)
            WRITE(UNIT_stdout,'(A,3(X,G0))') '     | PartTrajectory:        ',PartTrajectory(1:3)
            WRITE(UNIT_stdout,'(A,(G0))')    '     | Length PartTrajectory: ',lengthPartTrajectory
          END IF
        END IF
      END IF
#endif /*CODE_ANALYZE*/
    END DO ! PartisDone=.FALSE.
    IF(markTol)THEN
      IF(.NOT.PDM%ParticleInside(iPart))THEN
#ifdef IMPA
        DoParticle=.FALSE.
#endif /*IMPA*/
        CYCLE !particle is outside cell
      END IF
      CALL PartInElemCheck(PartState(iPart,1:3),iPart,ElemID,isHit)
      PEM%Element(iPart)=ElemID
      IF(.NOT.isHit) THEN
        IPWRITE(UNIT_stdOut,'(I0,A)') '     | Relocating....'
        CALL SingleParticleToExactElementNoMap(iPart,doHALO=.TRUE.,doRelocate=.FALSE.)!debug=.TRUE.)
      END IF
      PartIsDone=.TRUE.
      IF(.NOT.PDM%ParticleInside(iPart))THEN
        !WRITE(UNIT_stdOut,'(20(=))')
        IPWRITE(UNIT_stdOut,'(I0,A)') '     | Tolerance Issue during tracing! '
        IPWRITE(UNIT_stdOut,'(I0,2(A,I0))') '     | Proc: ',MyRank,' lost particle with ID', iPart
        IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') '     | LastPartPos: ',LastPartPos(ipart,1:3)
        IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') '     |     PartPos: ',PartState(ipart,1:3)
        IPWRITE(UNIT_stdOut,'(I0,A)') '     | Computing PartRefPos ... '
        CALL GetPositionInRefElem(LastPartPos(iPart,1:3),refpos(1:3),PEM%lastElement(ipart))
        IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') '     | LastPartRefPos: ',refpos
        CALL GetPositionInRefElem(PartState(iPart,1:3),refpos(1:3),PEM%lastElement(ipart))
        IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') '     |     PartRefPos: ',refpos
        !WRITE(UNIT_stdOut,'(20(=))')
#if USE_MPI
        InElem=PEM%Element(iPart)
        IF(InElem.LE.PP_nElems)THEN
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') '     | ElemID       ', InElem+offSetElem
        ELSE
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') '     | ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                                + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
        END IF
#else
        !IPWRITE(UNIT_stdOut,*) ' ElemID         ', InElem+offSetElem  ! old
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') '     | ElemID         ', ElemID+offSetElem   ! new
#endif
#if USE_MPI
        InElem=PEM%LastElement(iPart)
        IF(InElem.LE.PP_nElems)THEN
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') '     | Last-ElemID  ', InElem+offSetElem
        ELSE
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') '     | Last-ElemID  ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                    + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
        END IF
#else
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') '     | Last-ElemID    ', ElemID+offSetElem
#endif
        IF(CountNbOfLostParts) nLostParts=nLostParts+1
      END IF
    END IF ! markTol
#if USE_LOADBALANCE
    IF (PEM%Element(iPart).LE.PP_nElems) CALL LBElemPauseTime(PEM%Element(iPart),tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF ! Part inside
  ! Particle treatment for an axisymmetric simulation (cloning/deleting particles)
  ! IF(RadialWeighting%DoRadialWeighting) THEN
  !   IF(PDM%ParticleInside(iPart)) THEN
  !     CALL DSMC_2D_RadialWeighting(iPart,PEM%Element(iPart))
  !   END IF
  ! END IF
END DO ! iPart

#ifdef CODE_ANALYZE
! check if particle is still inside of bounding box of domain and in element
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/
DO iPart=1,PDM%ParticleVecLength
#ifdef IMPA
  IF(doPartInExists)THEN
    DoParticle=PDM%ParticleInside(iPart).AND.DoParticle_In(iPart)
  ELSE
    DoParticle=PDM%ParticleInside(iPart)
  END IF
  IF(DoParticle)THEN
#else
  IF (PDM%ParticleInside(iPart)) THEN
#endif /*IMPA*/
    IF( (PartState(iPart,1).GT.GEO%xmaxglob) &
    .OR.(PartState(iPart,1).LT.GEO%xminglob) &
    .OR.(PartState(iPart,2).GT.GEO%ymaxglob) &
    .OR.(PartState(iPart,2).LT.GEO%yminglob) &
    .OR.(PartState(iPart,3).GT.GEO%zmaxglob) &
    .OR.(PartState(iPart,3).LT.GEO%zminglob) ) THEN
#ifdef IMPA
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' DoParticle ', DoParticle
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(iPart)
      IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(iPart)
      IF(PRESENT(nInnerNewton_In))THEN
        IPWRITE(UNIt_stdOut,'(I0,A18,I0)')                           ' nInnerNewton_In ', nInnerNewton_In
      END IF
#endif /*IMPA*/
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
      IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' LastPosition   ', LastPartPos(iPart,1:3)
      IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' Velocity       ', PartState(iPart,4:6)
      IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
      IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(iPart,1), GEO%xmaxglob
      IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(iPart,2), GEO%ymaxglob
      IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(iPart,3), GEO%zmaxglob
      CALL abort(&
     __STAMP__ &
     ,' PartPos outside of mesh AFTER tracking. iPart= ,iStage= ',iPart,REAL(iStage))
    END IF
    ! caution: reuse of variable, isHit=TRUE == inside
    ElemID=PEM%Element(iPart)
    CALL PartInElemCheck(PartState(iPart,1:3),iPart,ElemID,isHit,IntersectionPoint,CodeAnalyze_Opt=.TRUE.)
    IF(.NOT.isHit)THEN  ! particle not inside
     IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
     IF(ElemID.LE.PP_nElems)THEN
       IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', ElemID+offSetElem
     ELSE
#if USE_MPI
       IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,ElemID)) &
                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,ElemID)
#endif /*USE_MPI*/
     END IF
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,ElemID)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' IntersectionPoint: ', IntersectionPoint
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos:       ', LastPartPos(iPart,1:3)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartPos:           ', PartState(iPart,1:3)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartTrajectory:    ', PartTrajectory
     IPWRITE(UNIT_stdOut,'(I0,A,E15.8)')      ' lengthPT:          ', lengthPartTrajectory
     CALL abort(&
     __STAMP__ &
     ,'iPart=. ',iPart)
    END IF
  END IF ! Part inside
END  DO ! iPart=1,PDM%ParticleVecLength
#endif

#ifdef IMPA
IF(PRESENT(nInnerNewton_In))THEN
  IF(nInnerNewton_In.EQ.-3)THEN
    IPWRITE(UNIt_stdOut,'(I0,A18,I0)') ' nInnerNewton_In ', nInnerNewton_In
  END IF
END IF
#endif /*IMPA*/

END SUBROUTINE ParticleTracing


#ifdef IMPA
SUBROUTINE ParticleRefTracking(doParticle_In)
#else
SUBROUTINE ParticleRefTracking()
#endif
!===================================================================================================================================
! Reference Tracking for particle without treatment of each inner faces
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef,LastPartPos
USE MOD_Mesh_Vars              ,ONLY: OffSetElem,useCurveds,NGeo,ElemBaryNGeo
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Tracking_Vars ,ONLY: nTracks,Distance,ListDistance,CartesianPeriodic
USE MOD_Particle_Mesh_Vars     ,ONLY: Geo,IsTracingBCElem,BCElem,epsOneCell
USE MOD_Utils                  ,ONLY: BubbleSortID,InsertionSort
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps2
USE MOD_Particle_Mesh          ,ONLY: SingleParticleToExactElement,PartInElemCheck
#if USE_MPI
USE MOD_MPI_Vars               ,ONLY: offsetElemMPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartHaloElemToProc
#endif /*USE_MPI*/
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Vars          ,ONLY: PartStateN
USE MOD_TimeDisc_Vars          ,ONLY: iStage
#endif /*IMPA OR ROS*/
#if defined(IMPA)
USE MOD_Particle_Vars          ,ONLY: PartIsImplicit
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: nTracksPerElem
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime, LBElemPauseTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef IMPA
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef IMPA
LOGICAL                           :: doParticle
LOGICAL                           :: doPartInExists
#endif
INTEGER                           :: iPart, ElemID,oldElemID,newElemID
INTEGER                           :: CellX,CellY,CellZ,iBGMElem,nBGMElems
REAL                              :: oldXi(3),newXi(3), LastPos(3),vec(3)!,loc_distance
!REAL                              :: epsOne
#if USE_MPI
INTEGER                           :: InElem
#endif
INTEGER                           :: TestElem,LastElemID
!LOGICAL                           :: ParticleFound(1:PDM%ParticleVecLength),PartisDone
LOGICAL                           :: PartisDone,PartIsMoved
!LOGICAL                           :: HitBC(1:PDM%ParticleVecLength)
REAL                              :: lengthPartTrajectory0, epsElement
! load balance
#if USE_LOADBALANCE
REAL                              :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#ifdef IMPA
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.
#endif /*IMPA*/

DO iPart=1,PDM%ParticleVecLength
#ifdef IMPA
  IF(doPartInExists)THEN
    DoParticle=PDM%ParticleInside(iPart).AND.DoParticle_In(iPart)
  ELSE
    DoParticle=PDM%ParticleInside(iPart)
  END IF
  IF(DoParticle)THEN
#else
  IF (PDM%ParticleInside(iPart)) THEN
#endif /*IMPA*/
    LastElemID = PEM%lastElement(iPart)
    ElemID=LastElemID
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    nTracks=nTracks+1
    ! sanity check
    PartIsDone=.FALSE.
    IF(IsTracingBCElem(ElemID))THEN
      lengthPartTrajectory0=0.
      !IF(GEO%nPeriodicVectors.GT.0.)THEN
      !  lengthPartTrajectory0=BCELEM(ElemID)%ElemToSideDistance(BCElem(ElemID)%lastSide)
      !END IF
      CALL ParticleBCTracking(lengthPartTrajectory0 &
                             ,ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,iPart,PartIsDone,PartIsMoved,1)
      IF(PartIsDone) THEN
#ifdef IMPA
        IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
#endif /*IMPA*/
        CYCLE ! particle has left domain by a boundary condition
      END IF
      IF(PartIsMoved)THEN ! particle is reflected at a wall
        CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
      ELSE
        ! particle has not encountered any boundary condition
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
        CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
#else
#if defined(IMPA) || defined(ROS)
        ! check if particle can be located within the last element
        !loc_distance = ((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
        !               +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
        !               +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )
        !IF(loc_distance.GT.ElemRadius2NGeo(ElemID))THEN
        !  PartPosRef(1:3,iPart) = (/2.,2.,2./)
        !ELSE
          CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
        !END IF
#else
        CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
#endif /*IMPA*/
#endif
      END IF
!      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell(ElemID)) THEN ! particle is inside
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle is inside
         PEM%Element(iPart)=ElemID
#if USE_LOADBALANCE
         IF(ElemID.GT.PP_nElems)THEN
           IF(LastElemID.LE.PP_nElems)THEN
             CALL LBElemPauseTime(LastElemID,tLBStart)
           ELSE
             CALL LBElemPauseTime(ElemID,tLBStart)
           END IF
         END IF
#endif /*USE_LOADBALANCE*/
        CYCLE
      END IF
    ELSE ! no bc elem, therefore, no bc interaction possible
      IF(GEO%nPeriodicVectors.GT.0.AND.CartesianPeriodic)THEN
        ! call here function for mapping of partpos and lastpartpos
        LastPos=PartState(iPart,1:3)
        CALL PeriodicMovement(iPart)
        IF(.NOT.IsTracingBCElem(ElemID))THEN
          DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(iPart,1)) &
              .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(iPart,2)) &
              .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(iPart,3)) )
            LastPos=PartState(iPart,1:3)
            ! call here function for mapping of partpos and lastpartpos
            CALL PeriodicMovement(iPart)
          END DO
        END IF
      END IF
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
      CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
#else
      CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
#endif
      !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle inside
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle inside
        PEM%Element(iPart)  = ElemID
#if USE_LOADBALANCE
         IF(ElemID.LE.PP_nElems)THEN
           CALL LBElemPauseTime(ElemID,tLBStart)
         ELSE IF(PEM%LastElement(iPart).LE.PP_nElems)THEN
           CALL LBElemPauseTime(PEM%LastElement(iPart),tLBStart)
         END IF
#endif /*USE_LOADBALANCE*/
        CYCLE
      !ELSE IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.5) THEN
      !  IPWRITE(UNIT_stdOut,*) ' partposref to large!',iPart
      END IF
    END IF ! initial check
#if USE_LOADBALANCE
    IF(ElemID.LE.PP_nElems)THEN
      CALL LBElemPauseTime(ElemID,tLBStart)
    ELSE IF(PEM%LastElement(iPart).LE.PP_nElems)THEN
      CALL LBElemPauseTime(PEM%LastElement(iPart),tLBStart)
    END IF
#endif /*USE_LOADBALANCE*/
    ! still not located
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    ! relocate particle
    oldElemID = PEM%lastElement(iPart) ! this is not!  a possible elem
    ! get background mesh cell of particle
    CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
    CellX = MAX(MIN(GEO%TFIBGMimax,CellX),GEO%TFIBGMimin)
    CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
    CellY = MAX(MIN(GEO%TFIBGMjmax,CellY),GEO%TFIBGMjmin)
    CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
    CellZ = MAX(MIN(GEO%TFIBGMkmax,CellZ),GEO%TFIBGMkmin)

    ! check all cells associated with this background mesh cell
    nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem
    IF(nBGMElems.GT.1)THEN
      ! get closest element barycenter by looping over all elements in BGMcell
      Distance=-1.
      ListDistance=-1
      DO iBGMElem = 1, nBGMElems
        ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
        ListDistance(iBGMElem)=ElemID
        IF(ElemID.EQ.-1)CYCLE
        IF(ElemID.EQ.OldElemID)THEN
          Distance(iBGMElem)=-1.0
        ELSE
          !Distance(iBGMElem)=SQRT((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))  &
          !                       +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
          !                       +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )
          Distance(iBGMElem)=    ((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
                                 +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
                                 +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )

          IF(Distance(iBGMElem).GT.ElemRadius2NGeo(ElemID))THEN
            Distance(iBGMElem)=-1.0
          END IF
        END IF
      END DO ! nBGMElems

      !CALL BubbleSortID(Distance,ListDistance,nBGMElems)
      CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)
    ELSE IF(nBGMElems.EQ.1)THEN
      Distance(1)=0.
      ListDistance(1)=GEO%TFIBGM(CellX,CellY,CellZ)%Element(1)
    END IF

    OldXi=PartPosRef(1:3,iPart)
    newXi=HUGE(1.0)
    newElemID=-1
    ! loop through sorted list and start by closest element
    DO iBGMElem=1,nBGMElems
      IF(ALMOSTEQUAL(Distance(iBGMELem),-1.0)) CYCLE
      ElemID=ListDistance(iBGMElem)
#if USE_LOADBALANCE
      IF(ElemID.LE.PP_nElems) nTracksPerElem(ElemID)=nTracksPerElem(ElemID)+1
#endif /*USE_LOADBALANCE*/
      CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle inside
      !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle inside
        PEM%Element(iPart) = ElemID
        PartIsDone=.TRUE.
        EXIT
      END IF
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.MAXVAL(ABS(newXi))) THEN
        newXi=PartPosRef(1:3,iPart)
        newElemID=ElemID
      END IF
    END DO ! iBGMElem
    IF(.NOT.PartIsDone)THEN
      ! use best xi
      IF(MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi)))THEN
        PartPosRef(1:3,iPart)=OldXi
        PEM%Element(iPart)=oldElemID
      ELSE
        PartPosRef(1:3,iPart)=NewXi
        PEM%Element(iPart)=NewElemID
        oldElemID=NewElemID
      END IF

      ! set stetelement
      TestElem=PEM%Element(iPart)
      IF(TestElem.EQ.0.)THEN
        epsElement=MAXVAL(epsOneCell)
#if defined(ROS) || defined(IMPA)
        TestElem=PEM%ElementN(iPart)
#else
        TestElem=PEM%Element(iPart)
#endif
      ELSE
        epsElement=epsOneCell(TestElem)
      END IF
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsElement) THEN
        PartIsDone=.FALSE.
        IF(.NOT.IsTracingBCElem(TestElem))THEN
          ! ausgabe
          IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with internal element '
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' xi                     ', PartPosRef(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' epsOneCell             ', epsElement
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' oldxi                  ', oldXi
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' newxi                  ', newXi
          IPWRITE(UNIT_stdOut,'(I0,A)')             ' PartPos:           '
          IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
          IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(iPart,1), GEO%xmaxglob
          IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(iPart,2), GEO%ymaxglob
          IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(iPart,3), GEO%zmaxglob
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos            ', LastPartPos(iPart,1:3)
          Vec=PartState(iPart,1:3)-LastPartPos(iPart,1:3)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
#ifdef IMPA
          IPWRITE(UNIT_stdOut,'(I0,A,X,L)') ' Implicit                ', PartIsImplicit(iPart)
#endif
#if defined(ROS) || defined(IMPA)
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' CurrentStage:    ', iStage
          Vec=PartState(iPart,1:3)-PartStateN(iPart,1:3)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacementN/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartStateN             ', PartStateN(iPart,1:3)
#if USE_MPI
          inelem=PEM%ElementN(ipart)
          IF(inelem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
            IPWRITE(UNIT_stdout,'(I0,A,I0)') ' elemid-N               ', inelem+offsetelem
          ELSE
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid-N              ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
                                                             + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
          END IF
#else
          IPWRITE(UNIt_stdOut,'(I0,A,I0)') ' elemid-N                 ', pem%element(ipart)+offsetelem
#endif
#endif /*ROS or IMPA*/
#if USE_MPI
          InElem=PEM%Element(iPart)
          IF(InElem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID                ', InElem+offSetElem
          ELSE
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                   + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
          END IF
#else
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID       ', PEM%Element(iPart)+offSetElem
#endif
#if USE_MPI
          InElem=PEM%LastElement(iPart)
          IF(InElem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Last-ElemID         ', InElem+offSetElem
          ELSE
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Last-ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                   + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
          END IF
#else
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Last-ElemID  ', PEM%LastElement(iPart)+offSetElem
#endif
CALL abort(&
__STAMP__ &
,'Particle Not inSide of Element, iPart',iPart)
        ELSE ! BCElem
          IPWRITE(UNIT_stdOut,'(I0,A,X,I0)') ' fallback for particle', iPart
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' particlepos            ', partstate(ipart,1:3)
          Vec=PartState(iPart,1:3)-LastPartPos(iPart,1:3)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
          !CALL RefTrackFaceIntersection(ElemID,1,BCElem(ElemID)%nInnerSides,BCElem(ElemID)%nInnerSides,iPart)
          IF(useCurveds)THEN
            IF(NGeo.GT.1)THEN
              CALL FallBackFaceIntersection(TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart)
            END IF
          END IF
          ! no fall back algorithm
          !LastPos=PartState(iPart,1:3)
          lengthPartTrajectory0=0.
          CALL ParticleBCTracking(lengthPartTrajectory0 &
                                 ,TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart,PartIsDone,PartIsMoved,1)
          IF(PartIsDone)THEN
#ifdef IMPA
            IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
#endif /*IMPA*/
            CYCLE
          END IF
          CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),TestElem)
          ! false, reallocate particle
          IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsOneCell(TestElem))THEN
            IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with BC element, relocating!! '
            CALL SingleParticleToExactElement(iPart,doHalo=.TRUE.,initFix=.FALSE.,doRelocate=.FALSE.)
            IF(.NOT.PDM%ParticleInside(iPart)) THEN
              IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with BC element '
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,I0))')    ' iPart                  ', ipart
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' xi                     ', partposref(1:3,ipart)
              IPWRITE(UNIT_stdOut,'(I0,A,1(X,E15.8))') ' epsonecell             ', epsonecell(TestElem)
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' oldxi                  ', oldxi
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' newxi                  ', newxi
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos            ', LastPartPos(iPart,1:3)
              IPWRITE(UNIT_stdOut,'(I0,A)')             ' PartPos:           '
              IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
              IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(iPart,1), GEO%xmaxglob
              IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(iPart,2), GEO%ymaxglob
              IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(iPart,3), GEO%zmaxglob
              Vec=PartState(iPart,1:3)-LastPartPos(iPart,1:3)
              IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
#if defined(ROS) || defined(IMPA)
              IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' CurrentStage:    ', iStage
              Vec=PartState(iPart,1:3)-PartStateN(iPart,1:3)
              IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacementN/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartStateN             ', PartStateN(iPart,1:3)
#if USE_MPI
              inelem=PEM%ElementN(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' elemid-N               ', inelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid-N              ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
                                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' testelem-N            ', testelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' testelem-N             ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
                                                               + PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#else
              IPWRITE(UNIt_stdOut,'(I0,A,I0)') ' elemid-N                 ', pem%element(ipart)+offsetelem
#endif
#endif /*ROS or IMPA*/
#if defined(IMPA)
              IPWRITE(UNIT_stdOut,'(I0,A,X,L)') ' Implicit               ', PartIsImplicit(iPart)
#endif
#if USE_MPI
              inelem=PEM%Element(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' elemid               ', inelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid               ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
                                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' testelem             ', testelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' testelem             ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
                                                               + PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#else
              IPWRITE(UNIt_stdOut,'(I0,A,I0)') ' elemid                 ', pem%element(ipart)+offsetelem
#endif
              CALL abort(&
    __STAMP__ &
    ,'particle noT inside of element, ipart',ipart)
            END IF ! inside
          ELSE
            PEM%Element(iPart)=TestElem
          END IF ! epsCell
        END IF ! BCElem
      END IF ! inner eps to large
    END IF
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
END DO ! iPart

END SUBROUTINE ParticleRefTracking


RECURSIVE SUBROUTINE ParticleBCTracking(lengthPartTrajectory0 &
                                       ,ElemID,firstSide,LastSide,nlocSides,PartId,PartisDone,PartisMoved,iCount)
!===================================================================================================================================
! Calculate intersection with boundary and choose boundary interaction type for reference tracking routine
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Particle_Mesh_Vars,          ONLY:PartBCSideList
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Mesh_Vars,          ONLY:BCElem,GEO,ElemRadiusNGeo
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Intersection,       ONLY:ComputeCurvedIntersection
USE MOD_Particle_Intersection,       ONLY:ComputePlanarRectInterSection
USE MOD_Particle_Intersection,       ONLY:ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection,       ONLY:ComputeBiLinearIntersection
USE MOD_Particle_Tracking_Vars,      ONLY:CartesianPeriodic
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars,      ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,firstSide,LastSide,nlocSides
INTEGER,INTENT(IN)            :: iCount
LOGICAL,INTENT(INOUT)         :: PartisDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES!
INTEGER,INTENT(INOUT)         :: ElemID
LOGICAL,INTENT(INOUT)         :: PartisMoved
REAL,INTENT(INOUT)            :: lengthPartTrajectory0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID, locSideList(firstSide:lastSide), hitlocSide
LOGICAL                       :: ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip,BCSideID,OldElemID
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL                       :: DoTracing,PeriMoved,Reflected
REAL                          :: alphaOld
LOGICAL                       :: doubleCheck
!===================================================================================================================================


PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )

IF(.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(ElemID)))THEN
  PEM%Element(PartID)=ElemID
  PartisDone=.TRUE.
  RETURN
END IF

PartTrajectory=PartTrajectory/lengthPartTrajectory

PartisMoved=.FALSE.
DoTracing=.TRUE.
lengthPartTrajectory0=MAX(lengthPartTrajectory0,lengthPartTrajectory)
! init variables for double check if lastpartpos is close to side and first intersection is found for this position (negative alpha)
doubleCheck = .FALSE.
alphaOld = -1.0

DO WHILE(DoTracing)
  IF(GEO%nPeriodicVectors.GT.0.AND.CartesianPeriodic)THEN
    ! call here function for mapping of partpos and lastpartpos
    CALL PeriodicMovement(PartID,PeriMoved)
    ! the position and trajectory has to be recomputed
    IF(PeriMoved)THEN
      IF(GEO%nPeriodicVectors.EQ.3) CYCLE
      PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
      lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                               +PartTrajectory(2)*PartTrajectory(2) &
                               +PartTrajectory(3)*PartTrajectory(3) )
    ELSE
      IF(GEO%nPeriodicVectors.EQ.3) RETURN
    END IF
  ELSE
    PeriMoved=.FALSE.
  END IF
  locAlpha=-1.0
  nInter=0
  DO iLocSide=firstSide,LastSide
    ! track particle vector until the final particle position is achieved ! check if particle can intersect wit current side
    IF(BCElem(ElemID)%ElemToSideDistance(ilocSide).GT.lengthPartTrajectory0) EXIT
    SideID=BCElem(ElemID)%BCSideID(ilocSide)
    BCSideID=PartBCSideList(SideID)
    locSideList(ilocSide)=ilocSide
    ! get correct flip, wrong for inner sides!!!
    flip  = 0
    !flip  =PartElemToSide(E2S_FLIP,ilocSide,ElemID)
    IF (doublecheck) THEN
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(110("="))')
          WRITE(UNIT_stdout,'(A)')    '     | Particle is double checked: '
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      SELECT CASE(SideType(BCSideID))
      CASE(PLANAR_RECT)
        CALL ComputePlanarRectInterSection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                      ,xi (ilocSide)            &
                                                                                      ,eta(ilocSide)   ,PartID,flip,BCSideID)
      CASE(BILINEAR,PLANAR_NONRECT)
        CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                         ,xi (ilocSide)      &
                                                                                         ,eta(ilocSide)      &
                                                                                         ,PartID,BCSideID &
                                                                                         ,alpha2=alphaOld)
      CASE(PLANAR_CURVED)
        CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                      ,xi (ilocSide)      &
                                                                                      ,eta(ilocSide)   ,PartID,flip,BCSideID)
      CASE(CURVED)
        CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                ,xi (ilocSide)      &
                                                                                ,eta(ilocSide)      ,PartID,BCSideID)
      END SELECT
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (DoubleCheck dorefmapping): '
          WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(BCSideID),' | SideID: ',BCSideID,' | Hit: ',isHit
          WRITE(UNIT_stdout,'(2(A,G0))') '     | Alpha: ',locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
          WRITE(UNIT_stdout,'((A,G0))') '     | AlphaOld: ',alphaOld
          WRITE(UNIT_stdout,'(A,2(X,G0))') '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
    ELSE
      SELECT CASE(SideType(BCSideID))
      CASE(PLANAR_RECT)
        CALL ComputePlanarRectInterSection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                      ,xi (ilocSide)            &
                                                                                      ,eta(ilocSide)   ,PartID,flip,BCSideID)
      CASE(BILINEAR,PLANAR_NONRECT)
        CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                         ,xi (ilocSide)      &
                                                                                         ,eta(ilocSide)      &
                                                                                         ,PartID,BCSideID)
      CASE(PLANAR_CURVED)
        CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                      ,xi (ilocSide)      &
                                                                                      ,eta(ilocSide)   ,PartID,flip,BCSideID)
      CASE(CURVED)
        CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                ,xi (ilocSide)      &
                                                                                ,eta(ilocSide)      ,PartID,BCSideID)
      END SELECT
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (dorefmapping, BCTracing): '
          WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(BCSideID),' | SideID: ',BCSideID,' | Hit: ',isHit
          WRITE(UNIT_stdout,'(2(A,G0))') '     | Alpha: ',locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
          WRITE(UNIT_stdout,'(A,2(X,G0))') '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
    END IF
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      nInter=nInter+1
    END IF
  END DO ! ilocSide

  IF(nInter.EQ.0)THEN
    IF(.NOT.PeriMoved) DoTracing=.FALSE.
  ELSE
    ! take first possible intersection
    !CALL BubbleSortID(locAlpha,locSideList,6)
    PartIsMoved=.TRUE.
    CALL InsertionSort(locAlpha,locSideList,nlocSides)
    DO ilocSide=1,nlocSides
      IF(locAlpha(ilocSide).GT.-1)THEN
        alphaOld = locAlpha(ilocSide)
        EXIT
      END IF
    END DO
    DO ilocSide=1,nlocSides
      IF(locAlpha(ilocSide).GT.-1)THEN
        hitlocSide=locSideList(ilocSide)
        SideID=BCElem(ElemID)%BCSideID(hitlocSide)
        flip  =0 !PartElemToSide(E2S_FLIP,hitlocSide,ElemID) !wrong for inner sides!!!
        BCSideID=PartBCSideList(SideID)
        OldElemID=ElemID
        CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                          ,xi(hitlocSide)     &
                                                                          ,eta(hitlocSide)    &
                                                                          ,PartId,SideID,flip,ElemID,reflected)
        !IF(PEM%Element(PartID).NE.OldElemID)THEN
        IF(ElemID.NE.OldElemID )THEN
          IF (iCount.GE.1000 .AND. MOD(iCount,1000).EQ.0) THEN !threshold might be changed...
            IPWRITE(*,'(I4,A,I0,A,3(x,I0))') ' WARNING: proc has called BCTracking ',iCount  &
              ,'x recursively! Part, Side, Elem:',PartId,SideID,ElemID
          END IF
          IF(GEO%nPeriodicVectors.GT.0)THEN
            lengthPartTrajectory0=BCElem(OldElemID)%ElemToSideDistance(LastSide)
          END IF
          CALL ParticleBCTracking(lengthPartTrajectory0&
                 ,ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,PartID,PartIsDone,PartIsMoved,iCount+1)
          PartisMoved=.TRUE.
          RETURN
        END IF
        IF(reflected) EXIT
      END IF
    END DO
    IF(.NOT.PDM%ParticleInside(PartID)) THEN
      PartisDone = .TRUE.
       RETURN
    END IF
    IF(.NOT.reflected) THEN
      IF (.NOT.doubleCheck) THEN
        doubleCheck = .TRUE.
      ELSE
        DoTracing=.FALSE.
      END IF
    END IF
  END IF ! nInter>0
END DO

END SUBROUTINE ParticleBCTracking


SUBROUTINE SelectInterSectionType(PartIsDone,crossedBC,doLocSide,flip,hitlocSide,ilocSide,PartTrajectory,lengthPartTrajectory &
                                 ,xi,eta,alpha,PartID,SideID,SideType,ElemID,TriNum)
!===================================================================================================================================
! Checks which type of interaction (BC,Periodic,innerSide) has to be applied for the face on the traced particle path
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Tracking_Vars,      ONLY:TriaTracking,TrackInfo
USE MOD_Particle_Surfaces_Vars,      ONLY:SideNormVec
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction,PARTSWITCHELEMENT
USE MOD_Particle_Intersection,       ONLY:IntersectionWithWall
USE MOD_Particle_Vars,               ONLY:PDM
USE MOD_Particle_Surfaces,           ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Mesh_Vars,                   ONLY:BC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: PartID,SideID,hitlocSide,ilocSide,SideType,flip
REAL,INTENT(INOUT)                :: Xi,Eta,Alpha
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(INOUT)             :: PartIsDone
LOGICAL,INTENT(OUT)               :: crossedBC
LOGICAL,INTENT(INOUT)             :: DoLocSide(1:6)
INTEGER,INTENT(INOUT)             :: ElemID
REAL,INTENT(INOUT),DIMENSION(1:3) :: PartTrajectory
REAL,INTENT(INOUT)                :: lengthPartTrajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: Moved(2)
REAL                              :: n_loc(3)
INTEGER                           :: TriNumTemp
!===================================================================================================================================

IF(BC(SideID).GT.0)THEN
  IF (PRESENT(TriNum)) THEN
    TriNumTemp = TriNum
  ELSE
    TriNumTemp = 0
  END IF
  IF (TriaTracking) THEN
    IF (TriNumTemp.NE.1 .AND. TriNumTemp.NE.2) CALL abort(&
__STAMP__ &
,'ERROR in SelectInterSectionType for TriaTracking. TriNum is:',TriNumTemp)
    CALL IntersectionWithWall(PartTrajectory,alpha,PartID,hitlocSide,ElemID,TriNumtemp)
  END IF
  CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha &
                                                                 ,xi    &
                                                                 ,eta   ,PartID,SideID,flip,hitlocSide,ElemID,crossedBC&
                                                                 ,TriNumTemp)
  TrackInfo%CurrElem=ElemID
  IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
  dolocSide=.TRUE.
  !dolocSide(hitlocSide)=.FALSE.
ELSE
  ! DO NOT move particle on edge
  ! issues with periodic grids
  !! move particle ON cell-edge
  !LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+alpha*PartTrajectory(1:3)
  !! recompute remaining particle trajectory
  !lengthPartTrajectory=lengthPartTrajectory-alpha
  ! check if particle leaves element
  IF (.NOT.TriaTracking) THEN
    SELECT CASE(SideType)
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
    IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0) RETURN
  END IF
  ! update particle element
  dolocSide=.TRUE.
  Moved = PARTSWITCHELEMENT(xi,eta,hitlocSide,SideID,ElemID)
  IF (Moved(1).LT.1 .OR. Moved(2).LT.1) CALL abort(&
__STAMP__ &
,'ERROR in SelectInterSectionType. No Neighbour Elem or Neighbour Side found --> increase haloregion')
  ElemID=Moved(1)
  TrackInfo%CurrElem=ElemID
  dolocSide(Moved(2))=.FALSE.
END IF

IF(1.EQ.2)THEN
  moved(1)=ilocSide
END IF

END SUBROUTINE SelectInterSectionType


SUBROUTINE PeriodicMovement(PartID,isMovedOut)
!----------------------------------------------------------------------------------------------------------------------------------!
! move particle in the periodic direction, if particle is outside of the box
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,          ONLY:GEO
USE MOD_Particle_Tracking_Vars,      ONLY:FastPeriodic
#if USE_MPI
USE MOD_Particle_MPI_Vars,           ONLY:PartShiftVector
#endif /*USE_MPI*/
#ifdef IMPA
USE MOD_Particle_Vars,               ONLY: PEM
#endif /*IMPA*/
#ifdef ROS
USE MOD_Particle_Vars,               ONLY: PEM
#endif /*ROS*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: PartID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL   :: isMovedOut
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPV
REAL                            :: MoveVector(1:3)
LOGICAL                         :: isMoved
!===================================================================================================================================

#if USE_MPI
PartShiftVector(1:3,PartID)=PartState(PartID,1:3)
#endif /*USE_MPI*/
isMoved=.FALSE.
IF(FastPeriodic)THEN
  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,1)-GEO%xmaxglob)/ABS(GEO%PeriodicVectors(1,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,1).LT.GEO%xminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,1)-GEO%xminglob)/ABS(GEO%PeriodicVectors(1,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(PartID,2).GT.GEO%ymaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,2)-GEO%ymaxglob)/ABS(GEO%PeriodicVectors(2,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,2).LT.GEO%yminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,2)-GEO%yminglob)/ABS(GEO%PeriodicVectors(2,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(PartID,3).GT.GEO%zmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,3)-GEO%zmaxglob)/ABS(GEO%PeriodicVectors(3,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,3).LT.GEO%zminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,3)-GEO%zminglob)/ABS(GEO%PeriodicVectors(3,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF


  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside x+, PartID',PartID)
    END IF
    IF(PartState(PartID,1).LT.GEO%xminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside x-, PartID',PartID)
    END IF
  END IF
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(PartID,2).GT.GEO%ymaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside y+, PartID',PartID)
    END IF
    IF(PartState(PartID,2).LT.GEO%yminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside y-, PartID',PartID)
    END IF
  END IF
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(PartID,3).GT.GEO%zmaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside z+, PartID',PartID)
    END IF
    IF(PartState(PartID,3).LT.GEO%zminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside z-, PartID',PartID)
    END IF
  END IF
ELSE
  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,1).LT.GEO%xminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF

  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(PartID,2).GT.GEO%ymaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,2).LT.GEO%yminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF

  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(PartID,3).GT.GEO%zmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,3).LT.GEO%zminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF
END IF

#if USE_MPI
PartShiftVector(1:3,PartID)=-PartState(PartID,1:3)+PartShiftvector(1:3,PartID)
#endif /*USE_MPI*/

IF(isMoved)THEN
#if defined(IMPA) || defined(ROS)
  PEM%PeriodicMoved(PartID)=.TRUE.
#endif
END IF

IF(PRESENT(isMovedOut)) isMovedOut=isMoved

END SUBROUTINE PeriodicMovement


SUBROUTINE FallBackFaceIntersection(ElemID,firstSide,LastSide,nlocSides,PartID)
!===================================================================================================================================
! checks if lost particle intersected with face and left Element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Particle_Mesh_Vars,          ONLY:PartBCSideList
USE MOD_Mesh_Vars,                   ONLY:ElemBaryNGeo
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Mesh_Vars,          ONLY:BCElem
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Intersection,       ONLY:ComputeCurvedIntersection
USE MOD_Particle_Intersection,       ONLY:ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection,       ONLY:ComputePlanarRectInterSection
USE MOD_Particle_INtersection,       ONLY:ComputeBiLinearIntersection
USE MOD_Particle_Vars,               ONLY:PartPosRef
USE MOD_Eval_xyz,                    ONLY:TensorProductInterpolation
USE MOD_Mesh_Vars,                   ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID, locSideList(firstSide:lastSide), hitlocSide
LOGICAL                       :: dolocSide(firstSide:lastSide),ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip,BCSideID
REAL                          :: tmpPos(3), tmpLastPartPos(3),tmpVec(3)
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
!===================================================================================================================================

!IPWRITE(*,*) ' Performing fallback algorithm. PartID: ', PartID
tmpPos=PartState(PartID,1:3)
tmpLastPartPos(1:3)=LastPartPos(PartID,1:3)
PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
tmpVec=PartTrajectory

LastPartPos(PartID,1:3)=PartState(PartID,1:3)
!PartState(PartID,1:3)=ElemBaryNGeo(:,ElemID)
LastPartPos(PartID,1:3)=ElemBaryNGeo(:,ElemID)

PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory

locAlpha=-1.0
nInter=0
dolocSide=.TRUE.
!nlocSides=lastSide-firstSide+1
DO iLocSide=firstSide,LastSide
  ! track particle vector until the final particle position is achieved
  SideID=BCElem(ElemID)%BCSideID(ilocSide)
  BCSideID=PartBCSideList(SideID)
  locSideList(ilocSide)=ilocSide
  ! get correct flip, wrong for inner sides!!!
  flip  = 0
  !flip  =PartElemToSide(E2S_FLIP,ilocSide,ElemID)
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT)
    CALL ComputePlanarRectInterSection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                  ,xi (ilocSide)            &
                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)
  CASE(BILINEAR,PLANAR_NONRECT)
    CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                      ,xi (ilocSide)      &
                                                                                      ,eta(ilocSide)      &
                                                                                      ,PartID,BCSideID)
  CASE(PLANAR_CURVED)
    CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                  ,xi (ilocSide)      &
                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)
  CASE(CURVED)
    CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                            ,xi (ilocSide)      &
                                                                            ,eta(ilocSide)      ,PartID,BCSideID)
  END SELECT
  IF(locAlpha(ilocSide).GT.-1.0)THEN
    nInter=nInter+1
  END IF
END DO ! ilocSide

IF(nInter.EQ.0) THEN
  !IPWRITE(*,*) 'not found',PartID
  !IPWRITE(*,*) 'ElemBary',LastPartPos(PartID,1:3)
  !IPWRITE(*,*) 'Part-Pos',tmpPos
  !IPWRITE(*,*) 'LastPart-Pos',tmpLastPartPos
  PartState(PartID,1:3)=tmpPos
  LastPartPos(PartID,1:3)=tmpLastPartPos(1:3)
  IF(PartPosRef(1,PartID).GT. 1.) PartPosRef(1,PartID)= 0.99
  IF(PartPosRef(1,PartID).LT.-1.) PartPosRef(1,PartID)=-0.99
  IF(PartPosRef(2,PartID).GT. 1.) PartPosRef(2,PartID)= 0.99
  IF(PartPosRef(2,PartID).LT.-1.) PartPosRef(2,PartID)=-0.99
  IF(PartPosRef(3,PartID).GT. 1.) PartPosRef(3,PartID)= 0.99
  IF(PartPosRef(3,PartID).LT.-1.) PartPosRef(3,PartID)=-0.99
  CALL TensorProductInterpolation(PartPosRef(:,PartID),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),PartState(PartID,1:3))
  ! crash
  RETURN
ELSE
  ! take first possible intersection
  !CALL BubbleSortID(locAlpha,locSideList,6)
  CALL InsertionSort(locAlpha,locSideList,nlocSides)
  DO ilocSide=firstSide,LastSide
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      hitlocSide=locSideList(ilocSide)
      !SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
      SideID=BCElem(ElemID)%BCSideID(hitlocSide)
      BCSideID=PartBCSideList(SideID)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+0.97*locAlpha(ilocSide)*PartTrajectory
      PartState(PartID,1:3)  =LastPartPos(PartID,1:3)!+tmpVec
      !PartState(PartID,1:3)  =PartState(PartID,1:3)+locAlpha(ilocSide)*PartTrajectory
      !PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
      !lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
      !                         +PartTrajectory(2)*PartTrajectory(2) &
      !                         +PartTrajectory(3)*PartTrajectory(3) )
      !PartTrajectory=PartTrajectory/lengthPartTrajectory
      !locAlpha(ilocSide)=0.
      !CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
      !                                                                  ,xi(hitlocSide)     &
      !                                                                  ,eta(hitlocSide)    &
      !                                                                  ,PartId,SideID)
      !IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
    END IF ! locAlpha>-1.0
  END DO ! ilocSide
END IF ! nInter>0

END SUBROUTINE FallBackFaceIntersection


SUBROUTINE ParticleThroughSideCheck3DFast(PartID,PartTrajectory,iLocSide,Element,ThroughSide,TriNum, IsMortar)
!===================================================================================================================================
!> Routine to check whether a particle crossed the given triangle of a side. The determinant between the normalized trajectory
!> vector and the vectors from two of the three nodes to the old particle position is calculated. If the determinants for the three
!> possible combinations are greater than zero, then the particle went through this triangle of the side.
!> Note that if this is a mortar side, the side of the small neighbouring mortar element has to be checked. Thus, the orientation
!> is reversed.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars, ONLY : GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
REAL,   INTENT(IN)               :: PartTrajectory(1:3)
LOGICAL,INTENT(OUT)              :: ThroughSide
LOGICAL, INTENT(IN), OPTIONAL    :: IsMortar
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: n, NodeID
REAL                             :: Px, Py, Pz
REAL                             :: Vx, Vy, Vz!, Vall
REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)
REAL                             :: det(3)
REAL                             :: eps
!===================================================================================================================================
eps = 0.

ThroughSide = .FALSE.

Px = lastPartPos(PartID,1)
Py = lastPartPos(PartID,2)
Pz = lastPartPos(PartID,3)

! Normalized particle trajectory (PartPos - lastPartPos)/ABS(PartPos - lastPartPos)
Vx = PartTrajectory(1)
Vy = PartTrajectory(2)
Vz = PartTrajectory(3)
! Get the coordinates of the first node and the vector from the particle position to the node
xNode(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
yNode(1) = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
zNode(1) = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
Ax(1) = xNode(1) - Px
Ay(1) = yNode(1) - Py
Az(1) = zNode(1) - Pz
! Get the vectors to the other two nodes, depending on the triangle number
IF(PRESENT(IsMortar)) THEN
  ! Note: reverse orientation in the mortar case, as the side is treated from the perspective of the smaller neighbouring element
  !       (TriNum=1: NodeID=3,2; TriNum=2: NodeID=4,3)
  xNode(2) = GEO%NodeCoords(1,GEO%ElemSideNodeID(2+TriNum,iLocSide,Element))
  yNode(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(2+TriNum,iLocSide,Element))
  zNode(2) = GEO%NodeCoords(3,GEO%ElemSideNodeID(2+TriNum,iLocSide,Element))

  Ax(2) = xNode(2) - Px
  Ay(2) = yNode(2) - Py
  Az(2) = zNode(2) - Pz

  xNode(3) = GEO%NodeCoords(1,GEO%ElemSideNodeID(1+TriNum,iLocSide,Element))
  yNode(3) = GEO%NodeCoords(2,GEO%ElemSideNodeID(1+TriNum,iLocSide,Element))
  zNode(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(1+TriNum,iLocSide,Element))

  Ax(3) = xNode(3) - Px
  Ay(3) = yNode(3) - Py
  Az(3) = zNode(3) - Pz
ELSE
  DO n = 2,3
    NodeID = n+TriNum-1       ! m = true node number of the sides (TriNum=1: NodeID=2,3; TriNum=2: NodeID=3,4)
    xNode(n) = GEO%NodeCoords(1,GEO%ElemSideNodeID(NodeID,iLocSide,Element))
    yNode(n) = GEO%NodeCoords(2,GEO%ElemSideNodeID(NodeID,iLocSide,Element))
    zNode(n) = GEO%NodeCoords(3,GEO%ElemSideNodeID(NodeID,iLocSide,Element))

    Ax(n) = xNode(n) - Px
    Ay(n) = yNode(n) - Py
    Az(n) = zNode(n) - Pz
  END DO
END IF
!--- check whether v and the vectors from the particle to the two edge nodes build
!--- a right-hand-system. If yes for all edges: vector goes potentially through side
det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
          (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
          (Ax(1) * Vy - Ay(1) * Vx) * Az(3))

det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
          (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
          (Ax(2) * Vy - Ay(2) * Vx) * Az(1))

det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
          (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
          (Ax(3) * Vy - Ay(3) * Vx) * Az(2))

! Comparison of the determinants with eps, where a zero is stored (due to machine precision)
IF ((det(1).ge.-eps).AND.(det(2).ge.-eps).AND.(det(3).ge.-eps)) THEN
  ThroughSide = .TRUE.
END IF

RETURN

END SUBROUTINE ParticleThroughSideCheck3DFast


SUBROUTINE ParticleThroughSideLastPosCheck(i,iLocSide,Element,InElementCheck,TriNum,det,isMortarSide,detPartPos)
!===================================================================================================================================
!> Routine used in the case of a particle crossing multipe sides. Calculates the determinant of the three vectors from the last
!> (and current in the case of mortar sides) particle position to the nodes of the triangle in order to determine if the particle
!> is inside the element. Output of determinant is used to determine which of the sides was crossed first.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars,  ONLY : GEO
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: i, Element, iLocSide, TriNum
LOGICAL, INTENT(IN),OPTIONAL     :: isMortarSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: InElementCheck
REAL   ,INTENT(OUT)              :: det
REAL   ,INTENT(OUT), OPTIONAL    :: detPartPos
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NodeNum, ind, iNode
REAL                             :: Ax(3),Ay(3),Az(3)
REAL                             :: NodeCoord(1:3,1:3)
!===================================================================================================================================

InElementCheck = .TRUE.

!--- coords of first node:
DO ind = 1,3
  NodeCoord(ind,1) = GEO%NodeCoords(ind,GEO%ElemSideNodeID(1,iLocSide,Element))
END DO

!--- coords of other two nodes (depending on triangle):
IF(PRESENT(isMortarSide)) THEN
  ! Note: reversed orientation as the triangle is treated from the perspective of the smaller neighbouring mortar element
  NodeCoord(1:3,2)  = GEO%NodeCoords(1:3,GEO%ElemSideNodeID(2+TriNum,iLocSide,Element))
  NodeCoord(1:3,3) = GEO%NodeCoords(1:3,GEO%ElemSideNodeID(1+TriNum,iLocSide,Element))
ELSE
  DO iNode = 2,3
    NodeNum = iNode + TriNum - 1
    DO ind = 1,3
      NodeCoord(ind,iNode) = GEO%NodeCoords(ind,GEO%ElemSideNodeID(NodeNum,iLocSide,Element))
    END DO
  END DO
END IF
!--- vector from lastPos(!) to triangle nodes
DO ind = 1,3
  Ax(ind) = NodeCoord(1,ind) - lastPartPos(i,1)
  Ay(ind) = NodeCoord(2,ind) - lastPartPos(i,2)
  Az(ind) = NodeCoord(3,ind) - lastPartPos(i,3)
END DO

!--- determine whether particle is on inner side (rel. to element) of triangle
!--- set corresponding "flag" (see below)
det = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
       (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
       (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

IF ((det.lt.0).OR.(det.NE.det)) THEN
  InElementCheck = .FALSE.
END IF

IF(PRESENT(isMortarSide).AND.PRESENT(detPartPos)) THEN
  DO ind = 1,3
    Ax(ind) = NodeCoord(1,ind) - PartState(i,1)
    Ay(ind) = NodeCoord(2,ind) - PartState(i,2)
    Az(ind) = NodeCoord(3,ind) - PartState(i,3)
  END DO

  detPartPos = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
                (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
                (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))
END IF

RETURN

END SUBROUTINE ParticleThroughSideLastPosCheck


!SUBROUTINE CheckPlanarInside(PartID,ElemID,lengthPartTrajectory,PartisDone)
!!===================================================================================================================================
!! checks if particle is inside of linear element with planar faces
!!===================================================================================================================================
!! MODULES
!USE MOD_Preproc
!USE MOD_Globals
!USE MOD_Particle_Vars,               ONLY:PartState
!USE MOD_Particle_Surfaces_Vars,      ONLY:SideNormVec,BezierControlPoints3D,epsilontol
!USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToSide,ElemRadiusNGeo
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!INTEGER,INTENT(IN)            :: PartID,ElemID
!REAL,INTENT(IN)               :: lengthPartTrajectory
!LOGICAL,INTENT(INOUT)         :: PartisDone
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                       :: ilocSide, SideID, flip, PlanarSideNum
!REAL                          :: NormVec(1:3), vector_face2particle(1:3), Direction, eps
!!===================================================================================================================================
!PartisDone = .TRUE.
!PlanarSideNum = 0
!eps = ElemRadiusNGeo(ElemID) / lengthPartTrajectory * epsilontol * 1000. !value can be further increased, so far "semi-empirical".
!
!DO ilocSide=1,6
!  SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID)
!  flip =   PartElemToSide(E2S_FLIP,ilocSide,ElemID)
!  ! new with flip
!  IF(flip.EQ.0)THEN
!    NormVec = SideNormVec(1:3,SideID)
!  ELSE
!    NormVec = -SideNormVec(1:3,SideID)
!  END IF
!  vector_face2particle(1:3) = PartState(PartID,1:3) - BezierControlPoints3D(1:3,0,0,SideID)
!  Direction = DOT_PRODUCT(NormVec,vector_face2particle)
!
!  !IF ( (Direction.GE.0.) .OR. (ALMOSTZERO(Direction)) ) THEN
!  IF ( Direction.GE.-eps ) THEN !less rigorous check for planar-assumed sides: they can still be planar-nonrect for which the
!                                !bilin-algorithm will be used which might give a different result for very small distances!
!    PartisDone = .FALSE.
!  END IF
!END DO
!
!END SUBROUTINE CheckPlanarInside


SUBROUTINE ParticleCollectCharges(initialCall_opt)
!===================================================================================================================================
! Routine for communicating and evaluating collected charges on surfaces as floating potential
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars,  ONLY : BCdata_auxSF,SurfMeshSubSideData,SurfFluxSideSize,TriaSurfaceFlux,SideType
USE MOD_Equation_Vars,           ONLY : eps0
USE MOD_TimeDisc_Vars,           ONLY : iter,IterDisplayStep,DoDisplayIter, dt,RKdtFrac
USE MOD_Particle_Mesh_Vars,      ONLY : GEO,NbrOfRegions
USE MOD_Particle_Vars,           ONLY : RegionElectronRef
USE MOD_Mesh_Vars,               ONLY : SideToElem
USE MOD_Particle_Vars,           ONLY : nCollectChargesBCs,CollectCharges
USE MOD_Particle_Boundary_Vars,  ONLY : PartBound
#if USE_MPI
USE MOD_Particle_MPI_Vars,       ONLY : PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL      :: initialCall_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iCC, currentBC
#if USE_MPI
REAL, ALLOCATABLE                :: NewRealChargesGlob(:)
#endif /*USE_MPI*/
REAL, ALLOCATABLE                :: NewRealChargesLoc(:)
INTEGER                          :: iSide, SideID, RegionID, iSample,jSample
REAL                             :: source_e, area
LOGICAL                          :: InitialCall
!===================================================================================================================================
IF (nCollectChargesBCs .GT. 0) THEN
  IF(PRESENT(initialCall_opt))THEN
    InitialCall=initialCall_opt
  ELSE
    InitialCall=.FALSE.
  END IF
  IF (.NOT.InitialCall) THEN !-- normal call during timedisc
    ALLOCATE( NewRealChargesLoc(1:nCollectChargesBCs) )
    NewRealChargesLoc=0.
#if USE_MPI
    ALLOCATE( NewRealChargesGlob(1:nCollectChargesBCs) )
    NewRealChargesGlob=0.
#endif /*USE_MPI*/
    DO iCC=1,nCollectChargesBCs
      NewRealChargesLoc(iCC)=CollectCharges(iCC)%NumOfNewRealCharges
      !--adding BR-electrons if corresponding regions exist: go through sides/trias if present in proc
      IF (NbrOfRegions .GT. 0) THEN
        currentBC = CollectCharges(iCC)%BC
        IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) THEN
          CYCLE
        ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
          CALL abort(__STAMP__,&
            'ERROR in ParticleBoundary: Something is wrong with SideNumber of BC ',currentBC)
        END IF
        DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
          SideID=BCdata_auxSF(currentBC)%SideList(iSide)
          IF (SideToElem(1,SideID).LT.1) CALL abort(&
__STAMP__,&
'ERROR in ParticleBoundary: Something is wrong with SideToElem') !would need SideToElem(2,SideID) as in SurfFlux
          RegionID=GEO%ElemToRegion(SideToElem(1,SideID))
          IF (RegionID .NE. 0) THEN
            source_e = PartBound%Voltage(currentBC)+PartBound%Voltage_CollectCharges(currentBC)-RegionElectronRef(2,RegionID) !dU
            IF (source_e .LT. 0.) THEN
              source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
              * EXP( (source_e) / RegionElectronRef(3,RegionID) )
            ELSE
              source_e = RegionElectronRef(1,RegionID) &         !--- linearized boltzmann relation at positive exponent
              * (1. + ((source_e) / RegionElectronRef(3,RegionID)) )
            END IF
            source_e = source_e*SQRT(RegionElectronRef(3,RegionID))*1.67309567E+05 !rho*v_flux, const.=sqrt(q/(2*pi*m_el))
            area=0.
            IF (TriaSurfaceFlux) THEN
              IF (SideType(SideID).NE.PLANAR_RECT .AND. SideType(SideID).NE.PLANAR_NONRECT) CALL abort(&
__STAMP__&
,'SurfMeshSubSideData has been triangulated for surfaceflux. It can be used only for planar_rect/nonrect!')
            END IF
            DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
              area = area &
                + SurfMeshSubSideData(iSample,jSample,SideID)%area
            END DO; END DO
            NewRealChargesLoc(iCC) = NewRealChargesLoc(iCC) - dt*RKdtFrac*source_e*area !dQ=dt*rho*v_flux*dA
          END IF !RegionID .NE. 0
        END DO ! iSide
      END IF ! NbrOfRegions .GT. 0
      !--
    END DO
#if USE_MPI
    CALL MPI_ALLREDUCE(NewRealChargesLoc,NewRealChargesGlob,nCollectChargesBCs,MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
    DO iCC=1,nCollectChargesBCs
      CollectCharges(iCC)%NumOfNewRealCharges=NewRealChargesGlob(iCC)
    END DO
    DEALLOCATE(NewRealChargesGlob)
#else
    DO iCC=1,nCollectChargesBCs
      CollectCharges(iCC)%NumOfNewRealCharges=NewRealChargesLoc(iCC)
    END DO
#endif /*USE_MPI*/
    DEALLOCATE(NewRealChargesLoc)
    DO iCC=1,nCollectChargesBCs
      CollectCharges(iCC)%NumOfRealCharges=CollectCharges(iCC)%NumOfRealCharges+CollectCharges(iCC)%NumOfNewRealCharges
      CollectCharges(iCC)%NumOfNewRealCharges=0.
      currentBC = CollectCharges(iCC)%BC
      PartBound%Voltage_CollectCharges(currentBC) = ABS(CollectCharges(iCC)%ChargeDist)/eps0 &
        *CollectCharges(iCC)%NumOfRealCharges/BCdata_auxSF(currentBC)%GlobalArea
      IF(BCdata_auxSF(currentBC)%LocalArea.GT.0.) THEN
        IF(DoDisplayIter)THEN
          IF(MOD(iter,IterDisplayStep).EQ.0) THEN
            IPWRITE(*,'(I4,A,I2,2(x,E16.9))') 'Floating Potential and charges',iCC &
              ,PartBound%Voltage_CollectCharges(CollectCharges(iCC)%BC),CollectCharges(iCC)%NumOfRealCharges
          END IF
        END IF
      END IF
    END DO !iCC=1,nCollectChargesBCs
  ELSE ! InitialCall: NumOfRealCharges only (contains initial value)
    DO iCC=1,nCollectChargesBCs
      currentBC = CollectCharges(iCC)%BC
      PartBound%Voltage_CollectCharges(currentBC) = ABS(CollectCharges(iCC)%ChargeDist)/eps0 &
        *CollectCharges(iCC)%NumOfRealCharges/BCdata_auxSF(currentBC)%GlobalArea
      IF(BCdata_auxSF(currentBC)%LocalArea.GT.0.) THEN
        IPWRITE(*,'(I4,A,I2,2(x,E16.9))') 'Initial Phi_float and charges:',iCC &
          ,PartBound%Voltage_CollectCharges(CollectCharges(iCC)%BC),CollectCharges(iCC)%NumOfRealCharges
      END IF
    END DO !iCC=1,nCollectChargesBCs
  END IF
END IF !ChargeCollecting

END SUBROUTINE ParticleCollectCharges

PURE FUNCTION PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo)
!================================================================================================================================
! check if particle has moved significantly within an element
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: lengthPartTrajectory
REAL,INTENT(IN)                      :: ElemRadiusNGeo
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: PARTHASMOVED
!================================================================================================================================

IF(ALMOSTZERO(lengthPartTrajectory/ElemRadiusNGeo))THEN
  PARTHASMOVED=.FALSE.
ELSE
  PARTHASMOVED=.TRUE.
END IF

END FUNCTION PARTHASMOVED


SUBROUTINE ParticleSanityCheck(PartID)
!===================================================================================================================================
! this routine checks the LastPartPos and PartPosition for sanity
! 1) check if LastPartPos is within globalminmax of proc
! 2) check if ParticlePosition is within globalminmax
! 3) check if PartPosRef is within the element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars,          ONLY:PEM,PDM,LastPartPos,PartState
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_TimeDisc_Vars,          ONLY:iStage
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Mesh,          ONLY:PartInElemCheck
#ifdef IMPA
USE MOD_Particle_Vars,          ONLY:PartIsImplicit,PartDtFrac
#endif /*IMPA*/
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY:PartHaloElemToProc
USE MOD_MPI_Vars,               ONLY:offsetElemMPI
#endif /*USE_MPI*/
USE MOD_Mesh_Vars,              ONLY:offsetelem,ElemBaryNGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: ElemID
LOGICAL                          :: IsHit
REAL                             :: IntersectionPoint(1:3)
!===================================================================================================================================

IF(   (LastPartPos(PartID,1).GT.GEO%xmaxglob) &
  .OR.(LastPartPos(PartID,1).LT.GEO%xminglob) &
  .OR.(LastPartPos(PartID,2).GT.GEO%ymaxglob) &
  .OR.(LastPartPos(PartID,2).LT.GEO%yminglob) &
  .OR.(LastPartPos(PartID,3).GT.GEO%zmaxglob) &
  .OR.(LastPartPos(PartID,3).LT.GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(PartID)
#ifdef IMPA
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(PartID)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, LastPartPos(PartID,1), GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, LastPartPos(PartID,2), GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, LastPartPos(PartID,3), GEO%zmaxglob
  CALL abort(&
         __STAMP__ &
         ,' LastPartPos outside of mesh. PartID=, iStage',PartID,REAL(iStage))
END IF
IF(   (PartState(PartID,1).GT.GEO%xmaxglob) &
  .OR.(PartState(PartID,1).LT.GEO%xminglob) &
  .OR.(PartState(PartID,2).GT.GEO%ymaxglob) &
  .OR.(PartState(PartID,2).LT.GEO%yminglob) &
  .OR.(PartState(PartID,3).GT.GEO%zmaxglob) &
  .OR.(PartState(PartID,3).LT.GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(PartID)
#ifdef IMPA
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                        ' PartIsImplicit ', PartIsImplicit(PartID)
      IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                   ' PartDtFrac ', PartDtFrac(PartID)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' LastPartPos    ', LastPartPos(PartID,1:3)
  IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' Velocity       ', PartState(PartID,4:6)
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(PartID,1), GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(PartID,2), GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(PartID,3), GEO%zmaxglob
  CALL abort(&
     __STAMP__ &
     ,' PartPos outside of mesh. PartID=, iStage',PartID,REAL(iStage))
END IF
IF(.NOT.DoRefMapping)THEN
  ElemID=PEM%Element(PartID)
#ifdef CODE_ANALYZE
  CALL PartInElemCheck(PartState(PartID,1:3),PartID,ElemID,isHit,IntersectionPoint,CodeAnalyze_Opt=.TRUE.)
#else
  CALL PartInElemCheck(PartState(PartID,1:3),PartID,ElemID,isHit,IntersectionPoint)
#endif /*CODE_ANALYZE*/
  IF(.NOT.isHit)THEN  ! particle not inside
    IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
    IF(ElemID.LE.PP_nElems)THEN
      IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', ElemID+offSetElem
    ELSE
#if USE_MPI
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,ElemID)) &
                                                    + PartHaloElemToProc(NATIVE_ELEM_ID,ElemID)
#endif /*USE_MPI*/
    END IF
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,ElemID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' IntersectionPoint: ', IntersectionPoint
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos:       ', LastPartPos(PartID,1:3)
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartPos:           ', PartState(PartID,1:3)
    CALL abort(&
    __STAMP__ &
    ,'PartID=. ',PartID)
  END IF
END IF

END SUBROUTINE ParticleSanityCheck

END MODULE MOD_Particle_Tracking
