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
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

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

TYPE,PUBLIC :: tIntersectLink
  REAL                          :: alpha=HUGE(1.)
  REAL                          :: alpha2=HUGE(1.)
  REAL                          :: xi=-1
  REAL                          :: eta=-1
  INTEGER                       :: Side=0
  INTEGER                       :: IntersectCase = 0
  TYPE(tIntersectLink), POINTER :: prev => null()
  TYPE(tIntersectLink), POINTER :: next => null()
END TYPE tIntersectLink
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
USE MOD_Particle_Vars               ,ONLY: PEM,PDM,PartSpecies
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Particle_Mesh_Tools         ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Tracking_vars      ,ONLY: ntracks,MeasureTrackTime,CountNbrOfLostParts, NbrOfLostParticles, TrackInfo
USE MOD_Particle_Tracking_vars      ,ONLY: DisplayLostParticles
USE MOD_Part_Tools                  ,ONLY: StoreLostParticleProperties
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound
USE MOD_Particle_Intersection       ,ONLY: IntersectionWithWall
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_DSMC_Vars                   ,ONLY: RadialWeighting
USE MOD_DSMC_Symmetry2D             ,ONLY: DSMC_2D_RadialWeighting, DSMC_2D_SetInClones
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime, LBElemSplitTime, LBElemPauseTime
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
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
INTEGER                          :: i, NblocSideID, NbElemID, ind, nbSideID, nMortarElems,BCType
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
REAL                             :: PartTrajectory(1:3),lengthPartTrajectory
REAL                             :: xi = -1. , eta = -1. , alpha = -1.
REAL, PARAMETER                  :: eps = 0
!-----------------------------------------------------------------------------------------------------------------------------------
! variabes necessary for tracking with macroparticles inside domain
!LOGICAL                          :: sphereHitExists, sideHitExists
!REAL                             :: alphaPart, alphaPartTria
!REAL                             :: locAlphaPart
!INTEGER                          :: hitSideNbr
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#ifdef IMPA
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.
#endif /*IMPA*/

IF(RadialWeighting%PerformCloning) CALL DSMC_2D_SetInClones()

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
        PartTrajectory=PartState(1:3,i) - LastPartPos(1:3,i)
        lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                                 +PartTrajectory(2)*PartTrajectory(2) &
                                 +PartTrajectory(3)*PartTrajectory(3) )
        IF(ABS(lengthPartTrajectory).GT.0.) PartTrajectory=PartTrajectory/lengthPartTrajectory
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
                CALL ParticleThroughSideCheck3DFast(i,PartTrajectory,NblocSideID,NbElemID,ThroughSide,TriNum,.TRUE.)
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
                CALL ParticleThroughSideCheck3DFast(i,PartTrajectory,localSideID,ElemID,ThroughSide,TriNum)
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
            PDM%ParticleInside(i) = .FALSE.
            IF(CountNbrOfLostParts) THEN
              CALL StoreLostParticleProperties(i, ElemID)
              NbrOfLostParticles=NbrOfLostParticles+1
            END IF
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
              PDM%ParticleInside(i) = .FALSE.
              IF(CountNbrOfLostParts) THEN
                CALL StoreLostParticleProperties(i, ElemID)
                NbrOfLostParticles=NbrOfLostParticles+1
              END IF
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
          IF(BCType.NE.1) CALL IntersectionWithWall(PartTrajectory,alpha,i,LocalSide,ElemID,TriNum)
          CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha &
                                                                       ,xi    &
                                                                       ,eta   ,i,SideID,flip,ElemID,crossedBC&
                                                                       ,TriNum=TriNum,LocSideID=LocalSide)
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
        IF (ElemID.LT.1) THEN
          IPWRITE(UNIT_stdout,*) 'Particle Velocity: ',SQRT(DOTPRODUCT(PartState(4:6,i)))
          CALL abort(&
           __STAMP__ &
           ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
        END IF
        TrackInfo%CurrElem = ElemID
      END IF  ! InElementCheck = T/F
    END DO  ! .NOT.PartisDone
#if USE_LOADBALANCE
    IF (PEM%GlobalElemID(i).GE.offsetElem+1 .AND.PEM%GlobalElemID(i).LE.offsetElem+PP_nElems) &
      CALL LBElemPauseTime(PEM%GlobalElemID(i)-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
  ! Particle treatment for an axisymmetric simulation (cloning/deleting particles)
  IF(RadialWeighting%PerformCloning) THEN
    IF(PDM%ParticleInside(i)) THEN
      IF ((PEM%GlobalElemID(i).GE.1+offSetElem).AND.(PEM%GlobalElemID(i).LE.PP_nElems+offSetElem)) &
          CALL DSMC_2D_RadialWeighting(i,PEM%GlobalElemID(i))
    END IF
  END IF
END DO ! i = 1,PDM%ParticleVecLength

END SUBROUTINE ParticleTriaTracking


#ifdef IMPA
SUBROUTINE ParticleTracing(doParticle_In)
#else
SUBROUTINE ParticleTracing()
#endif /*NOT IMPA*/
!===================================================================================================================================
!> Routine for tracking of moving particles using polynomial description of sides.
!> Routine calculates intersection and boundary interaction for (dorefmapping = false) and (TriaTracking = false)
!> Time is analyzed for LoadBalancing purposes for each element independently because elements with e.g. surface are more costly
!> ---------------------------------------------------------------------------------------------------------------------------------
!> - Loop over all particles, which are in own proc --> PDM%ParticleInside(1:PDM%ParticleVecLength)
!> -- 1. Initialize particle path and tracking info
!> -- 2. Track particle vector up to final particle position
!> -- 3. special check if some double check has to be performed (only necessary for bilinear sides and macrospheres)
!> -- 4. Check if particle intersected a side and also which side (also MacroSpheres and AuxBCs)
!>         For each side only one intersection is chosen, but particle might insersect more than one side. Assign pointer list
!> -- 5. Loop over all intersections in pointer list and check intersection type: inner side, BC, auxBC or MacroSphere
!>       and calculate interaction
!> -- 6. Update particle position and decide if double check might be necessary
!> -- 7. Correct intersection list if double check will be performed and leave loop to do double check
!> -- 8. Reset interscetion list if no double check is performed
!> -- 9. If tolerance was marked, check if particle is inside of proc volume and try to find it in case it was lost
!> ---------------------------------------------------------------------------------------------------------------------------------
!> - DoubleCheck:
!> -- If a tracked particle hits a bilinear side but the PartTrajectory points inside of the element,
!>    then the second alpha for this side might have been the actual intersection, which has been dropped in intersection routine.
!> -- Consequently, alpha for doublecheck side is saved (moved to the last position in intersectionlist)
!>    and neglected during the second check of for the appropriate sideID.
!> -- This occurs after surfaceflux, reflection, or for periodic particles moving almost in tangential direction to bilinear side.
!> -- The DoubleCheck replaces the need of tolerances
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemRadiusNGeo,ElemHasAuxBCs!,PartElemToSide
USE MOD_Particle_Boundary_Vars      ,ONLY: nAuxBCs,UseAuxBCs
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteractionAuxBC
USE MOD_Particle_Tracking_vars      ,ONLY: ntracks, MeasureTrackTime, CountNbrOfLostParts, NbrOfLostParticles, DisplayLostParticles
USE MOD_Particle_Mesh_Tools         ,ONLY: GetGlobalElemID,GetCNElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Localization       ,ONLY: LocateParticleInElement
USE MOD_Particle_Localization       ,ONLY: PartInElemCheck
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectInterSection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeAuxBCIntersection
USE MOD_Eval_xyz                    ,ONLY: GetPositionInRefElem
USE MOD_MacroBody_Vars              ,ONLY: nMacroBody, UseMacroBody, ElemHasMacroBody
USE MOD_MacroBody_tools             ,ONLY: INSIDEMACROBODY
USE MOD_MacroBody_tools             ,ONLY: ComputeMacroSphereIntersection
USE MOD_MacroBody_tools             ,ONLY: GetInteractionWithMacroBody
USE MOD_Part_Tools                  ,ONLY: StoreLostParticleProperties
#ifdef CODE_ANALYZE
#ifdef IMPA
USE MOD_Particle_Vars               ,ONLY: PartIsImplicit,PartDtFrac
USE MOD_Particle_Vars               ,ONLY: PartStateN
#endif /*IMPA*/
USE MOD_Mesh_Vars                   ,ONLY: ElemBaryNGeo
USE MOD_Particle_Intersection       ,ONLY: OutputTrajectory
USE MOD_Particle_Tracking_Vars      ,ONLY: PartOut,MPIRankOut
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO
USE MOD_TimeDisc_Vars               ,ONLY: iStage
#endif /*CODE_ANALYZE*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime,LBElemPauseTime,LBElemSplitTime
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef IMPA
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
#endif /*IMPA*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef IMPA
LOGICAL                       :: doParticle
LOGICAL                       :: doPartInExists
#endif /*IMPA*/
INTEGER                       :: iPart
INTEGER                       :: ElemID,CNElemID,OldElemID,firstElem
INTEGER                       :: ilocSide,SideID,flip,iAuxBC
LOGICAL                       :: dolocSide(1:6)
LOGICAL                       :: PartisDone,foundHit,markTol,crossedBC,SwitchedElement,isCriticalParallelInFace
REAL                          :: localpha,xi,eta
REAL                          :: locAlphaSphere
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL                       :: onlyMacroPart
INTEGER                       :: iMB
REAL                          :: alphaDoneRel, oldLengthPartTrajectory
#if USE_LOADBALANCE
REAL                          :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
LOGICAL                       :: moveList, PartDoubleCheck
#ifdef CODE_ANALYZE
REAL                          :: refpos(1:3)
INTEGER                       :: nIntersections
#endif /*CODE_ANALYZE*/
! intersection info list
TYPE(tIntersectLink),POINTER    :: firstIntersect   => NULL()
TYPE(tIntersectLink),POINTER    :: lastIntersect    => NULL()
TYPE(tIntersectLink),POINTER    :: currentIntersect => NULL()
TYPE(tIntersectLink),POINTER    :: tmp => NULL()
!===================================================================================================================================

#ifdef IMPA
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.
#endif /*IMPA*/

! initialize the first and last pointer in intersection info
IF (.NOT. ASSOCIATED(firstIntersect)) THEN
  ALLOCATE(firstIntersect)
  IF (.NOT. ASSOCIATED(firstIntersect%next)) ALLOCATE(firstIntersect%next)
  lastIntersect => firstIntersect%next
  lastIntersect%prev => firstIntersect
END IF

! IF(RadialWeighting%PerformCloning) CALL DSMC_2D_SetInClones()

DO iPart=1,PDM%ParticleVecLength
  PartDoubleCheck=.FALSE.
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

#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
    ! check if particle is inside domain bounding box
    IF(GEO%nPeriodicVectors.EQ.0)THEN
      IF(   (LastPartPos(1,iPart).GT.GEO%xmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,iPart),GEO%xmaxglob) &
        .OR.(LastPartPos(1,iPart).LT.GEO%xminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,iPart),GEO%xminglob) &
        .OR.(LastPartPos(2,iPart).GT.GEO%ymaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,iPart),GEO%ymaxglob) &
        .OR.(LastPartPos(2,iPart).LT.GEO%yminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,iPart),GEO%yminglob) &
        .OR.(LastPartPos(3,iPart).GT.GEO%zmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,iPart),GEO%zmaxglob) &
        .OR.(LastPartPos(3,iPart).LT.GEO%zminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,iPart),GEO%zminglob) ) THEN
        IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(iPart)
#ifdef IMPA
        IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(iPart)
        IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(iPart)
#endif /*IMPA*/
        IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
        IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
        IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, LastPartPos(1,iPart), GEO%xmaxglob
        IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, LastPartPos(2,iPart), GEO%ymaxglob
        IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, LastPartPos(3,iPart), GEO%zmaxglob
        CALL abort(&
           __STAMP__ &
           ,' LastPartPos outside of mesh. iPart=, iStage',iPart,REAL(iStage))
      END IF
    END IF
    ! caution: reuse of variable, foundHit=TRUE == inside
    ElemID   = PEM%lastElement(iPart)
    CNElemID = GetCNElemID(ElemID)
    CALL GetPositionInRefElem(LastPartPos(1:3,iPart),RefPos,ElemID)
    IF (MAXVAL(ABS(RefPos)).LE.1.0+1e-4) foundHit=.TRUE.
    IF(.NOT.foundHit)THEN  ! particle not inside
     IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')  ' PartID         ', iPart
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')  ' global ElemID  ', ElemID
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,CNElemID)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' LastPartPos:       ', LastPartPos(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' PartPos:           ', PartState(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' PartRefPos:        ', RefPos(1:3)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' Velocity:          ', PartState(4:6,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' PartTrajectory:    ', PartState(1:3,iPart) - LastPartPos(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,ES25.14E3)')      ' lengthPT:          ', SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
     CALL abort(&
         __STAMP__ &
         ,'ERROR: Lastpartpos in wrong element. PartID:',iPart)
    END IF
    IF (INSIDEMACROBODY(LastPartPos(1:3,iPart))) CALL abort(&
        __STAMP__&
        ,'ERROR: particle found inside macroscopic sphere. PartID: ',iPart)
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

! -- 1. Initialize particle path and tracking info
    IF (MeasureTrackTime) nTracks=nTracks+1
    PartisDone=.FALSE.
    ElemID   = PEM%LastGlobalElemID(iPart)
    CNElemID = GetCNElemID(ElemID)

    ! Calculate particle trajectory
    PartTrajectory=PartState(1:3,iPart) - LastPartPos(1:3,iPart)
    lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
    alphaDoneRel=0.
    oldLengthPartTrajectory=LengthPartTrajectory

    ! Check if the particle moved at all. If not, tracking is done
    IF(.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(CNElemID)) .OR. LengthPartTrajectory.EQ.0)THEN
      ! if Macroparticle are in element, they might move and consequently have to be treated although lengthparttrajectory is 0
      ! partvelo - macropartvelo might be > 0 --> Relative lengthPartTrajectory > 0
      IF (UseMacroBody) THEN
        onlyMacroPart=.TRUE.
      ELSE
        PEM%GlobalElemID(iPart)=ElemID
        PartisDone=.TRUE.
        CYCLE
      END IF
    ELSE
      PartTrajectory=PartTrajectory/lengthPartTrajectory
      OnlyMacroPart=.FALSE.
    END IF

#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
    IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A32)')  ' ---------------------------------------------------------------'
      WRITE(UNIT_stdout,'(A)')    '     | Output of Particle information '
      CALL OutputTrajectory(iPart,PartState(1:3,iPart),PartTrajectory,lengthPartTrajectory)
      WRITE(UNIT_stdOut,'(A,I0)') '     | global ElemID       ', PEM%LastGlobalElemID(iPart)
    END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

    ! track particle vector until the final particle position is achieved
    dolocSide=.TRUE.
    firstElem=ElemID
!      !removed CheckPlanarInside since it can be inconsistent for planar-assumed sides:
!      !they can still be planar-nonrect for which the bilin-algorithm will be used which might give a different result
!      !(anyway, this was a speed-up for completely planar meshes only, but those should be now calculated with triatracking)
    markTol =.FALSE.
! -- 2. Track particle vector up to the final particle position
    DO WHILE (.NOT.PartisDone)
      ! do not reset markTol after first intersection of for doublecheck.
      ! This prevents particles to get lost unnoticed in case any intersection has marked tolerance.
      ! markTol =.FALSE.
      IF (PartDoubleCheck) THEN
! -- 3. special check if some double check has to be performed (only necessary for bilinear sides and macrospheres)
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A)')    '     | Calculation of double check: '
        END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        currentIntersect => lastIntersect%prev
        IF (currentIntersect%IntersectCase.EQ.1) THEN
          iLocSide=currentIntersect%Side
          SideID = GetGlobalNonUniqueSideID(ElemID,iLocSide)
          ! TODO: missing!!! : mapping from GlobalNonUnique to CNtotalsides
          CALL ComputeBiLinearIntersection(foundHit,PartTrajectory,lengthPartTrajectory,locAlpha,xi,eta,iPart,SideID &
              ,alpha2=currentIntersect%alpha)
          currentIntersect%alpha=HUGE(1.)
          currentIntersect%IntersectCase=0
          IF(foundHit) THEN
            CALL AssignListPosition(currentIntersect,locAlpha,iLocSide,1,xi_IN=xi,eta_IN=eta)
            IF((ABS(xi).GE.0.99).OR.(ABS(eta).GE.0.99)) markTol=.TRUE.
          END IF
        ELSE IF (currentIntersect%IntersectCase.EQ.3) THEN
          iMB = currentIntersect%Side
          CALL ComputeMacroSphereIntersection(foundHit,PartTrajectory,lengthPartTrajectory,iMB&
              ,locAlpha,locAlphaSphere,alphaDoneRel,iPart,alpha2=currentIntersect%alpha)
          currentIntersect%alpha=HUGE(1.)
          currentIntersect%IntersectCase=0
          IF(foundHit) CALL AssignListPosition(currentIntersect,locAlpha,iMB,3,alpha2_IN=locAlphaSphere)
        END IF
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)')             '     | Output after compute intersection (tracing double check): '
          IF (currentIntersect%IntersectCase.EQ.1) THEN
            WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,' | Hit: ',foundHit
            WRITE(UNIT_stdout,'(A,2(X,G0))')   '     | Intersection xi/eta: ',xi,eta
            WRITE(UNIT_stdout,'((A,G0))')      '     | RelAlpha: ',locAlpha/lengthpartTrajectory
          ELSE IF (currentIntersect%IntersectCase.EQ.3) THEN
            WRITE(UNIT_stdout,'(A,I0,A,L)')    '     | MaroPartID: ',iMB,' | Hit: ',foundHit
            WRITE(UNIT_stdout,'(A,G0)')        '     | AlphaSphere: ',locAlphaSphere
          END IF
          WRITE(UNIT_stdout,'(2(A,G0))')       '     | Alpha: ',locAlpha,' | LengthPartTrajectory: ', lengthPartTrajectory
        END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        ! if double check found no intersection reset entry in list and adjust last entry pointer
        IF (.NOT.foundHit) THEN
          currentIntersect%alpha = HUGE(1.)
          currentIntersect%intersectCase = 0
          IF (ASSOCIATED(currentIntersect%prev) .AND. .NOT.ASSOCIATED(currentIntersect%prev,firstIntersect)) THEN
            lastIntersect => currentIntersect
            lastIntersect%prev => currentIntersect%prev
          END IF
        END IF

      ELSE ! NOT PartDoubleCheck
! -- 4. Check if particle intersected a side and also which side (also MacroSpheres and AuxBCs)
!       For each side only one intersection is chosen, but particle might insersect more than one side. Assign pointer list
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(110("="))')
          WRITE(UNIT_stdout,'(A)')    '     | Calculation of particle intersections: '
        END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        DO ilocSide=1,6
          locAlpha=-1.
          IF (UseMacroBody) THEN
            IF (OnlyMacroPart) CYCLE
          END IF
          IF(.NOT.dolocSide(ilocSide)) CYCLE
          SideID = GetGlobalNonUniqueSideID(ElemID,iLocSide)

          ! BezierControlPoints are now built in cell local system. Hence, sides have always the flip from the shared SideInfo
          flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

          ! TODO missing!!! : mapping from GlobalNonUnique to CNtotalsides
          isCriticalParallelInFace=.FALSE.

          SELECT CASE(SideType(SideID))
            CASE(PLANAR_RECT)
              CALL ComputePlanarRectInterSection(   foundHit,PartTrajectory,lengthPartTrajectory,locAlpha,xi,eta,iPart,flip,SideID  &
                                                ,   isCriticalParallelInFace)
            CASE(BILINEAR,PLANAR_NONRECT)
              CALL ComputeBiLinearIntersection(     foundHit,PartTrajectory,lengthPartTrajectory,locAlpha,xi,eta,iPart,     SideID)
            CASE(PLANAR_CURVED)
              CALL ComputePlanarCurvedIntersection( foundHit,PartTrajectory,lengthPartTrajectory,locAlpha,xi,eta,iPart,flip,SideID  &
                                                  , isCriticalParallelInFace)
            CASE(CURVED)
              CALL ComputeCurvedIntersection(       foundHit,PartTrajectory,lengthPartTrajectory,locAlpha,xi,eta,iPart,     SideID &
                                            ,       isCriticalParallelInFace)
            CASE DEFAULT
              CALL abort(__STAMP__,' Missing required side-data. Please increase halo region. ',SideID)
          END SELECT

#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(30("-"))')
            WRITE(UNIT_stdout,'(A)')           '     | Output after compute intersection (particle tracing): '
            WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,' | Hit: ',foundHit
            WRITE(UNIT_stdout,'(2(A,G0))')     '     | Alpha: ',locAlpha,' | LengthPartTrajectory: ', lengthPartTrajectory
            WRITE(UNIT_stdout,'((A,G0))')      '     | RelAlpha: ',locAlpha/lengthpartTrajectory
            WRITE(UNIT_stdout,'(A,2(X,G0))')   '     | Intersection xi/eta: ',xi,eta
          END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        ! Particle detected inside of face and PartTrajectory parallel to face
          IF(isCriticalParallelInFace)THEN
            IF(DisplayLostParticles)THEN
              IPWRITE(UNIT_stdOut,'(I0,A)') ' Warning: Particle located inside of face and moves parallel to side. Undefined position. '
              IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Removing particle with id: ',iPart
            END IF ! DisplayLostParticles
            PartIsDone=.TRUE.
            PDM%ParticleInside(iPart)=.FALSE.
#ifdef IMPA
            DoParticle=.FALSE.
#endif /*IMPA*/
            IF(CountNbrOfLostParts) THEN
              CALL StoreLostParticleProperties(iPart, ElemID)
              NbrOfLostParticles=NbrOfLostParticles+1
            END IF
            EXIT
          END IF
          IF(foundHit) THEN
            currentIntersect => lastIntersect
            CALL AssignListPosition(currentIntersect,locAlpha,iLocSide,1,xi_IN=xi,eta_IN=eta)
            currentIntersect => lastIntersect
            lastIntersect    => currentIntersect%next
            lastIntersect%prev => currentIntersect
            IF((ABS(xi).GE.0.99).OR.(ABS(eta).GE.0.99)) markTol=.TRUE.
            !IF(ALMOSTZERO(locAlpha)) markTol=.TRUE.
            !IF(locAlpha/lengthPartTrajectory.GE.0.99 .OR. locAlpha/lengthPartTrajectory.LT.0.01) markTol=.TRUE.
          END IF
        END DO ! ilocSide
        IF (UseAuxBCs) THEN
          DO iAuxBC=1,nAuxBCs
            locAlpha=-1
            IF (UseMacroBody) THEN
              IF (OnlyMacroPart) CYCLE
            END IF
            isCriticalParallelInFace=.FALSE.
            IF (ElemHasAuxBCs(ElemID,iAuxBC)) THEN
              CALL ComputeAuxBCIntersection(foundHit,PartTrajectory,lengthPartTrajectory &
                  ,iAuxBC,locAlpha,iPart,isCriticalParallelInFace)
            ELSE
              foundHit=.FALSE.
            END IF
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(30("-"))')
              WRITE(UNIT_stdout,'(A)')        '     | Output after compute AuxBC intersection (particle tracing): '
              WRITE(UNIT_stdout,'(A,I0,A,L)') '     | AuxBC: ',iAuxBC,' | Hit: ',foundHit
              WRITE(UNIT_stdout,'(2(A,G0))')  '     | Alpha: ',locAlpha,' | LengthPartTrajectory: ',lengthPartTrajectory
            END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
            ! Particle detected inside of face and PartTrajectory parallel to face
            IF(isCriticalParallelInFace)THEN
              IF(DisplayLostParticles)THEN
                IPWRITE(UNIT_stdOut,'(I0,A)') ' Warning: Particle located inside of BC and moves parallel to side. Undefined position. '
                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Removing particle with id: ',iPart
              END IF ! DisplayLostParticles
              PartIsDone=.TRUE.
              PDM%ParticleInside(iPart)=.FALSE.
#ifdef IMPA
              DoParticle=.FALSE.
#endif /*IMPA*/
              IF(CountNbrOfLostParts) THEN
                CALL StoreLostParticleProperties(iPart, ElemID)
                NbrOfLostParticles=NbrOfLostParticles+1
              END IF
              EXIT
            END IF
            IF(foundHit) THEN
              ! start from last intersection entry and place current intersection in correct entry position
              currentIntersect => lastIntersect
              CALL AssignListPosition(currentIntersect,locAlpha,iAuxBC,2)
              currentIntersect => lastIntersect
              lastIntersect    => currentIntersect%next
              lastIntersect%prev => currentIntersect
            END IF ! foundHit
          END DO !iAuxBC
        END IF !UseAuxBCs
        IF (UseMacroBody) THEN
          DO iMB=1,nMacroBody
            locAlpha=-1
            IF (ElemHasMacroBody(ElemID,iMB)) THEN
              CALL ComputeMacroSphereIntersection(foundHit,PartTrajectory,lengthPartTrajectory,iMB&
                                               ,locAlpha,locAlphaSphere,alphaDoneRel,iPart)
            ELSE
              foundHit=.FALSE.
            END IF
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(30("-"))')
              WRITE(UNIT_stdout,'(A)')        '     | Output after compute MacroPart intersection (Particle tracing): '
              WRITE(UNIT_stdout,'(A,I0,A,L)') '     | MaroPartID: ',iMB,' | Hit: ',foundHit
              WRITE(UNIT_stdout,'(2(A,G0))')  '     | Alpha: ',locAlpha,' | LengthPartTrajectory: ',lengthPartTrajectory
              WRITE(UNIT_stdout,'(A,G0)')     '     | AlphaSphere: ',locAlphaSphere
            END IF ;END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
            IF(foundHit) THEN
              currentIntersect => lastIntersect
              CALL AssignListPosition(currentIntersect,locAlpha,iMB,3,alpha2_IN=locAlphaSphere)
              currentIntersect => lastIntersect
              lastIntersect => currentIntersect%next
              lastIntersect%prev => currentIntersect
            END IF ! foundHit
          END DO
        END IF
      END IF

! -- 5. Loop over all intersections in pointer list and check intersection type: inner side, BC, auxBC or MacroSphere
!       and calculate interaction
#ifdef CODE_ANALYZE
      nIntersections = 0
#endif /*CODE_ANALYZE*/
      currentIntersect => firstIntersect
      DO WHILE(ASSOCIATED(currentIntersect))
        SwitchedElement=.FALSE.
        crossedBC=.FALSE.
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        nIntersections=nIntersections+1
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(45(":"))')
          WRITE(UNIT_stdout,'(A,I0)')  '     -> Check intersection: ', nIntersections
          WRITE(UNIT_stdout,'(A,I0)')  '     -> Case: '   ,currentIntersect%IntersectCase
          WRITE(UNIT_stdout,'(A,G0)')  '     -> alpha: '  ,currentIntersect%alpha
          WRITE(UNIT_stdout,'(A,I0)')  '     -> locSide: ',currentIntersect%Side
          IF (currentIntersect%IntersectCase.EQ.1) THEN
            WRITE(UNIT_stdout,'(A,I0)') '     -> SideID: ',PartElemToSide(E2S_SIDE_ID,currentIntersect%Side,ElemID)
          END IF
        END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        OldElemID=ElemID
        IF (currentIntersect%IntersectCase.EQ.0) THEN
          ! no intersection
          PEM%GlobalElemID(iPart)=ElemID
          PartisDone=.TRUE.
        ELSE
          SELECT CASE(currentIntersect%IntersectCase)
          !------------------------------------
          CASE(1) ! intersection with cell side
          !------------------------------------
            SideID = GetGlobalNonUniqueSideID(ElemID,currentIntersect%Side)
            flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

            ! missing!!! : mapping from GlobalNonUnique to CNtotalsides
            CALL SelectInterSectionType( PartIsDone                   &
                                       , crossedBC                    &
                                       , doLocSide                    &
                                       , flip                         &
                                       , currentIntersect%Side        &
                                       , PartTrajectory               &
                                       , lengthPartTrajectory         &
                                       , currentIntersect%xi          &
                                       , currentIntersect%eta         &
                                       , currentIntersect%alpha       &
                                       , iPart                        &
                                       , SideID                       &
                                       , SideType(SideID)             &
                                       , ElemID)

            IF (ElemID.NE.OldElemID) THEN
              IF (.NOT.crossedBC) SwitchedElement=.TRUE.
            END IF
          !------------------------------------
          CASE(2) ! AuxBC intersection
          !------------------------------------
            CALL GetBoundaryInteractionAuxBC( PartTrajectory          &
                                            , lengthPartTrajectory    &
                                            , currentIntersect%alpha  &
                                            , iPart                   &
                                            , currentIntersect%Side   &
                                            , crossedBC)
            IF (.NOT.PDM%ParticleInside(iPart)) PartisDone = .TRUE.
            dolocSide=.TRUE. !important when in previously traced portion an elemchange occured, check all sides again!
          !------------------------------------
          CASE(3) ! MacroPart intersection
          !------------------------------------
            CALL GetInteractionWithMacroBody( PartTrajectory           &
                                            , lengthPartTrajectory     &
                                            , currentIntersect%alpha   &
                                            , currentIntersect%alpha2  &
                                            , alphaDoneRel             &
                                            , currentIntersect%Side    &
                                            , iPart                    &
                                            , crossedBC)
            IF(.NOT.PDM%ParticleInside(iPart)) PartisDone = .TRUE.
            dolocSide=.TRUE. !important when in previously traced portion an elemchange occured, check all sides again!
            OnlyMacroPart=.FALSE. !important, since microscopic particle starts to move after collision !
          END SELECT
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
            IF (crossedBC) THEN
              SELECT CASE(currentIntersect%IntersectCase)
              CASE(1) ! intersection with cell side
                WRITE(UNIT_stdout,'(A,L)') '     -> BC was intersected on a side'
              CASE(2) ! AuxBC intersection
                WRITE(UNIT_stdout,'(A,L)') '     -> BC was intersected on an AuxBC'
              CASE(3) ! MacroPart intersection
                WRITE(UNIT_stdout,'(A,L)') '     -> BC was intersected on an MacroPart'
              END SELECT
            END IF
          END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

! -- 6. Update particle position and decide if double check might be necessary
! check what happened with particle (crossedBC or switched element) and set partisdone or double check
#if USE_LOADBALANCE
          IF (OldElemID.GE.offsetElem+1 .AND.OldElemID.LE.offsetElem+PP_nElems) CALL LBElemSplitTime(OldElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
          IF(crossedBC) THEN
            firstElem=ElemID
            IF (UseMacroBody) THEN
              IF (currentIntersect%IntersectCase.EQ.1 .OR. currentIntersect%IntersectCase.EQ.2) THEN
                alphaDoneRel = alphaDoneRel + (1.-alphaDoneRel)*(currentIntersect%alpha/oldLengthPartTrajectory)
              ELSE IF (currentIntersect%IntersectCase.EQ.3) THEN
                alphaDoneRel = alphaDoneRel + (1.-alphaDoneRel)*currentIntersect%alpha2
              END IF
              oldLengthPartTrajectory=LengthPartTrajectory
            END IF
#ifdef IMPA
            IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
#endif /*IMPA*/
          END IF
          IF(crossedBC .OR. SwitchedElement) THEN
            IF (PartDoubleCheck) THEN
              PartDoubleCheck=.FALSE.
            END IF
            EXIT
          ELSE !((.NOT.crossedBC).AND.(.NOT.SwitchedElement)) THEN
            IF (.NOT.PartDoubleCheck) THEN
              PartDoubleCheck = .TRUE.
              PartIsDone= .FALSE.
            END IF
          END IF
        END IF ! IntersectCase.EQ.0

! -- 7. Correct intersection list if double check will be performed and leave loop to do double check
        ! move current list entry to the end and the total list to the front. exit and check if the last is the correct intersection
        IF(.NOT.crossedBC .AND. .NOT.SwitchedElement .AND. .NOT.PartIsDone .AND. PartDoubleCheck) THEN
          moveList=.FALSE.
          SELECT CASE (currentIntersect%intersectCase)
          CASE(1)
            SideID = GetGlobalNonUniqueSideID(OldElemID,currentIntersect%Side)
            ! missing!!! : mapping from GlobalNonUnique to CNtotalsides
            SELECT CASE(SideType(SideID))
            CASE(BILINEAR,PLANAR_NONRECT)
              moveList=.TRUE.
            END SELECT
          CASE(3)
            moveList=.TRUE.
          END SELECT
          IF (moveList) THEN
            lastIntersect%alpha  = currentIntersect%alpha
            lastIntersect%alpha2 = currentIntersect%alpha2
            lastIntersect%xi     = currentIntersect%xi
            lastIntersect%eta    = currentIntersect%eta
            lastIntersect%Side   = currentIntersect%Side
            lastIntersect%intersectCase = currentIntersect%intersectCase
            tmp=>firstIntersect
            DO WHILE (.NOT.ASSOCIATED(tmp,lastIntersect))
              tmp%alpha  = tmp%next%alpha
              tmp%alpha2 = tmp%next%alpha2
              tmp%xi     = tmp%next%xi
              tmp%eta    = tmp%next%eta
              tmp%Side   = tmp%next%Side
              tmp%intersectCase = tmp%next%intersectCase
              tmp=>tmp%next
            END DO
            EXIT
          END IF
        END IF

        ! leave loop because particle is found to remain in element (none of the found intersections is valid)
        currentIntersect=>currentIntersect%next
        IF (ASSOCIATED(currentIntersect,LastIntersect)) THEN
          PartDoubleCheck=.FALSE.
          PartIsDone= .TRUE.
          EXIT
        END IF
      END DO ! ASSOCIATED(currentIntersect)
#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN ; IF(iPart.EQ.PARTOUT)THEN
        WRITE(UNIT_stdout,'(128("="))')
        WRITE(UNIT_stdout,'(A)')           '     | Output of tracking information after the check of number of intersections: '
        WRITE(UNIT_stdout,'(4(A,L))')      '     | crossed Side: ',crossedBC,' switched Element: ',SwitchedElement,&
          ' Particle tracking done: ',PartisDone,' Particle is double checked: ',PartDoubleCheck
        IF(SwitchedElement) THEN
          WRITE(UNIT_stdout,'(A,I0,A,I0)') '     | First_ElemID: ',PEM%LastGlobalElemID(iPart),' | new Element: ',ElemID
          WRITE(UNIT_stdOut,'(A,I0)')      '     | first global ElemID     ', PEM%LastGlobalElemID(iPart)
          WRITE(UNIT_stdOut,'(A,I0)')      '     | new global ElemID       ', ElemID
        END IF
        IF( crossedBC) THEN
          WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Last    PartPos:       ',lastPartPos(1:3,iPart)
          WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Current PartPos:       ',PartState(1:3,iPart)
          WRITE(UNIT_stdout,'(A,3(X,G0))') '     | PartTrajectory:        ',PartTrajectory(1:3)
          WRITE(UNIT_stdout,'(A,(G0))')    '     | Length PartTrajectory: ',lengthPartTrajectory
        END IF
        WRITE(UNIT_stdout,'(128("="))')
      END IF ; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
! -- 8. Reset interscetion list if no double check is performed
      ! reset intersection list because no intersections where found or no double check is performed or no interacions occured
      currentIntersect=>firstIntersect
      IF (currentIntersect%intersectCase.GT.0 .AND. .NOT.PartDoubleCheck)THEN
        DO WHILE (ASSOCIATED(currentIntersect))
          currentIntersect%alpha = HUGE(1.)
          currentIntersect%intersectCase = 0
          IF(ASSOCIATED(currentIntersect,lastIntersect)) THEN
            lastIntersect => firstIntersect%next
            lastIntersect%prev => firstIntersect
            EXIT
          END IF
          currentIntersect => currentIntersect%next
        END DO
      END IF
    END DO ! PartisDone=.FALSE.

#if USE_LOADBALANCE
    IF (PEM%GlobalElemID(iPart).GE.offsetElem+1.AND.PEM%GlobalElemID(iPart).GE.offsetElem+PP_nElems) &
      CALL LBElemPauseTime(PEM%GlobalElemID(iPart),tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF ! Part inside
END DO ! iPart


#ifdef CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
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
    IF( (PartState(1,iPart).GT.GEO%xmaxglob) &
    .OR.(PartState(1,iPart).LT.GEO%xminglob) &
    .OR.(PartState(2,iPart).GT.GEO%ymaxglob) &
    .OR.(PartState(2,iPart).LT.GEO%yminglob) &
    .OR.(PartState(3,iPart).GT.GEO%zmaxglob) &
    .OR.(PartState(3,iPart).LT.GEO%zminglob) ) THEN
#ifdef IMPA
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' DoParticle ', DoParticle
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(iPart)
      IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(iPart)
#endif /*IMPA*/
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
      IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' LastPosition   ', LastPartPos(1:3,iPart)
      IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' Velocity       ', PartState(4:6,iPart)
      IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
      IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(1,iPart), GEO%xmaxglob
      IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(2,iPart), GEO%ymaxglob
      IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(3,iPart), GEO%zmaxglob
      CALL abort(&
     __STAMP__ &
     ,' PartPos outside of mesh AFTER tracking. iPart= ,iStage= ',iPart,REAL(iStage))
    END IF
    ! caution: reuse of variable, foundHit=TRUE == inside
    ElemID=PEM%GlobalElemID(iPart)
    CALL GetPositionInRefElem(PartState(1:3,iPart),RefPos,ElemID)
    IF (MAXVAL(ABS(RefPos)).LE.1.0+1e-4) foundHit=.TRUE.
    IF(.NOT.foundHit)THEN  ! particle not inside
     IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')  ' PartID         ', iPart
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')  ' gloabal ElemID ', ElemID
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,ElemID)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' LastPartPos:       ', LastPartPos(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' PartPos:           ', PartState(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' PartRefPos:        ', RefPos(1:3)
     IPWRITE(UNIT_stdOut,'(I0,A,3(X,ES25.14E3))') ' PartTrajectory:    ', PartTrajectory
     IPWRITE(UNIT_stdOut,'(I0,A,ES25.14E3)')      ' lengthPT:          ', lengthPartTrajectory
     CALL abort(&
       __STAMP__ &
       ,'iPart=. ',iPart)
    END IF

  END IF ! Part inside
END  DO ! iPart=1,PDM%ParticleVecLength
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

END SUBROUTINE ParticleTracing


SUBROUTINE AssignListPosition(inLink,alpha_IN,sideID_IN,IntersectCase_IN,xi_IN,eta_IN,alpha2_IN)
!===================================================================================================================================
!> Checks the given intersection linked list starting from the last to the first entry and compares with a given alpha.
!> adds the given alpha at the correct position extending the list if necessary.
!> exits the search if position was assigned.
!> -----------------------
!> first list entry is the smallest found alpha and last is the largest.
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tIntersectLink),POINTER,INTENT(INOUT) :: inLink
REAL,INTENT(IN)              :: alpha_IN
INTEGER,INTENT(IN)           :: sideID_IN
INTEGER,INTENT(IN)           :: IntersectCase_IN
REAL,INTENT(IN),OPTIONAL     :: xi_IN
REAL,INTENT(IN),OPTIONAL     :: eta_IN
REAL,INTENT(IN),OPTIONAL     :: alpha2_IN
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tIntersectLink),POINTER :: tmpLink
!===================================================================================================================================
! start from last intersection entry and place current intersection in correct entry position
tmpLink => inLink
DO WHILE(ASSOCIATED(tmpLink))
  IF (alpha_IN.LE.tmpLink%alpha) THEN
    IF (.NOT. ASSOCIATED(inLink%next)) THEN
      ALLOCATE(inLink%next)
      inLink%next%prev => inLink
    END IF
    tmpLink%next%alpha = tmpLink%alpha
    tmpLink%next%alpha2 = tmpLink%alpha2
    tmpLink%next%xi = tmpLink%xi
    tmpLink%next%eta = tmpLink%eta
    tmpLink%next%Side = tmpLink%Side
    tmpLink%next%intersectCase = tmpLink%intersectCase
    IF (ASSOCIATED(tmpLink%prev)) THEN
      IF (alpha_IN.GT.tmpLink%prev%alpha) THEN
        ! assign new values
        tmpLink%alpha = alpha_IN
        tmpLink%Side = sideID_IN
        tmpLink%intersectCase = IntersectCase_IN
        IF (PRESENT(xi_IN)) tmpLink%xi = xi_IN
        IF (PRESENT(eta_IN)) tmpLink%eta = eta_IN
        IF (PRESENT(alpha2_IN)) tmpLink%alpha2 = alpha2_IN
        EXIT
      END IF
    ELSE
      ! assign new values
      tmpLink%alpha = alpha_IN
      tmpLink%Side = sideID_IN
      tmpLink%intersectCase = IntersectCase_IN
      IF (PRESENT(xi_IN)) tmpLink%xi = xi_IN
      IF (PRESENT(eta_IN)) tmpLink%eta = eta_IN
      IF (PRESENT(alpha2_IN)) tmpLink%alpha2 = alpha2_IN
      EXIT
    END IF
  END IF
  tmpLink => tmpLink%prev
END DO

END SUBROUTINE AssignListPosition


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
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars              ,ONLY: OffSetElem,useCurveds,NGeo
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_Particle_Localization  ,ONLY: PartInElemCheck
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems,FIBGM_Element,FIBGM_offsetElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps2
USE MOD_Particle_Tracking_Vars ,ONLY: nTracks,Distance,ListDistance,CartesianPeriodic
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef,LastPartPos,PartSpecies
USE MOD_Utils                  ,ONLY: InsertionSort
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
! Counters
INTEGER                           :: iPart
! Elements
INTEGER                           :: ElemID,oldElemID,newElemID,LastElemID
INTEGER                           :: CNElemID
! Background mesh
INTEGER                           :: CellX,CellY,CellZ,iBGMElem,nBGMElems
! Particles
REAL                              :: oldXi(3),newXi(3),LastPos(3)
REAL                              :: vec(3),lengthPartTrajectory0
#if USE_MPI
INTEGER                           :: InElem
#endif
! Tracking
INTEGER                           :: TestElem,CNTestElem
LOGICAL                           :: PartisDone,PartIsMoved
REAL                              :: epsElement
#if USE_LOADBALANCE
REAL                              :: tLBStart
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
    LastElemID = PEM%LastGlobalElemID(iPart)
    ElemID     = LastElemID
    CNElemID   = GetCNElemID(ElemID)
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    nTracks=nTracks+1
    ! sanity check
    PartIsDone=.FALSE.

    ! check if element is a BC element. If yes, handle with Tracing instead of RefMapping
    IF (ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID).GT.0) THEN
      lengthPartTrajectory0 = 0.
      CALL ParticleBCTracking(lengthPartTrajectory0                                                               &
                             ,ElemID                                                                              &
                             ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + 1                                           &
                             ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID)    &
                             ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNElemID)                                               &
                             ,iPart                                                                               &
                             ,PartIsDone                                                                          &
                             ,PartIsMoved                                                                         &
                             ,1)
      IF(PartIsDone) THEN
#ifdef IMPA
        IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
#endif /*IMPA*/
        CYCLE ! particle has left domain by a boundary condition
      END IF

#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
      ! particle has not encountered any boundary condition
      IF (.NOT.PartIsMoved) THEN
        CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)
      ! particle is reflected at a wall
      ELSE
#endif
        CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
      END IF
#endif

      ! particle is inside
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN
         PEM%GlobalElemID(iPart) = ElemID
#if USE_LOADBALANCE
         ! Particle is on current proc, assign load to new cell
         IF (ElemID.GT.offsetElem+1.AND.ElemID.LE.offsetElem+PP_nElems) THEN
           CALL LBElemPauseTime(ElemID-offsetElem,tLBStart)
           END IF
#endif /*USE_LOADBALANCE*/
        CYCLE
      END IF

    ! No boundary element, therefore, no boundary interaction possible
    ELSE
      ! Account for periodic displacement
      IF (GEO%nPeriodicVectors.GT.0.AND.CartesianPeriodic) THEN
        LastPos=PartState(1:3,iPart)
        CALL PeriodicMovement(iPart)
        IF (ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID).EQ.-1) THEN
          DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(1,iPart)) &
              .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(2,iPart)) &
              .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(3,iPart)) )
            LastPos=PartState(1:3,iPart)
            CALL PeriodicMovement(iPart)
          END DO
        END IF
      END IF
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
#else
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)
#endif

      ! Position in reference space smaller unity, particle inside
      IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN
        PEM%GlobalElemID(iPart) = ElemID
#if USE_LOADBALANCE
         ! Particle is on current proc, assign load to new cell
         IF (ElemID.GT.offsetElem+1.AND.ElemID.LE.offsetElem+PP_nElems) THEN
           CALL LBElemPauseTime(ElemID-offsetElem,tLBStart)
         ! Dp not assign load in halo region yet
!         ELSE IF(PEM%LastGlobalElemID(iPart).LE.nComputeNodeTotalElems) THEN
!           CALL LBElemPauseTime(PEM%LastGlobalElemID(iPart),tLBStart)
         END IF
#endif /*USE_LOADBALANCE*/
        CYCLE
      END IF
    END IF ! initial check

#if USE_LOADBALANCE
    ! Particle is on current proc, assign load to new cell
    IF (ElemID.GT.offsetElem+1.AND.ElemID.LE.offsetElem+PP_nElems) THEN
      CALL LBElemPauseTime(ElemID-offsetElem,tLBStart)
    ! Particle moved into halo region, so blame the last cell on proc
    ! Dp not assign load in halo region yet
!    ELSE IF(PEM%LastGlobalElemID(iPart).LE.nComputeNodeTotalElems)THEN
!      CALL LBElemPauseTime(PEM%LastGlobalElemID(iPart),tLBStart)
    END IF
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

    ! relocate particle
    oldElemID = PEM%LastGlobalElemID(iPart) ! this is not!  a possible elem

    ! get background mesh cell of particle
    ! FLOOR might give the wrong element if we are right on the edge
    CellX = MAX(FLOOR((PartState(1,iPart)-GEO%xminglob)/GEO%FIBGMdeltas(1)),0) + 1
    CellY = MAX(FLOOR((PartState(2,iPart)-GEO%yminglob)/GEO%FIBGMdeltas(2)),0) + 1
    CellZ = MAX(FLOOR((PartState(3,iPart)-GEO%zminglob)/GEO%FIBGMdeltas(3)),0) + 1

    ! check all cells associated with this background mesh cell
    nBGMElems = FIBGM_nElems(CellX,CellY,CellZ)

    ! Multiple potential cells found in BGM. Get closest element barycenter by looping over all elements in BGMcell
    IF(nBGMElems.GT.1)THEN
      Distance     = -1.
      ListDistance = -1
      DO iBGMElem = 1,nBGMElems
        ElemID   = FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+iBGMElem)
        CNElemID = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+iBGMElem))
        ListDistance(iBGMElem) = CNElemID

        ! no element associated with BGM elelemt
        IF (ElemID.EQ.-1) &
          CALL ABORT(__STAMP__,'Error during RefMapping: unable to find element associated with BGM element!')

        ! oldElemID was previously checked and particle not found inside
        IF (ElemID.EQ.OldElemID) THEN
          Distance(iBGMElem) = -1.0
        ELSE
          Distance(iBGMElem) = ( (PartState(1,iPart)-ElemBaryNGeo(1,CNElemID))*(PartState(1,iPart)-ElemBaryNGeo(1,CNElemID)) &
                               + (PartState(2,iPart)-ElemBaryNGeo(2,CNElemID))*(PartState(2,iPart)-ElemBaryNGeo(2,CNElemID)) &
                               + (PartState(3,iPart)-ElemBaryNGeo(3,CNElemID))*(PartState(3,iPart)-ElemBaryNGeo(3,CNElemID)))

          ! Do not consider the element if it is too far away
          IF(Distance(iBGMElem).GT.ElemRadius2NGeo(CNElemID))THEN
            Distance(iBGMElem) = -1.0
          END IF
        END IF
      END DO ! nBGMElems
      CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

    ! Only one potential cell found in BGM. No Need to sort
    ELSE IF(nBGMElems.EQ.1)THEN
      Distance(1)     = 0.
      ListDistance(1) = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+1))
    END IF

    OldXi = PartPosRef(1:3,iPart)
    newXi = HUGE(1.0)
    newElemID = -1

    ! loop through sorted list and start by closest element
    DO iBGMElem=1,nBGMElems
      ! ignore old element and elements out of range
      IF (ALMOSTEQUAL(Distance(iBGMELem),-1.0)) CYCLE

      ElemID = GetGlobalElemID(ListDistance(iBGMElem))
#if USE_LOADBALANCE
      ! Cell is on current proc, assign load to new cell
      IF (ElemID.GT.offsetElem+1.AND.ElemID.LE.offsetElem+PP_nElems) THEN
        nTracksPerElem(ElemID-offsetElem)=nTracksPerElem(ElemID-offsetElem)+1
      END IF
#endif /*USE_LOADBALANCE*/
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)

      ! Position in reference space smaller unity, particle inside
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN
        PEM%GlobalElemID(iPart) = ElemID
        PartIsDone         = .TRUE.
        EXIT
      END IF

      ! Position in reference space greater then unity, so particle is not in the current cell. Use it as new guess if it is better than the old guess
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.MAXVAL(ABS(newXi))) THEN
        newXi     = PartPosRef(1:3,iPart)
        newElemID = ElemID
      END IF
    END DO ! iBGMElem

    ! Particle not located, continue by using the best xi
    IF(.NOT.PartIsDone)THEN
      ! Old guess was the best, keep it
      IF (MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi))) THEN
        PartPosRef(1:3,iPart) = OldXi
        PEM%GlobalElemID(iPart)    = oldElemID
      ! New guess is better, so overwrite the old one
      ELSE
        PartPosRef(1:3,iPart) = NewXi
        PEM%GlobalElemID(iPart)    = NewElemID
        oldElemID             = NewElemID
      END IF

      ! Set test element
      TestElem = PEM%GlobalElemID(iPart)
      CNTestElem = GetCNElemID(TestElem)
      IF(CNTestElem.LE.0.)THEN
        epsElement = MAXVAL(ElemEpsOneCell)
#if defined(ROS) || defined(IMPA)
        TestElem = PEM%GlobalElemID(iPart)
        CNTestElem = GetCNElemID(TestElem)
#endif
      ELSE
        epsElement = ElemEpsOneCell(CNTestElem)
      END IF

      ! Position in reference space is outside tolerance of the test element
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsElement) THEN
        PartIsDone=.FALSE.
        ! Element is not a boundary element
        IF (ElemToBCSides(ELEM_NBR_BCSIDES,CNTestElem).LE.0) THEN
          ! ausgabe
          IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with internal element '
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' xi                     ', PartPosRef(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' epsOneCell             ', epsElement
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' oldxi                  ', oldXi
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' newxi                  ', newXi
          IPWRITE(UNIT_stdOut,'(I0,A)')             ' PartPos:           '
          IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
          IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(1,iPart), GEO%xmaxglob
          IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(2,iPart), GEO%ymaxglob
          IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(3,iPart), GEO%zmaxglob
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos            ', LastPartPos(1:3,iPart)
          Vec=PartState(1:3,iPart)-LastPartPos(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
#ifdef IMPA
          IPWRITE(UNIT_stdOut,'(I0,A,X,L)') ' Implicit                ', PartIsImplicit(iPart)
#endif
#if defined(ROS) || defined(IMPA)
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' CurrentStage:    ', iStage
          Vec=PartState(1:3,iPart)-PartStateN(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacementN/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartStateN             ', PartStateN(1:3,iPart)
#if USE_MPI
          InElem=PEM%GlobalElemID(ipart)
          IF(InElem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID                ', InElem+offSetElem
          ELSE
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
!            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid-N              ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
!                                                             + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
          END IF
#else
          IPWRITE(UNIt_stdOut,'(I0,A,I0)') ' elemid-N                 ', PEM%GlobalElemID(ipart)+offsetelem
#endif
#endif /*ROS or IMPA*/
#if USE_MPI
          InElem=PEM%GlobalElemID(iPart)
          IF(InElem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID                ', InElem+offSetElem
          ELSE
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
!            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
!                                                   + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
          END IF
#else
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID       ', PEM%GlobalElemID(iPart)+offSetElem
#endif
#if USE_MPI
          InElem=PEM%LastGlobalElemID(iPart)
          IF(InElem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Last-ElemID         ', InElem+offSetElem
          ELSE
            IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
!            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Last-ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
!                                                   + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
          END IF
#else
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Last-ElemID  ', PEM%LastGlobalElemID(iPart)+offSetElem
#endif
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' PartSpecies  ', PartSpecies(iPart)
          CALL abort(&
              __STAMP__ &
              ,'Particle not inside of Element, ipart',iPart)
        ELSE ! BCElem
          IPWRITE(UNIT_stdOut,'(I0,A,I0,A,I0,A,I0,A)')' Fallback for Particle [', iPart, '] in Element [', TestElem,'] Species [',&
              PartSpecies(iPart),']'
          IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' ParticlePos          ' , partstate(1:3,iPart)
          Vec=PartState(1:3,iPart)-LastPartPos(1:3,iPart)

          ! Tracking on curved meshes with NGeo>1. Check if the particle intersected with face and left the element
          IF(useCurveds)THEN
            IF(NGeo.GT.1)THEN
              CALL FallBackFaceIntersection(TestElem                                                                                &
                                           ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + 1                                             &
                                           ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + ElemToBCSides(ELEM_NBR_BCSIDES,CNTestElem)    &
                                           ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNTestElem)                                                 &
                                           ,iPart)
            END IF
          END IF

          ! No fall back algorithm algorithm, try the normal tracking for boundary elements
          lengthPartTrajectory0 = 0.
          CALL ParticleBCTracking(lengthPartTrajectory0                                                                   &
                                 ,TestElem                                                                                &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + 1                                             &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + ElemToBCSides(ELEM_NBR_BCSIDES,CNTestElem)    &
                                 ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNTestElem)                                                 &
                                 ,iPart                                                                                   &
                                 ,PartIsDone                                                                              &
                                 ,PartIsMoved                                                                             &
                                 ,1)
          IF (PartIsDone) THEN
#ifdef IMPA
            IF(.NOT.PDM%ParticleInside(iPart)) DoParticle=.FALSE.
#endif /*IMPA*/
            CYCLE
          END IF

          CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),TestElem)
          ! false, reallocate particle
          IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.ElemEpsOneCell(CNTestElem))THEN
            IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with BC element, relocating!! '
            CALL LocateParticleInElement(iPart,doHALO=.TRUE.)

            ! Re-localization to a new element failed
            IF(.NOT.PDM%ParticleInside(iPart)) THEN
              IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with BC element '
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,I0))')    ' iPart                  ', ipart
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' xi                     ', partposref(1:3,ipart)
              IPWRITE(UNIT_stdOut,'(I0,A,1(X,E15.8))') ' EpsOneCell             ', ElemEpsOneCell(CNTestElem)
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' oldxi                  ', oldxi
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' newxi                  ', newxi
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos            ', LastPartPos(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A)')             ' PartPos:           '
              IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
              IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(1,iPart), GEO%xmaxglob
              IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(2,iPart), GEO%ymaxglob
              IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(3,iPart), GEO%zmaxglob
              Vec=PartState(1:3,iPart)-LastPartPos(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
#if defined(ROS) || defined(IMPA)
              IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' CurrentStage:    ', iStage
              Vec=PartState(1:3,iPart)-PartStateN(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)') ' displacementN/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
              IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartStateN             ', PartStateN(1:3,iPart)
#if USE_MPI
              inelem=PEM%GlobalElemID(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' elemid-N               ', inelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid-N         ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
!                                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' testelem-N            ', testelem+offsetelem
              ELSE
!                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' testelem-N       ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
                                                               !+ PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#else
              IPWRITE(UNIt_stdOut,'(I0,A,I0)') ' elemid-N                 ', PEM%GlobalElemID(ipart)+offsetelem
#endif
#endif /*ROS or IMPA*/
#if defined(IMPA)
              IPWRITE(UNIT_stdOut,'(I0,A,X,L)') ' Implicit               ', PartIsImplicit(iPart)
#endif
#if USE_MPI
              inelem=PEM%GlobalElemID(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' elemid               ', inelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid         ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
!                                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' testelem             ', testelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' testelem         ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
!                                                               + PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#else
              IPWRITE(UNIt_stdOut,'(I0,A,I0)') ' elemid                 ', PEM%GlobalElemID(ipart)+offsetelem
#endif
              IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' PartSpecies  ', PartSpecies(iPart)
              CALL abort(&
                  __STAMP__ &
                  ,'Particle not inside of Element, ipart',ipart)
            END IF ! inside
          ELSE
            PEM%GlobalElemID(iPart)=TestElem
          END IF ! epsCell
        END IF ! BCElem
      END IF ! inner eps too large
    END IF
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
END DO ! iPart

END SUBROUTINE ParticleRefTracking


RECURSIVE SUBROUTINE ParticleBCTracking(lengthPartTrajectory0 &
                                       ,ElemID,firstSide,lastSide,nlocSides,PartId,PartisDone,PartisMoved,iCount)
!===================================================================================================================================
! Calculate intersection with boundary and choose boundary interaction type for reference tracking routine
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Mesh_Vars          ,ONLY: SideBCMetrics,ElemToBCSides
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO,ElemRadiusNGeo
USE MOD_Particle_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Utils                       ,ONLY: InsertionSort
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectInterSection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Tracking_Vars      ,ONLY: CartesianPeriodic
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars      ,ONLY: PartOut,MPIRankOut
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
! Counters
INTEGER                       :: ilocSide
INTEGER                       :: nInter
! Elements
INTEGER                       :: OldElemID
INTEGER                       :: CNElemID,CNOldElemID
! Sides
INTEGER                       :: SideID,flip
! Particles
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL                       :: DoTracing,PeriMoved,Reflected
! Tracking
INTEGER                       :: locSideList(firstSide:lastSide),hitlocSide
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
REAL                          :: alphaOld
LOGICAL                       :: isHit,doubleCheck
!===================================================================================================================================

! Calculate particle trajectory
PartTrajectory       = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory = SQRT( PartTrajectory(1)*PartTrajectory(1) &
                           + PartTrajectory(2)*PartTrajectory(2) &
                           + PartTrajectory(3)*PartTrajectory(3))

! Check if the particle moved at all. If not, tracking is done
IF (.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(ElemID))) THEN
  PEM%GlobalElemID(PartID) = ElemID
  PartisDone          = .TRUE.
  RETURN
END IF

PartTrajectory=PartTrajectory/lengthPartTrajectory

! Init variables. Double check if LastPartPos is close to side and first intersection is found for this position (negative alpha)
PartisMoved = .FALSE.
DoTracing   = .TRUE.
lengthPartTrajectory0 = MAX(lengthPartTrajectory0,lengthPartTrajectory)
doubleCheck = .FALSE.
alphaOld    = -1.0

DO WHILE(DoTracing)
  IF(GEO%nPeriodicVectors.GT.0.AND.CartesianPeriodic) THEN
    ! Account for periodic displacement
    CALL PeriodicMovement(PartID,PeriMoved)
    ! Position and trajectory has to be recomputed after periodic displacement
    IF(PeriMoved)THEN
      IF(GEO%nPeriodicVectors.EQ.3) CYCLE
      PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
      lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                               +PartTrajectory(2)*PartTrajectory(2) &
                               +PartTrajectory(3)*PartTrajectory(3) )
    ELSE
      IF(GEO%nPeriodicVectors.EQ.3) RETURN
    END IF
  ELSE
    PeriMoved = .FALSE.
  END IF
  locAlpha = -1.0
  nInter   = 0

  ! track particle vector until the final particle position is achieved
  ! check if particle can intersect with current side
  DO iLocSide=firstSide,LastSide

    ! SideBCMetrics is sorted by distance. stop if the first side is out of range
    IF (SideBCMetrics(BCSIDE_DISTANCE,ilocSide).GT.lengthPartTrajectory0) EXIT

    ! side potentially in range (halo_eps)
    SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
    locSideList(ilocSide) = ilocSide

    ! BezierControlPoints are now built in cell local system. Hence, sides have always the flip from the shared SideInfo
    flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

    ! double check
    IF (doublecheck) THEN
#ifdef CODE_ANALYZE
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
        IF (PartID.EQ.PARTOUT) THEN
          WRITE(UNIT_stdout,'(110("="))')
          WRITE(UNIT_stdout,'(A)')    '     | Particle is double checked: '
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      SELECT CASE(SideType(SideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectInterSection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                            ,  xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
        CASE(BILINEAR,PLANAR_NONRECT)
          CALL ComputeBiLinearIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                          ,    xi (ilocSide),eta(ilocSide),PartID,    SideID                &
                                          ,    alpha2=alphaOld)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                          ,    xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(      isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                        ,      xi(ilocSide),eta(ilocSide),PartID,     SideID)
      END SELECT

#ifdef CODE_ANALYZE
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
        IF (PartID.EQ.PARTOUT) THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)')           '     | Output after compute intersection (DoubleCheck dorefmapping): '
          WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,' | Hit: ',isHit
          WRITE(UNIT_stdout,'(2(A,G0))')     '     | Alpha: ',locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
          WRITE(UNIT_stdout,'((A,G0))')      '     | AlphaOld: ',alphaOld
          WRITE(UNIT_stdout,'(A,2(X,G0))')   '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
        END IF
      END IF
#endif /*CODE_ANALYZE*/

    ! not double check
    ELSE
      SELECT CASE(SideType(SideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectInterSection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                            ,  xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
        CASE(BILINEAR,PLANAR_NONRECT)
          CALL ComputeBiLinearIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                          ,    xi(ilocSide),eta(ilocSide),PartID,    SideID)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                              ,xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(      isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                        ,      xi(ilocSide),eta(ilocSide),PartID,     SideID)
      END SELECT

#ifdef CODE_ANALYZE
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
        IF (PartID.EQ.PARTOUT) THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (dorefmapping, BCTracing): '
          WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,' | Hit: ',isHit
          WRITE(UNIT_stdout,'(2(A,G0))') '     | Alpha: ',locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
          WRITE(UNIT_stdout,'(A,2(X,G0))') '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
    END IF

    IF (locAlpha(ilocSide).GT.-1.0) THEN
      nInter = nInter+1
    END IF
  END DO ! ilocSide

  ! No periodic movement in first intersection
  IF (nInter.EQ.0) THEN
    IF(.NOT.PeriMoved) DoTracing=.FALSE.
  ! Take first possible intersection
  ELSE
    PartIsMoved = .TRUE.
    CALL InsertionSort(locAlpha,locSideList,nlocSides)

    DO ilocSide = firstSide,lastSide
      ! Stop on first intersection
      IF (locAlpha(ilocSide).GT.-1) THEN
        alphaOld = locAlpha(ilocSide)
        EXIT
      END IF
    END DO

    DO ilocSide = firstSide,lastSide
      IF (locAlpha(ilocSide).GT.-1) THEN
        hitlocSide = locSideList(ilocSide)
        SideID     = INT(SideBCMetrics(BCSIDE_SIDEID,hitlocSide))
        flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)
        OldElemID = ElemID
        CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                   ,xi(hitlocSide),eta(hitlocSide),PartId,SideID,flip      &
                                   ,ElemID,reflected)

        ! particle moved to a new element in boundary interaction
        IF(ElemID.NE.OldElemID )THEN
          ! Try to recursively calculate the intersection 1000 times. Threshold might be changed...
          IF (iCount.GE.1000 .AND. MOD(iCount,1000).EQ.0) THEN
            IPWRITE(*,'(I4,A,I0,A,3(x,I0))') ' WARNING: proc has called BCTracking ',iCount &
                                             ,'x recursively! Part, Side, Elem:'    ,PartId,SideID,ElemID
          END IF

          ! check if a periodic boundary was crossed during boundary interaction
          CNOldElemID = GetCNElemID(OldElemID)

          IF (GEO%nPeriodicVectors.GT.0) THEN
            lengthPartTrajectory0 = MAXVAL(SideBCMetrics(BCSIDE_DISTANCE,                     &
                                           ElemToBCSides(ELEM_FIRST_BCSIDE,CNOldElemID) + 1:      &
                                           ElemToBCSides(ELEM_FIRST_BCSIDE,CNOldElemID)+ElemToBCSides(ELEM_NBR_BCSIDES,CNOldElemID)))
          END IF

          CNElemID = GetCNElemID(ElemID)
          CALL ParticleBCTracking(lengthPartTrajectory0                                                               &
                                 ,ElemID                                                                              &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + 1                                         &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID)  &
                                 ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNElemID)                                               &
                                 ,PartID                                                                              &
                                 ,PartIsDone                                                                          &
                                 ,PartIsMoved                                                                         &
                                 ,iCount+1)
          PartisMoved = .TRUE.
          RETURN
        END IF

        ! Particle got reflected and stays in the same element. Done with this particle
        IF(reflected) EXIT
      END IF
    END DO ! ilocSide

    ! Particle left the domain through a boundary
    IF(.NOT.PDM%ParticleInside(PartID)) THEN
      PartisDone = .TRUE.
       RETURN
    END IF

    ! Particle did not leave the domain and did not get reflected
    IF (.NOT.reflected) THEN
      ! Double check the particle if not done already
      IF (.NOT.doubleCheck) THEN
        doubleCheck = .TRUE.
      ! Stop tracing if already double checked
      ELSE
        DoTracing=.FALSE.
      END IF
    END IF
  END IF ! nInter>0
END DO

END SUBROUTINE ParticleBCTracking


SUBROUTINE SelectInterSectionType(PartIsDone,crossedBC,doLocSide,flip,hitlocSide,PartTrajectory,lengthPartTrajectory &
                               ,xi,eta,alpha,PartID,SideID,hitSideType,ElemID)
!===================================================================================================================================
!> Use only for TrackingMethod = TRACING or TracingElement with RefMapping
!> Checks which type of interaction (BC,Periodic,innerSide) has to be applied for the face on the traced particle path
!> - If face is BC-side BoundaryInteraction routine is called
!>   - for triatracking the intersection location of partice trajectory with face is calculated first
!> - If face is innerside switch to respective Element
!>   - for tracing check if path for considered intersection point into current element.
!>     Can happen for particles inserted during surface flux at bilinear faces (double checks filters those intersections after)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectInterSection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Tools         ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces           ,ONLY: CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideNormVec
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Tracking_Vars      ,ONLY: TrackInfo
USE MOD_Particle_Vars               ,ONLY: PDM
#if CODE_ANALYZE
USE MOD_Mesh_Vars                   ,ONLY: NGeo
USE MOD_Particle_Localization       ,ONLY: SinglePointToElement
USE MOD_Particle_Surfaces_Vars      ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemBaryNGeo
USE MOD_Particle_Vars               ,ONLY: PartState
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: PartID                   !< Index of Considered Particle
INTEGER,INTENT(IN)                :: SideID                   !< SideID particle intersected with
INTEGER,INTENT(IN)                :: hitlocSide               !< local side of considered element where intersection occured
INTEGER,INTENT(IN)                :: hitSideType              !< type of SideID (planar,bilinear,...)
INTEGER,INTENT(IN)                :: flip                     !< flip of SideID
REAL,INTENT(INOUT)                :: Xi                       !<
REAL,INTENT(INOUT)                :: Eta                      !<
REAL,INTENT(INOUT)                :: Alpha                    !< portion of PartTrajectory until hit with face
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(INOUT)             :: PartIsDone               !< Flag indicating if tracking of PartID is finished
LOGICAL,INTENT(OUT)               :: crossedBC                !< Flag indicating if BC has been hit
LOGICAL,INTENT(INOUT)             :: DoLocSide(1:6)           !<
INTEGER,INTENT(INOUT)             :: ElemID                   !< Element ID particle is currently in
REAL,INTENT(INOUT),DIMENSION(1:3) :: PartTrajectory           !< normalized particle trajectory (x,y,z)
REAL,INTENT(INOUT)                :: lengthPartTrajectory     !< length of particle trajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: isHit
INTEGER                           :: iMortar,nMortarElems
INTEGER                           :: NbElemID,NbSideID
INTEGER                           :: iLocalSide
INTEGER                           :: locFlip
REAL                              :: locAlpha,locXi,locEta
REAL                              :: n_loc(3)
#if CODE_ANALYZE
REAL                              :: v1(3),v2(3)
#endif /* CODE_ANALYZE */
!===================================================================================================================================
! Side is a boundary side
IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
  CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,ElemID,crossedBC &
                                                                 ,locSideID=hitLocSide)
  TrackInfo%CurrElem = ElemID
  IF (.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
  dolocSide = .TRUE.

! Side is not a boundary side
ELSE
  ! get side normal vector. Must be calculated with current point if the side is bilinear or curved
  SELECT CASE(hitSideType)
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc = SideNormVec(1:3,SideID)
      ! BezierControlPoints are now built in cell local system. Hence, side always have the flip from the shared SideInfo
      IF (flip.NE.0) n_loc = -n_loc

    ! bilinear sides have no valid side normal vector. Calculate with current intersection point and do not flip
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    ! curved sides have no valid side normal vector. Calculate with current intersection point and do not flip
    CASE(CURVED)
      CALL CalcNormAndTangBezier(  nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
  END SELECT

#if CODE_ANALYZE
  ! check if normal vector points outwards
  v1 = 0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)  &
           + BezierControlPoints3D(:,NGeo,0   ,SideID)  &
           + BezierControlPoints3D(:,0   ,NGeo,SideID)  &
           + BezierControlPoints3D(:,NGeo,NGeo,SideID))
  v2 = v1  - ElemBaryNGeo(:,ElemID)

  IF (DOT_PRODUCT(v2,n_loc).LT.0) THEN
    IPWRITE(UNIT_stdout,*) 'Obtained wrong side orientation from flip. SideID:',SideID,'flip:',flip,'PartID:',PartID
    IPWRITE(UNIT_stdout,*) 'n_loc (flip)', n_loc,'n_loc (estimated):',v2
    CALL ABORT(__STAMP__,'SideID',SideID)
  END IF
#endif /* CODE_ANALYZE */

  IF (DOT_PRODUCT(n_loc,PartTrajectory).LE.0) RETURN

  ! update particle element
  ! check if the side is a big mortar side
  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

#if CODE_ANALYZE
  WRITE(UNIT_stdout,'(30("-"))')
  WRITE(UNIT_stdout,*) 'ElemID:',ElemID,'PartID',PartID,'SideID:',SideID,'Move rel. to Side:',DOT_PRODUCT(n_loc,PartTrajectory),'NbElemID:',NbElemID, 'PartElem (w/o refl.)', SinglePointToElement(PartState(1:3,PartID),doHalo=.TRUE.)
  WRITE(UNIT_stdout,*) 'PartPos',PartState(1:3,PartID), 'PartVel:',PartState(4:6,PartID)
#endif /* CODE_ANALYZE */

  IF (NbElemID.LT.0) THEN ! Mortar side
  nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

    DO iMortar = 1,nMortarElems
      NbSideID = SideInfo_Shared(SIDE_NBSIDEID,SideID + iMortar)
      ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
      ! no way to recover it during runtime
      IF (NbSideID.LT.1) CALL ABORT(__STAMP__,'Small mortar side not defined!',SideID + iMortar)

      NbElemID = SideInfo_Shared(SIDE_ELEMID,nbSideID)
      ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
      ! no way to recover it during runtime
      IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)

      ! BezierControlPoints are now built in cell local system. We are checking mortar sides, so everything is reversed
      locFlip = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,nbSideID),10),SideInfo_Shared(SIDE_ID,nbSideID).GT.0)

      SELECT CASE(SideType(NbSideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectIntersection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha &
                                            ,  locXi,locEta,PartID,locFlip,NbSideID)
        CASE(BILINEAR,PLANAR_NONRECT)
          CALL ComputeBiLinearIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha &
                                          ,    locXi,locEta,PartID,        NbSideID)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha &
                                          ,    locXi,locEta,PartID,locFlip,NbSideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(      isHit,PartTrajectory,lengthPartTrajectory,locAlpha &
                                        ,      locXi,locEta,PartID,        NbSideID)
      END SELECT

      IF (isHit) THEN
        ElemID = SideInfo_Shared(SIDE_ELEMID,NbSideID)
        TrackInfo%CurrElem = ElemID

        ! flag the side the particle passed through, so it does not have to be checked in the next tracing step
        dolocSide = .TRUE.
        DO iLocalSide = 1,6
          IF (NbSideID.EQ.GetGlobalNonUniqueSideID(NbElemID,iLocalSide)) THEN
            dolocSide(iLocalSide) = .FALSE.
            EXIT
          END IF
        END DO
        RETURN
      END IF
    END DO

    ! passed none of the mortar elements. Keep particle inside current element and warn
    IPWRITE(UNIT_stdOut,*) 'Boundary issue with inner mortar element', ElemID

  ! regular side
  ELSE
    ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
    IF (ElemID.LT.1) &
      CALL abort(__STAMP__,'ERROR in SelectInterSectionType. No Neighbour Elem found!')
!      CALL abort(__STAMP__,'ERROR in SelectInterSectionType. No Neighbour Elem found --> increase haloregion')

    TrackInfo%CurrElem = ElemID

    ! flag the side the particle passed through, so it does not have to be checked in the next tracing step
    dolocSide = .TRUE.
    DO iLocalSide = 1,6
      IF (ABS(SideInfo_Shared(SIDE_ID,SideID)).EQ.ABS(SideInfo_Shared(SIDE_ID,GetGlobalNonUniqueSideID(ElemID,iLocalSide)))) THEN
        dolocSide(iLocalSide) = .FALSE.
        EXIT
      END IF
    END DO

  END IF
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
PartShiftVector(1:3,PartID)=PartState(1:3,PartID)
#endif /*USE_MPI*/

isMoved = .FALSE.

! Routine if particle domain is rectangular. Periodic displacement can be determined by the bounding box
IF(FastPeriodic)THEN
  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(1,PartID).GT.GEO%xmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(1,PartID)-GEO%xmaxglob)/ABS(GEO%PeriodicVectors(1,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(1,PartID).LT.GEO%xminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(1,PartID)-GEO%xminglob)/ABS(GEO%PeriodicVectors(1,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(2,PartID).GT.GEO%ymaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(2,PartID)-GEO%ymaxglob)/ABS(GEO%PeriodicVectors(2,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(2,PartID).LT.GEO%yminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(2,PartID)-GEO%yminglob)/ABS(GEO%PeriodicVectors(2,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(3,PartID).GT.GEO%zmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(3,PartID)-GEO%zmaxglob)/ABS(GEO%PeriodicVectors(3,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(3,PartID).LT.GEO%zminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(3,PartID)-GEO%zminglob)/ABS(GEO%PeriodicVectors(3,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -MoveVector
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF

  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(1,PartID).GT.GEO%xmaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(:,PartID)
      CALL abort(&
      __STAMP__ &
      ,' particle outside x+, PartID',PartID)
    END IF
    IF(PartState(1,PartID).LT.GEO%xminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(:,PartID)
      CALL abort(&
      __STAMP__ &
      ,' particle outside x-, PartID',PartID)
    END IF
  END IF
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(2,PartID).GT.GEO%ymaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(:,PartID)
      CALL abort(&
      __STAMP__ &
      ,' particle outside y+, PartID',PartID)
    END IF
    IF(PartState(2,PartID).LT.GEO%yminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(:,PartID)
      CALL abort(&
      __STAMP__ &
      ,' particle outside y-, PartID',PartID)
    END IF
  END IF
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(3,PartID).GT.GEO%zmaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(:,PartID)
      CALL abort(&
      __STAMP__ &
      ,' particle outside z+, PartID',PartID)
    END IF
    IF(PartState(3,PartID).LT.GEO%zminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(:,PartID)
      CALL abort(&
      __STAMP__ &
      ,' particle outside z-, PartID',PartID)
    END IF
  END IF

! No fast periodic displacement
ELSE
  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(1,PartID).GT.GEO%xmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(1,PartID).LT.GEO%xminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF

  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(2,PartID).GT.GEO%ymaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(2,PartID).LT.GEO%yminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF

  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(3,PartID).GT.GEO%zmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(3,PartID).LT.GEO%zminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(1:3,PartID)  =PartState(1:3,PartID)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(1:3,PartID)  =PartState(1:3,PartID)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(1:3,PartID)=LastPartPos(1:3,PartID)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF
END IF

#if USE_MPI
PartShiftVector(1:3,PartID)=-PartState(1:3,PartID)+PartShiftvector(1:3,PartID)
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
! Checks if lost particle intersected with face and left Element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Eval_xyz,                    ONLY: TensorProductInterpolation
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
USE MOD_Particle_Localization,       ONLY:LocateParticleInElement
USE MOD_Particle_Intersection,       ONLY:ComputeCurvedIntersection
USE MOD_Particle_Intersection,       ONLY:ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection,       ONLY:ComputePlanarRectInterSection
USE MOD_Particle_Intersection,       ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Vars,          ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars,          ONLY: SideBCMetrics
USE MOD_Particle_Mesh_Vars,          ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Tools,         ONLY: GetCNElemID
USE MOD_Particle_Surfaces_Vars,      ONLY: SideType
USE MOD_Utils                       ,ONLY: InsertionSort
USE MOD_Particle_Vars,               ONLY: PDM,PartState,LastPartPos
USE MOD_Particle_Vars,               ONLY:PartPosRef
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID, locSideList(firstSide:lastSide), hitlocSide
LOGICAL                       :: dolocSide(firstSide:lastSide)
LOGICAL                       :: ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip
REAL                          :: tmpPos(3), tmpLastPartPos(3),tmpVec(3)
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
!===================================================================================================================================

! Get particle position and trajectory
tmpPos              = PartState(1:3,PartID)
tmpLastPartPos(1:3) = LastPartPos(1:3,PartID)
tmpVec              = PartTrajectory
LastPartPos(1:3,PartID) = ElemBaryNGeo(:,GetCNElemID(ElemID))

PartTrajectory = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
IF (lengthPartTrajectory.GT.0) PartTrajectory = PartTrajectory/lengthPartTrajectory

locAlpha  = -1.0
nInter    = 0
dolocSide =.TRUE.

! Loop over all sides
DO iLocSide=firstSide,LastSide
  ! track particle vector until the final particle position is achieved
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
  locSideList(ilocSide) = ilocSide
  flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectInterSection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                        ,  xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
    CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                      ,    xi (ilocSide),eta(ilocSide),PartID,    SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                      ,    xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(      isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                    ,      xi(ilocSide),eta(ilocSide),PartID,     SideID)
  END SELECT

  IF (locAlpha(ilocSide).GT.-1.0) THEN
    nInter=nInter+1
  END IF
END DO ! ilocSide

! no intersection found. Try to locate particle manually
IF(nInter.EQ.0) THEN
  PartState(1:3,PartID)   = tmpPos
  LastPartPos(1:3,PartID) = tmpLastPartPos(1:3)
  IF(PartPosRef(1,PartID).GT. 1.) PartPosRef(1,PartID)= 0.99
  IF(PartPosRef(1,PartID).LT.-1.) PartPosRef(1,PartID)=-0.99
  IF(PartPosRef(2,PartID).GT. 1.) PartPosRef(2,PartID)= 0.99
  IF(PartPosRef(2,PartID).LT.-1.) PartPosRef(2,PartID)=-0.99
  IF(PartPosRef(3,PartID).GT. 1.) PartPosRef(3,PartID)= 0.99
  IF(PartPosRef(3,PartID).LT.-1.) PartPosRef(3,PartID)=-0.99
  CALL LocateParticleInElement(PartID,doHalo=.FALSE.)

  ! particle successfully located
  IF (PDM%ParticleInside(PartID)) THEN
    RETURN
  ELSE
    CALL ABORT(__STAMP__,'FallBackFaceIntersection failed for particle!')
  END IF

! One or more intersection
ELSE
  ! take first possible intersection and place particle "just a little" further back
  CALL InsertionSort(locAlpha,locSideList,nlocSides)
  DO ilocSide=firstSide,LastSide
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      hitlocSide = locSideList(ilocSide)
      SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
      LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID)+0.97*locAlpha(ilocSide)*PartTrajectory
      PartState  (1:3,PartID) = LastPartPos(1:3,PartID)
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
USE MOD_Particle_Mesh_Tools,         ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
REAL,   INTENT(IN)               :: PartTrajectory(1:3)
LOGICAL, INTENT(IN), OPTIONAL    :: IsMortar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: ThroughSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
INTEGER                          :: n, NodeID
REAL                             :: Px, Py, Pz
REAL                             :: Vx, Vy, Vz
REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)
REAL                             :: det(3)
REAL                             :: eps
!===================================================================================================================================

eps = 0.
CNElemID = GetCNElemID(Element)

ThroughSide = .FALSE.

Px = lastPartPos(1,PartID)
Py = lastPartPos(2,PartID)
Pz = lastPartPos(3,PartID)

! Normalized particle trajectory (PartPos - lastPartPos)/ABS(PartPos - lastPartPos)
Vx = PartTrajectory(1)
Vy = PartTrajectory(2)
Vz = PartTrajectory(3)

! Get the coordinates of the first node and the vector from the particle position to the node
xNode(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
yNode(1) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
zNode(1) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)

Ax(1) = xNode(1) - Px
Ay(1) = yNode(1) - Py
Az(1) = zNode(1) - Pz
! Get the vectors to the other two nodes, depending on the triangle number
IF(PRESENT(IsMortar)) THEN
  ! Note: reverse orientation in the mortar case, as the side is treated from the perspective of the smaller neighbouring element
  !       (TriNum=1: NodeID=3,2; TriNum=2: NodeID=4,3)
  xNode(2) = NodeCoords_Shared(1,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  yNode(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  zNode(2) = NodeCoords_Shared(3,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)

  Ax(2) = xNode(2) - Px
  Ay(2) = yNode(2) - Py
  Az(2) = zNode(2) - Pz

  xNode(3) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
  yNode(3) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
  zNode(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)

  Ax(3) = xNode(3) - Px
  Ay(3) = yNode(3) - Py
  Az(3) = zNode(3) - Pz
ELSE
  DO n = 2,3
    NodeID = n+TriNum-1       ! m = true node number of the sides (TriNum=1: NodeID=2,3; TriNum=2: NodeID=3,4)
    xNode(n) = NodeCoords_Shared(1,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
    yNode(n) = NodeCoords_Shared(2,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
    zNode(n) = NodeCoords_Shared(3,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)

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

! Comparison of the determinants with against machine precision (eps)
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
USE MOD_Particle_Mesh_Tools,         ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Vars
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
INTEGER                          :: CNElemID
INTEGER                          :: NodeNum, ind, iNode
REAL                             :: Ax(3),Ay(3),Az(3)
REAL                             :: NodeCoord(1:3,1:3)
!===================================================================================================================================

CNElemID = GetCNElemID(Element)

InElementCheck = .TRUE.

!--- coords of first node:
DO ind = 1,3
  NodeCoord(ind,1) = NodeCoords_Shared(ind,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
END DO

!--- coords of other two nodes (depending on triangle):
IF(PRESENT(isMortarSide)) THEN
  ! Note: reversed orientation as the triangle is treated from the perspective of the smaller neighbouring mortar element
  NodeCoord(1:3,2) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  NodeCoord(1:3,3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
ELSE
  DO iNode = 2,3
    NodeNum = iNode + TriNum - 1
    DO ind = 1,3
      NodeCoord(ind,iNode) = NodeCoords_Shared(ind,ElemSideNodeID_Shared(NodeNum,iLocSide,CNElemID)+1)
    END DO
  END DO
END IF
!--- vector from lastPos(!) to triangle nodes
DO ind = 1,3
  Ax(ind) = NodeCoord(1,ind) - lastPartPos(1,i)
  Ay(ind) = NodeCoord(2,ind) - lastPartPos(2,i)
  Az(ind) = NodeCoord(3,ind) - lastPartPos(3,i)
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
    Ax(ind) = NodeCoord(1,ind) - PartState(1,i)
    Ay(ind) = NodeCoord(2,ind) - PartState(2,i)
    Az(ind) = NodeCoord(3,ind) - PartState(3,i)
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
!  vector_face2particle(1:3) = PartState(1:3,PartID) - BezierControlPoints3D(1:3,0,0,SideID)
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
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF,SurfMeshSubSideData,SurfFluxSideSize,TriaSurfaceFlux,SideType
USE MOD_Globals_Vars           ,ONLY: eps0
USE MOD_TimeDisc_Vars          ,ONLY: iter,IterDisplayStep,DoDisplayIter, dt,RKdtFrac
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,NbrOfRegions
USE MOD_Particle_Vars          ,ONLY: RegionElectronRef
USE MOD_Mesh_Vars              ,ONLY: SideToElem
USE MOD_Particle_Vars          ,ONLY: nCollectChargesBCs,CollectCharges
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
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
USE MOD_Mesh_Vars,              ONLY:offsetelem
USE MOD_Particle_Localization,  ONLY:PartInElemCheck
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNGeo_Shared
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Vars,          ONLY:PEM,PDM,LastPartPos,PartState
USE MOD_TimeDisc_Vars,          ONLY:iStage
#ifdef IMPA
USE MOD_Particle_Vars,          ONLY:PartIsImplicit,PartDtFrac
#endif /*IMPA*/
#if USE_MPI
USE MOD_MPI_Vars,               ONLY:offsetElemMPI
#endif /*USE_MPI*/
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

IF(   (LastPartPos(1,PartID).GT.GEO%xmaxglob) &
  .OR.(LastPartPos(1,PartID).LT.GEO%xminglob) &
  .OR.(LastPartPos(2,PartID).GT.GEO%ymaxglob) &
  .OR.(LastPartPos(2,PartID).LT.GEO%yminglob) &
  .OR.(LastPartPos(3,PartID).GT.GEO%zmaxglob) &
  .OR.(LastPartPos(3,PartID).LT.GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(PartID)
#ifdef IMPA
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(PartID)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, LastPartPos(1,PartID), GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, LastPartPos(2,PartID), GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, LastPartPos(3,PartID), GEO%zmaxglob
  CALL abort(&
         __STAMP__ &
         ,' LastPartPos outside of mesh. PartID=, iStage',PartID,REAL(iStage))
END IF
IF(   (PartState(1,PartID).GT.GEO%xmaxglob) &
  .OR.(PartState(1,PartID).LT.GEO%xminglob) &
  .OR.(PartState(2,PartID).GT.GEO%ymaxglob) &
  .OR.(PartState(2,PartID).LT.GEO%yminglob) &
  .OR.(PartState(3,PartID).GT.GEO%zmaxglob) &
  .OR.(PartState(3,PartID).LT.GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(PartID)
#ifdef IMPA
      IPWRITE(UNIt_stdOut,'(I0,A18,L)')                        ' PartIsImplicit ', PartIsImplicit(PartID)
      IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                   ' PartDtFrac ', PartDtFrac(PartID)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' LastPartPos    ', LastPartPos(1:3,PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,3(X,E27.16))')                  ' Velocity       ', PartState(4:6,PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' x', GEO%xminglob, PartState(1,PartID), GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' y', GEO%yminglob, PartState(2,PartID), GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,E27.16,x,E27.16,x,E27.16)') ' z', GEO%zminglob, PartState(3,PartID), GEO%zmaxglob
  CALL abort(&
     __STAMP__ &
     ,' PartPos outside of mesh. PartID=, iStage',PartID,REAL(iStage))
END IF
IF(.NOT.DoRefMapping)THEN
  ElemID=PEM%GlobalElemID(PartID)
#ifdef CODE_ANALYZE
  CALL PartInElemCheck(PartState(1:3,PartID),PartID,ElemID,isHit,IntersectionPoint,CodeAnalyze_Opt=.TRUE.)
#else
  CALL PartInElemCheck(PartState(1:3,PartID),PartID,ElemID,isHit,IntersectionPoint)
#endif /*CODE_ANALYZE*/
  IF(.NOT.isHit)THEN  ! particle not inside
    IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
    IF(ElemID.LE.PP_nElems)THEN
      IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', ElemID+offSetElem
!    ELSE
!#if USE_MPI
!          IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' ElemID         ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,ElemID)) &
!                                                    + PartHaloElemToProc(NATIVE_ELEM_ID,ElemID)
!#endif /*USE_MPI*/
    END IF
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' ElemBaryNGeo:      ', ElemBaryNGeo_Shared(1:3,ElemID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' IntersectionPoint: ', IntersectionPoint
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' LastPartPos:       ', LastPartPos(1:3,PartID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartPos:           ', PartState(1:3,PartID)
    CALL abort(&
    __STAMP__ &
    ,'PartID=. ',PartID)
  END IF
END IF

END SUBROUTINE ParticleSanityCheck

END MODULE MOD_Particle_Tracking
