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

MODULE MOD_Particle_RefTracking
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC::ParticleRefTracking
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

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
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID,GetCNElemID
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

!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
!      ! particle has not encountered any boundary condition
!      IF (.NOT.PartIsMoved) THEN
!        CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)
!      ! particle is reflected at a wall
!      ELSE
!#endif
        CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)
!      END IF
!#endif

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
    ! get background mesh cell of particle. MAX might give the wrong element if we are right on the edge, so restrict to valid value
    CellX = MAX(FLOOR((PartState(1,iPart)-GEO%xminglob)/GEO%FIBGMdeltas(1)),0) + 1
    CellX = MIN(GEO%FIBGMimax,CellX)
    CellY = MAX(FLOOR((PartState(2,iPart)-GEO%yminglob)/GEO%FIBGMdeltas(2)),0) + 1
    CellY = MIN(GEO%FIBGMjmax,CellY)
    CellZ = MAX(FLOOR((PartState(3,iPart)-GEO%zminglob)/GEO%FIBGMdeltas(3)),0) + 1
    CellZ = MIN(GEO%FIBGMkmax,CellZ)

    ! check all cells associated with this background mesh cell
    nBGMElems = FIBGM_nElems(CellX,CellY,CellZ)

    ! Multiple potential cells found in BGM. Get closest element barycenter by looping over all elements in BGMcell
    IF (nBGMElems.GT.1) THEN
      Distance     = -HUGE(1.)
      ListDistance = -1
      DO iBGMElem = 1,nBGMElems
        ElemID   = FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+iBGMElem)
        CNElemID = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+iBGMElem))
        ListDistance(iBGMElem) = ElemID

        ! no element associated with BGM elelemt
        IF (ElemID.EQ.-1) &
          CALL ABORT(__STAMP__,'Error during RefMapping: unable to find element associated with BGM element!')

        ! oldElemID was previously checked and particle not found inside
        IF (ElemID.EQ.OldElemID) THEN
          Distance(iBGMElem) = -HUGE(1.)
        ELSE
          Distance(iBGMElem) = SUM((PartState(1:3,iPart)-ElemBaryNGeo(1:3,CNElemID))**2)

          ! Do not consider the element if it is too far away
          IF(Distance(iBGMElem).GT.ElemRadius2NGeo(CNElemID))THEN
            Distance(iBGMElem) = -HUGE(1.)
          END IF
        END IF
      END DO ! nBGMElems
      CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

    ! Only one potential cell found in BGM. No Need to sort
    ELSE IF (nBGMElems.EQ.1) THEN
      Distance(1)     = 0.
      ListDistance(1) = FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+1)
    END IF

    OldXi = PartPosRef(1:3,iPart)
    newXi = HUGE(1.0)
    newElemID = -1

    ! loop through sorted list and start by closest element
    DO iBGMElem=1,nBGMElems
      ! ignore old element and elements out of range
      IF(Distance(iBGMELem).EQ.-HUGE(1.)) CYCLE

      ElemID = ListDistance(iBGMElem)
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
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' xi                     ', PartPosRef(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,1X,E15.8)')    ' epsOneCell             ', epsElement
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' oldxi                  ', oldXi
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' newxi                  ', newXi
          IPWRITE(UNIT_stdOut,'(I0,A)')            ' PartPos:           '
          IPWRITE(UNIt_stdOut,'(I0,A18,1X,A18,1X,A18)') '    min ', ' value ', ' max '
          IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' x', GEO%xminglob, PartState(1,iPart), GEO%xmaxglob
          IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' y', GEO%yminglob, PartState(2,iPart), GEO%ymaxglob
          IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' z', GEO%zminglob, PartState(3,iPart), GEO%zmaxglob
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' LastPartPos            ', LastPartPos(1:3,iPart)
          Vec=PartState(1:3,iPart)-LastPartPos(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,1X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
#ifdef IMPA
          IPWRITE(UNIT_stdOut,'(I0,A,1X,L1)') ' Implicit                ', PartIsImplicit(iPart)
#endif
#if defined(ROS) || defined(IMPA)
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' CurrentStage:    ', iStage
          Vec=PartState(1:3,iPart)-PartStateN(1:3,iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,1X,E15.8)') ' displacementN/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' PartStateN             ', PartStateN(1:3,iPart)
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
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' ParticlePos          ' , partstate(1:3,iPart)
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
              IPWRITE(UNIT_stdOut,'(I0,A,3(1X,I0))')    ' iPart                  ', ipart
              IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' xi                     ', partposref(1:3,ipart)
              IPWRITE(UNIT_stdOut,'(I0,A,1(1X,E15.8))') ' EpsOneCell             ', ElemEpsOneCell(CNTestElem)
              IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' oldxi                  ', oldxi
              IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' newxi                  ', newxi
              IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' LastPartPos            ', LastPartPos(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A)')             ' PartPos:           '
              IPWRITE(UNIt_stdOut,'(I0,A18,1X,A18,1X,A18)')                  '    min ', ' value ', ' max '
              IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' x', GEO%xminglob, PartState(1,iPart), GEO%xmaxglob
              IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' y', GEO%yminglob, PartState(2,iPart), GEO%ymaxglob
              IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' z', GEO%zminglob, PartState(3,iPart), GEO%zmaxglob
              Vec=PartState(1:3,iPart)-LastPartPos(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A,1X,E15.8)') ' displacement /halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
#if defined(ROS) || defined(IMPA)
              IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' CurrentStage:    ', iStage
              Vec=PartState(1:3,iPart)-PartStateN(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A,1X,E15.8)') ' displacementN/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
              IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' PartStateN             ', PartStateN(1:3,iPart)
#if USE_MPI
              inelem=PEM%GlobalElemID(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' elemid-N               ', inelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' elemid-N         ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
!                                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') ' testelem-N            ', testelem
              ELSE
!                IPWRITE(UNIT_stdout,'(I0,A)') ' halo-elem-N = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' testelem-N       ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
                                                               !+ PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#endif
#endif /*ROS or IMPA*/
              IPWRITE(UNIT_StdOut,*) "PEM%GlobalElemID(ipart) =", PEM%GlobalElemID(ipart)
#if defined(IMPA)
              IPWRITE(UNIT_stdOut,'(I0,A,1X,L1)') ' Implicit               ', PartIsImplicit(iPart)
#endif
#if USE_MPI
              inelem=PEM%GlobalElemID(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)') ' PEM%GlobalElemID(ipart) <= PP_nElems: halo-elem = F'
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)') ' PEM%GlobalElemID(ipart)  > PP_nElems: halo-elem = T'
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,'(I0,A)')    '                testelem <= PP_nElems: halo-elem = F'
                IPWRITE(UNIT_stdout,'(I0,A,I0)') '                testelem             ', testelem
              ELSE
                IPWRITE(UNIT_stdout,'(I0,A)')    '                testelem  > PP_nElems: halo-elem = T'
!                IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' testelem         ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
!                                                               + PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#endif
              IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' PartSpecies  ', PartSpecies(iPart)
              CALL abort(__STAMP__ ,'Particle not inside of Element, ipart',ipart)
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
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID,GetCNSideID
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectInterSection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Vars          ,ONLY: SideBCMetrics,ElemToBCSides
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO,ElemRadiusNGeo
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Tracking_Vars      ,ONLY: CartesianPeriodic, TrackInfo
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Utils                       ,ONLY: InsertionSort
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars      ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
USE MOD_Particle_Localization       ,ONLY: PARTHASMOVED
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
INTEGER                       :: SideID,CNSideID,flip
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
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))

! Check if the particle moved at all. If not, tracking is done
CNElemID = GetCNElemID(ElemID)
IF (.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(CNElemID))) THEN
  PEM%GlobalElemID(PartID) = ElemID
  PartisDone               = .TRUE.
  RETURN
END IF

PartTrajectory = PartTrajectory/lengthPartTrajectory

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
    IF (PeriMoved) THEN
      IF(GEO%nPeriodicVectors.EQ.3) CYCLE
      PartTrajectory       = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
      lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
    ELSE
      IF (GEO%nPeriodicVectors.EQ.3) RETURN
    END IF
  ELSE
    PeriMoved = .FALSE.
  END IF
  locAlpha = -1.0
  nInter   = 0

  ! track particle vector until the final particle position is achieved
  ! check if particle can intersect with current side
  DO iLocSide = firstSide,LastSide

    ! SideBCMetrics is sorted by distance. stop if the first side is out of range
    IF (SideBCMetrics(BCSIDE_DISTANCE,ilocSide).GT.lengthPartTrajectory0) EXIT

    ! side potentially in range (halo_eps)
    SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
    CNSideID = GetCNSideID(SideID)
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
      SELECT CASE(SideType(CNSideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectInterSection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                            ,  xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
        CASE(BILINEAR,PLANAR_NONRECT)
          CALL ComputeBiLinearIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                          ,    xi (ilocSide),eta(ilocSide),PartID,    SideID                &
                                          ,    alpha2=alphaOld)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                        ,      xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                        ,      xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
      END SELECT

#ifdef CODE_ANALYZE
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
        IF (PartID.EQ.PARTOUT) THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)')           '     | Output after compute intersection (DoubleCheck refmapping): '
          WRITE(UNIT_stdout,'(2(A,I0),A,L1)') '     | SideType: ',SideType(CNSideID) ,' | SideID: ',SideID,' | Hit: ',isHit
          WRITE(UNIT_stdout,'(2(A,G0))')     '     | Alpha: '   , locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
          WRITE(UNIT_stdout,'((A,G0))')      '     | AlphaOld: ',alphaOld
          WRITE(UNIT_stdout,'(A,2(1X,G0))')   '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
        END IF
      END IF
#endif /*CODE_ANALYZE*/

    ! not double check
    ELSE
      SELECT CASE(SideType(CNSideID))
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
                                        ,      xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
    END SELECT

#ifdef CODE_ANALYZE
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
        IF (PartID.EQ.PARTOUT) THEN
          WRITE(UNIT_stdout,'(30("-"))')
          WRITE(UNIT_stdout,'(A)')           '     | Output after compute intersection (refmapping, BCTracing): '
          WRITE(UNIT_stdout,'(2(A,I0),A,L1)') '     | SideType: ',SideType(CNSideID) ,' | SideID: ',SideID,' | Hit: ',isHit
          WRITE(UNIT_stdout,'(2(A,G0))')     '     | Alpha: '   , locAlpha(ilocSide),' | LengthPartTrajectory: ', lengthPartTrajectory
          WRITE(UNIT_stdout,'(A,2(1X,G0))')   '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
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

        TrackInfo%PartTrajectory(1:3)  = PartTrajectory
        TrackInfo%lengthPartTrajectory = lengthPartTrajectory
        TrackInfo%alpha = locAlpha(ilocSide)
        TrackInfo%xi    = xi(hitlocSide)
        TrackInfo%eta   = eta(hitlocSide)
        CALL GetBoundaryInteraction(PartId,SideID,flip,ElemID,reflected)
        PartTrajectory = TrackInfo%PartTrajectory(1:3)
        lengthPartTrajectory = TrackInfo%lengthPartTrajectory
        locAlpha(ilocSide) = TrackInfo%alpha
        xi(hitlocSide) = TrackInfo%xi
        eta(hitlocSide) = TrackInfo%eta

        ! particle moved to a new element in boundary interaction
        IF(ElemID.NE.OldElemID)THEN
          ! Try to recursively calculate the intersection 1000 times. Threshold might be changed...
          IF (iCount.GE.1000 .AND. MOD(iCount,1000).EQ.0) THEN
            IPWRITE(*,'(I4,A,I0,A,3(1X,I0))') ' WARNING: proc has called BCTracking ',iCount &
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
USE MOD_Eval_xyz                    ,ONLY: TensorProductInterpolation
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID,GetCNSideID
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Localization       ,ONLY: LocateParticleInElement
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectInterSection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars          ,ONLY: SideBCMetrics
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemBaryNGeo
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Vars               ,ONLY: PDM,PartState,LastPartPos
USE MOD_Utils                       ,ONLY: InsertionSort
!USE MOD_Particle_Vars,               ONLY:PartPosRef
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID,CNSideID,locSideList(firstSide:lastSide),hitlocSide
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

PartTrajectory       = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
IF (lengthPartTrajectory.GT.0) PartTrajectory = PartTrajectory/lengthPartTrajectory

locAlpha  = -1.0
nInter    = 0
dolocSide =.TRUE.

! Loop over all sides
DO iLocSide=firstSide,LastSide
  ! track particle vector until the final particle position is achieved
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
  CNSideID = GetCNSideID(SideID)
  locSideList(ilocSide) = ilocSide
  flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                        ,   xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
    CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                           ,xi(ilocSide),eta(ilocSide),PartID     ,SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                      ,     xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                    ,       xi(ilocSide),eta(ilocSide),PartID,flip,SideID)
  END SELECT

  IF (locAlpha(ilocSide).GT.-1.0) THEN
    nInter=nInter+1
  END IF
END DO ! ilocSide

! no intersection found. Try to locate particle manually
IF (nInter.EQ.0) THEN
  PartState(1:3,PartID)   = tmpPos
  LastPartPos(1:3,PartID) = tmpLastPartPos(1:3)
!  IF(PartPosRef(1,PartID).GT. 1.) PartPosRef(1,PartID)= 0.99
!  IF(PartPosRef(1,PartID).LT.-1.) PartPosRef(1,PartID)=-0.99
!  IF(PartPosRef(2,PartID).GT. 1.) PartPosRef(2,PartID)= 0.99
!  IF(PartPosRef(2,PartID).LT.-1.) PartPosRef(2,PartID)=-0.99
!  IF(PartPosRef(3,PartID).GT. 1.) PartPosRef(3,PartID)= 0.99
!  IF(PartPosRef(3,PartID).LT.-1.) PartPosRef(3,PartID)=-0.99
  CALL LocateParticleInElement(PartID,doHalo=.TRUE.)

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


END MODULE MOD_Particle_RefTracking
