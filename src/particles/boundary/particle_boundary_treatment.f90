#include "boltzplatz.h"

MODULE                              MOD_part_boundary                                              !
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
  PRIVATE                                                                                          !
!--------------------------------------------------------------------------------------------------!
!!--------------------------------------------------------------------------------------------------!
!  PUBLIC                        :: ParticleBoundary                                                   !
!#ifdef MPI
!  PUBLIC                        :: Communicate_PIC                                                 !
!#endif
!!--------------------------------------------------------------------------------------------------!
                                                                                                   !
CONTAINS                                                                                           !
                                                                                                   !
!SUBROUTINE ParticleBoundary()                                                                          !
!#ifdef MPI
  !USE MOD_Globals,        ONLY : myRank
  !USE MOD_PICDepo_Vars,   ONLY : DepositionType
  !use mpi
!#endif
  !USE MOD_Particle_Vars !,  ONLY : PDM
  !USE MOD_part_MPI_Vars,  ONLY : PMPIVAR,PMPIExchange,ExtPartsAllocated,NbrOfextParticles
  !USE MOD_part_tools,     ONLY : UpdateNextFreePosition
!!--------------------------------------------------------------------------------------------------!
  !IMPLICIT NONE                                                                                    !
!!--------------------------------------------------------------------------------------------------!
!!#ifdef MPI
!!  INCLUDE 'mpif.h'                                                                                 !
!!#endif
!!--------------------------------------------------------------------------------------------------!
!! argument list declaration                                                                        !
!! Local variable declaration                                                                       !
  !INTEGER                          :: i, j                                                         !
!!--------------------------------------------------------------------------------------------------!
!#ifdef MPI
  !INTEGER                          :: k, ks                                                        !
  !INTEGER                          :: IERROR, allocstat                                            !
  !INTEGER                          :: buffsize                                                     !
  !INTEGER                          :: MPIProcPIC                                                   !
  !INTEGER                          :: MPIGROUPMAP                                                  !
  !INTEGER                          :: CellX,CellY,CellZ                                            !
  !INTEGER, ALLOCATABLE             :: indices(:), shape_indices(:)                                 !
  !LOGICAL                          :: isMPIPart                                                    !
!#endif
!!--------------------------------------------------------------------------------------------------!

!#ifdef MPI
   !IF (DepositionType .EQ. 'shape_function') THEN
       !NbrOfextParticles=0
       !extPartsAllocated=.FALSE.
   !END IF
!!   ALLOCATE(PMPIExchange%MPINbrOfParticles(0:PMPIVAR%nProcs-1), &
!!            PMPIExchange%MPIProcNbr(1:PDM%ParticleVecLength),     &
!!            PMPIExchange%MPITags(1:PDM%ParticleVecLength),        &
!!            indices(1:PDM%ParticleVecLength),      stat=allocstat)
   !ALLOCATE(indices(1:PDM%ParticleVecLength),                    &
            !shape_indices(1:PDM%ParticleVecLength), stat=allocstat)
   !IF (allocStat .NE. 0) THEN
      !WRITE(*,*)'ERROR in Boundary_PIC: cannot allocate MPI_Exchange variables!'; STOP
   !END IF
!! -------------------------------------------------------------------------------------------!
!! Diese Struktur wird so nicht mehr gebraucht. Wir könnten sie durch ein logical Array
!! ersetzen, dass uns sagt, ob es sich um ein MPI-Teilchen handelt oder nicht. 
!! Die Info an welchen Prozess das Teilchen geschickt wird, ist erst mal nicht
!! wichtig, da wir in der MPIParticleSend ja noch über das FIBGM rausfinden wohin
!! das Teilchen kommuniziert werden muss.
!!   PMPIExchange%MPINbrOfParticles = 0
!!   PMPIExchange%MPIProcNbr = -1  
!!   PMPIExchange%MPITags = 0
!! ggf. reicht hier auch nur ein Flag aus, der sagt, ob das Teilchen MPI ist oder nicht. Dann 
!! braucht man da kein Array mehr.
!! -------------------------------------------------------------------------------------------!
!!  indices wird weiterhin gebraucht um schnell an die Teilchen zu kommen, die MPI-Prozessiert 
!!  werden müssen
   !indices(:) = 0
   !shape_indices(:) = 0
!!  Wenn wir die Send-Message in MPI-Parts und Shape-Parts aufteilen wollen, sollten wir
!!  vielleicht ein indices_Shape einführen. Dann kann man die Liste nacheinander befüllen.
!! -------------------------------------------------------------------------------------------!
   !k = 1   ! pure MPI particles
   !ks = 1  ! shapefunction MPI particles
!#endif
!!   WRITE(*,*) PDM%ParticleVecLength
   !DO i = 1,PDM%ParticleVecLength
      !IF (PDM%ParticleInside(i)) THEN
!#ifndef MPI
         !CALL Boundary_PIC_Particle(i)
!#else
!!         CALL Boundary_PIC_Particle(i,PMPIExchange%MPIProcNbr(k),PMPIExchange%MPITags(k))
         !CALL Boundary_PIC_Particle(i,isMPIPart)
         !IF (isMPIPart) THEN     ! physically transferred particle
           !indices(k)=i
           !k=k+1
         !ELSE
           !IF(DepositionType .EQ. 'shape_function') THEN    ! only shape particle
             !IF (PDM%ParticleInside(i).EQV..FALSE.) CYCLE            ! Don't try sending particles that have disappeared (open BC)
             !IF (Species(PartSpecies(i))%ChargeIC.EQ.0) CYCLE        ! Don't deposite neutral particles!
             !CellX = INT((PartState(i,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
             !CellX = MIN(GEO%FIBGMimax,CellX)
             !CellX = MAX(GEO%FIBGMimin,CellX)
             !CellY = INT((PartState(i,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
             !CellY = MIN(GEO%FIBGMkmax,CellY)
             !CellY = MAX(GEO%FIBGMkmin,CellY)
             !CellZ = INT((PartState(i,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
             !CellZ = MIN(GEO%FIBGMlmax,CellZ)
             !CellZ = MAX(GEO%FIBGMlmin,CellZ)
             !IF(ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
               !IF(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1) .GT. 0) THEN
                 !shape_indices(ks) = i
                 !ks = ks + 1
               !END IF
             !END IF
           !END IF
         !END IF
!#endif /* MPI */
      !END IF
   !END DO
!#ifdef MPI
   !k  = k  - 1
   !ks = ks - 1
!!WRITE(*,*) PMPIVAR%iProc, k, ks,PDM%ParticleVecLength
   !IF ((k.GT.0).AND.(ks.GT.0)) THEN
     !! pure MPI- and shapefunction particles
     !CALL MPIParticleSend(k,ks,indices = indices(1:k),shape_indices = shape_indices(1:ks))
   !ELSE IF ((k.GT.0).AND.(.NOT.ks.GT.0)) THEN
     !! pure MPI-particles only 
     !CALL MPIParticleSend(k,ks,indices = indices(1:k))
   !ELSE IF ((.NOT.k.GT.0).AND.(ks.GT.0)) THEN
     !! shapefunction particles only
     !CALL MPIParticleSend(k,ks,shape_indices = shape_indices(1:ks))
   !ELSE
     !! neither pure MPI- nor shapefunction particles
     !CALL MPIParticleSend(k,ks)
   !END IF
   !DEALLOCATE(indices, shape_indices)
!#endif
 !RETURN
!END SUBROUTINE ParticleBoundary


!#ifdef MPI
!SUBROUTINE Communicate_PIC()                                                                       !
   !USE MOD_PICDepo_Vars          ,ONLY : DepositionType
   !USE MOD_part_tools            ,ONLY : UpdateNextFreePosition                                    !
   !USE MOD_part_MPI_Vars         ,ONLY : PMPIVAR,PMPIExchange
   !use mpi
!!--------------------------------------------------------------------------------------------------!
   !IMPLICIT NONE                                                                                   !
!!--------------------------------------------------------------------------------------------------!
   !!INCLUDE 'mpif.h'                                                                                !
!!--------------------------------------------------------------------------------------------------!
!! argument list declaration                                                                        !
!! Local variable declaration                                                                       !
!!--------------------------------------------------------------------------------------------------!
   !INTEGER                          :: iProc, nSend, nRecv                                         !
   !INTEGER                          :: send_status(1:MPI_STATUS_SIZE)                              !
   !INTEGER                          :: IERROR                                                      !
   !INTEGER                          :: totalNSendParticles,totalNRecvParticles                     !
   !REAL                             :: localDQ                                                     !
   !REAL                             :: myAnodeAbsorbtionPoint                                      !
!!--------------------------------------------------------------------------------------------------!
!!   nSend=PMPIExchange%nMPIParticles
   !CALL MPIParticleRecv(nRecv)
   !DO iProc = 0,PMPIVAR%nProcs-1
     !!--- Don't send messages to myself!
     !IF ((iProc.EQ.PMPIVAR%iProc).AND.(DepositionType.NE.'shape_function')) CYCLE
     !IF (.NOT.(PMPIVAR%MPINeighbor(iProc))) CYCLE
     !!--- Wait for number of communicated particles message to be complete
     !CALL MPI_WAIT(PMPIExchange%send_request(iProc,1),send_status(:),IERROR)
     !IF (SUM(PMPIExchange%nbrOfSendParticles(iProc,:)) .GT. 0) THEN
       !CALL MPI_WAIT(PMPIExchange%send_request(iProc,2),send_status(:),IERROR)
       !!--- be nice: deallocate the send buffer
       !DEALLOCATE( PMPIExchange%send_message(iProc)%content )
     !END IF
   !END DO
   !!--- deallocate message buffers and request buffers
   !DEALLOCATE(PMPIExchange%send_message,PMPIExchange%send_request,PMPIExchange%nbrOfSendParticles,PMPIExchange%NbrArray)
!!   !--- error checking


!!   CALL MPI_ALLREDUCE(nSend,totalNSendParticles,1,MPI_INTEGER,MPI_SUM,PMPIVAR%COMM,IERROR)
!!   IF ((PMPIVAR%iProc.EQ.0)) THEN
!!     WRITE(*,*) 'number of sent particles: ',totalNSendParticles
!!   END IF
   !!CALL UpdateNextFreePosition() ! now is done in the timedisc 
   !RETURN
!END SUBROUTINE Communicate_PIC
!#endif

!#ifdef MPI
!!SUBROUTINE Boundary_PIC_Particle(i,MPIBoundary,MPITag)                                             !
!SUBROUTINE Boundary_PIC_Particle(i,isMPIPart)                                                      ! 
!#else
!SUBROUTINE Boundary_PIC_Particle(i)                                                                !
!#endif
!!DEC$ ATTRIBUTES FORCEINLINE :: Boundary_PIC_Particle
!USE MOD_Particle_Vars!, ONLY: 
!USE MOD_Mesh_Vars,     ONLY : nInnerSides, nBCSides, BC, ElemToSide, SideToElem
!#ifdef MPI
!USE MOD_part_MPI_Vars, ONLY : MPIGEO,PMPIVAR
!#endif
  !USE MOD_BoundaryTools, ONLY : ParticleThroughSideCheck3DFast, ParticleThroughSideCheck3DFastold
!USE MOD_BoundaryTools, ONLY : ParticleInsideQuad3D
!USE MOD_BoundaryTools, ONLY : ParticleInsideQuad3Dold
!USE MOD_BoundaryTools, ONLY : SingleParticleToExactElement
!USE MOD_BoundaryTools, ONLY : PerfectReflection3D
!USE MOD_BoundaryTools, ONLY : DiffuseReflection3D 
!USE MOD_BoundaryTools, ONLY : PeriodicWallBnd3D 
!USE MOD_BoundaryTools, ONLY : ParticleThroughSideLastPosCheck
!USE MOD_part_MPI_Vars, ONLY : partShiftVector
!USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut,PartAnalyzeStep
!USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
!USE MOD_Particle_Surfaces_Vars, ONLY:nTriangles
!!--------------------------------------------------------------------------------------------------!
!!--------------------------------------------------------------------------------------------------!
  !IMPLICIT NONE                                                                                    !
!!--------------------------------------------------------------------------------------------------!
!! argument list declaration                                                                        !
  !INTEGER, INTENT(IN)              :: i                                                            !
!#ifdef MPI
!!  INTEGER, INTENT(OUT)             :: MPIBoundary, MPITag                                          !
  !LOGICAL, INTENT(OUT)             :: isMPIPart                                                    !
!#endif
!!--------------------------------------------------------------------------------------------------!
!! Local variable declaration                                                                       !
  !INTEGER                          :: n, k, Element, LocalSide, tempcount
  !INTEGER                          :: NrOfThroughSides, ind, ind2
  !INTEGER                          :: GlobSideID,TempSideID,iLocSide                    !
  !INTEGER                          :: TriNum, LocSidesTemp(nTriangles*6),TriNumTemp(nTriangles*6)                    !
  !INTEGER                          :: PV, PVsign,SecondNrOfThroughSides
  !INTEGER                          :: DoneSideID(1:2)  ! 1 = Side, 2 = TriNum
  !INTEGER                          :: DoneLastElem(1:3,1:2) ! 1:3: 1=Element,2=LocalSide,3=TriNum 1:2: 1=last 2=beforelast
  !INTEGER                          :: DoneLastElemHC(1:3) !see above, for haloCells 
  !LOGICAL                          :: ThroughSide, InElementCheck, DONE                            !
  !REAL                             :: det(6,nTriangles),detM,ratio,minRatio
  !REAL, PARAMETER                  :: eps = 0                                                      !
  !REAL                             :: RanNum
!!--------------------------------------------------------------------------------------------------!
!IF((GEO%nPeriodicVectors.GE.0).AND.(ALLOCATED(partShiftVector))) partShiftVector(1:3,i) = 0.
!#ifdef MPI
!!   MPIBoundary = -1
!!   MPITag = 0
   !isMPIPart = .FALSE.
!#endif
   !DONE = .FALSE.
   !Element = PEM%lastElement(i)
   !!print*,Element
   !GlobSideID = 0
   !DoneSideID(:) = 0
   !DoneLastElem(:,:) = 0
   !DoneLastElemHC(:) = 0
!!tempcount = 0

   !DO WHILE (.NOT.DONE)
!!!---this can/will be removed as soon as there are no more localization problems
!!tempcount = tempcount + 1
!!IF(tempcount.GT.100)THEN
!!WRITE(*,*) 'Part',i
!!STOP
!!END IF
!!---------------------------------------------------------------------------
      !!---- Check whether particle is in element
      !!CALL ParticleInsideQuad3Dold(i,Element,InElementCheck,det)
      !CALL ParticleInsideQuad3D(i,Element,InElementCheck,det)
  !!    print*,InElementCheck
      !!---- If it is, set new ElementNumber = Element and LocalizeOn = .FALSE. ->DONE
      !IF (InElementCheck) THEN
         !PEM%Element(i) = Element
         !DONE = .TRUE.
!!               WRITE(*,*)'Particle ',i,' -> element ',Element

      !!---- If it is not, check through which side it moved

      !ELSE
         !NrOfThroughSides = 0
         !LocSidesTemp(:) = 0
         !TriNumTemp(:) = 0
         !DO iLocSide=1,6
            !TempSideID=ElemToSide(E2S_SIDE_ID,iLocSide,Element)
            !DO TriNum = 1,nTriangles
              !IF (det(iLocSide,TriNum).le.-eps) THEN
                !IF((TempSideID.EQ.DoneSideID(1)).AND.(TriNum.EQ.DoneSideID(2))) CYCLE  !necessary??? test one day
                !ThroughSide = .FALSE.
                !!print*,'-----------------------------------------------------------------'
                !!print*,'old'
                !!CALL ParticleThroughSideCheck3DFastold(i,iLocSide,Element,ThroughSide,TriNum)
                !!print*,'-----------------------------------------------------------------'
                !!print*,'new'
                !CALL ParticleThroughSideCheck3DFast(i,iLocSide,Element,ThroughSide,TriNum)
                !IF (ThroughSide) THEN
                  !NrOfThroughSides = NrOfThroughSides + 1
                  !LocSidesTemp(NrOfThroughSides) = iLocSide
                  !TriNumTemp(NrOfThroughSides) = TriNum
                  !GlobSideID = TempSideID
                  !LocalSide = iLocSide
                !END IF
              !END IF
            !END DO
         !END DO
         !TriNum = TriNumTemp(1)
!!!$IF(i.EQ.676378)THEN
!!!$WRITE(*,*) 'el',Element
!!!$WRITE(*,*) PartState(i,1:3)
!!!$WRITE(*,*) lastPartPos(i,1:3)
!!!$WRITE(*,*) LocSidesTemp
!!!$WRITE(*,*) TriNumTemp
!!!$END IF
         !!--- if no side is found use the slower search method
         !!--- if more than one is found, figure out which one it is
         !IF (NrOfThroughSides.NE.1) THEN
           !IF (NrOfThroughSides.EQ.0) THEN    !no side
             !GlobSideID = 0
             !WRITE(*,*) 'Error in Iteration-Step ??? ! Particle Number',i,'lost. Searching for particle....'
             !WRITE(*,*) 'Element: ', Element
             !WRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
             !WRITE(*,*) 'Pos:     ', PartState(i,1:3)
             !WRITE(*,*) 'Velo:    ', PartState(i,4:6)
             !CALL SingleParticleToExactElement(i)
             !! Retrace to check through which side the particle went
             !DO iLocSide=1,6
               !TempSideID=ElemToSide(E2S_SIDE_ID,iLocSide,Element)
               !IF(ElemToSide(E2S_FLIP,iLocSide,Element).EQ.0) THEN
                 !IF(SideToElem(S2E_NB_ELEM_ID,TempSideID).EQ.PEM%Element(i)) GlobSideID = TempSideID
               !ELSE
                 !IF(SideToElem(S2E_ELEM_ID   ,TempSideID).EQ.PEM%Element(i)) GlobSideID = TempSideID
               !END IF
             !END DO
             !IF(.NOT.PDM%ParticleInside(i))THEN
               !WRITE(*,*)'Particle',i,' lost completely!'
               !WRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
               !WRITE(*,*) 'Pos:     ', PartState(i,1:3)
               !WRITE(*,*) 'Velo:    ', PartState(i,3:6)
               !PDM%ParticleInside(i) = .FALSE.
               !GlobSideID = 0
             !ELSE
              !WRITE(*,*) '...Particle found again'
              !WRITE(*,*) 'Element: ', PEM%Element(i)
             !END IF
             !DONE = .TRUE.
            !! dirty bug fix. This kills application if particle is found directly at edge with
            !! velocity vector tang. to edge
!!             IF (WeirdElems.EQ.0) STOP  ! temporary stop to identify problems:
           !ELSE IF (NrOfThroughSides.GT.1) THEN   ! more  than one side (possible for irregular hexagons
             !SecondNrOfThroughSides = 0
             !minRatio = 0
             !DO ind2 = 1, NrOfThroughSides
               !IF(.NOT.((DoneLastElem(1,2).EQ.Element).AND. &
                        !(DoneLastElem(2,2).EQ.LocSidesTemp(ind2)).AND. &
                        !(DoneLastElem(3,2).EQ.TriNumTemp(ind2)))) THEN
                 !CALL ParticleThroughSideLastPosCheck(i,LocSidesTemp(ind2),Element,InElementCheck,TriNumTemp(ind2),detM)
                 !IF (InElementCheck) THEN
                   !IF((detM.EQ.0).AND.(det(LocSidesTemp(ind2),TriNumTemp(ind2)).EQ.0)) CYCLE ! particle moves within side
                   !IF((detM.EQ.0).AND.(minRatio.EQ.0))THEN !safety measure
                     !SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                     !GlobSideID = ElemToSide(E2S_SIDE_ID,LocSidesTemp(ind2),Element)
                     !LocalSide = LocSidesTemp(ind2)
                     !TriNum = TriNumTemp(ind2)
                   !ELSE
                     !!--- compare ratio of spatial product of PartPos->Tri-Nodes and LastPartPos->Tri-Nodes
                     !ratio = det(LocSidesTemp(ind2),TriNumTemp(ind2))/detM
                     !IF (ratio.LT.minRatio) THEN ! ratio is always negative, i.e. maximum abs is wanted!
                       !minRatio = ratio
                       !SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                       !GlobSideID = ElemToSide(E2S_SIDE_ID,LocSidesTemp(ind2),Element)
                       !LocalSide = LocSidesTemp(ind2)
                       !TriNum = TriNumTemp(ind2)
                     !END IF
                   !END IF
                 !END IF
               !END IF
             !END DO
             !IF (SecondNrOfThroughSides.EQ.0) THEN
               !WRITE(*,*) 'Error in Boundary_treatment: Particle',i,'went through no Sides on second check'
               !WRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
               !WRITE(*,*) 'Pos:     ', PartState(i,1:3)
               !WRITE(*,*) 'Velo:    ', PartState(i,4:6)
               !STOP
             !END IF
           !END IF
         !END IF
!!!$IF(i.EQ.676378)THEN
!!!$WRITE(*,*) 'final',Element, LocalSide,TriNum
!!!$WRITE(*,*) 'side', GlobSideID, nBCSides, nInnerSides
!!!$END IF

         !!---- Is this side a boundary side?
         !IF (GlobSideID.LE.0) THEN
           !WRITE(*,*)'Boundary_PIC: Could not identify through which side Particle passed.'
           !WRITE(*,*)'              Skipping boundary condition for this particle.'
           !EXIT
         !END IF
         !IF ((GlobSideID.GT.0).AND.(GlobSideID.LE.nBCSides)) THEN
           !!---- Side is boundary side: Call corresponding boundary subroutine
           !IF (.NOT. ASSOCIATED(PartBound%Map)) THEN
             !WRITE(*,*)'ERROR: PartBound not associated!'
             !STOP
           !END IF
           !IF (PartBound%Map(BC(GlobSideID)).EQ.PartBound%OpenBC) THEN
             !IF(CalcPartBalance) THEN
               !!IF(MOD(iter+1,PartAnalyzeStep).EQ.0)THEN ! caution if correct
                 !nPartOut(PartSpecies(i))=nPartOut(PartSpecies(i)) + 1
                 !PartEkinOut(PartSpecies(i))=PartEkinOut(PartSpecies(i))+CalcEkinPart(i)
               !!END IF ! iter+1
             !END IF ! CalcPartBalance
             !PDM%ParticleInside(i) = .FALSE.
             !DONE = .TRUE.
           !ELSE IF (PartBound%Map(BC(GlobSideID)).EQ.PartBound%ReflectiveBC) THEN
!!IF(i.EQ.788729)THEN
!!WRITE(*,*) 'YESSSSSSSSSSS'
!!END IF
              !print*,'we hit boundary condition'
             !CALL RANDOM_NUMBER(RanNum)
             !IF(RanNum.GE.PartBound%MomentumACC(BC(GlobSideID))) THEN
!!.... Specular re-emission
               !CALL PerfectReflection3D(i,LocalSide,Element,TriNum,PartBound%WallVelo(1:3,BC(GlobSideID)))
               !DoneSideID(1) = GlobSideID
               !DoneSideID(2) = TriNum
               !DoneLastElem(:,:) = 0
               !DoneLastElemHC(:) = 0
             !ELSE
!!.... Diffuse re-emission
               !CALL DiffuseReflection3D (i,LocalSide,Element,TriNum,            &
                                         !PartBound%WallTemp(BC(GlobSideID)),    &
                                         !PartBound%TransACC(BC(GlobSideID)),    &
                                         !PartBound%VibACC(BC(GlobSideID)),      &
                                         !PartBound%RotACC(BC(GlobSideID)),       &
                                         !PartBound%WallVelo(1:3,BC(GlobSideID)) )
               !DoneSideID(1) = GlobSideID
               !DoneSideID(2) = TriNum
               !DoneLastElem(:,:) = 0
               !DoneLastElemHC(:) = 0
             !END IF
           !ELSE IF (PartBound%Map(BC(GlobSideID)).EQ.PartBound%SimpleAnodeBC) THEN
             !WRITE(*,*)'ERROR: SimpleAnodeBC not implemented yet';STOP
           !ELSE IF (PartBound%Map(BC(GlobSideID)).EQ.PartBound%SimpleCathodeBC) THEN
             !WRITE(*,*)'ERROR: SimpleCathodeBC not implemented yet';STOP
           !END IF
           
!#ifdef MPI
         !!---- If it is not a domain boundary side - it is a MPI boundary
         !ELSE IF (GlobSideID.GT.nInnerSides+nBCSides) THEN
           !! ---AS: do tracing of particle in halo cells
           !! ---TS: if this side is additionally a periodic side, we first shift the particle
           !IF(GEO%PeriodicElemSide(LocalSide,Element).NE.0) THEN
             !PV = abs(GEO%PeriodicElemSide(LocalSide,Element))
             !PVsign = INT(GEO%PeriodicElemSide(LocalSide,Element)/PV)
             !lastPartPos(i,1) = lastPartPos(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
             !lastPartPos(i,2) = lastPartPos(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
             !lastPartPos(i,3) = lastPartPos(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
             !PartState(i,1)   = PartState(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
             !PartState(i,2)   = PartState(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
             !PartState(i,3)   = PartState(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
             !partShiftVector(1:3,i) = partShiftVector(1:3,i) + GEO%PeriodicVectors(1:3,PV)*PVsign
           !END IF
           !DoneLastElem(1,:) = Element
           !DoneLastElem(2,:) = LocalSide
           !DoneLastElem(3,:) = TriNum
           !CALL Boundary_Particle_halocells(i,GlobSideID,isMPIPart,TriNum, DoneLastElemHC)
!!!$IF(i.EQ.676378)THEN
!!!$WRITE(*,*) PMPIVAR%iProc,isMPIPart, PDM%ParticleInside(i)
!!!$END IF
           !IF ((.NOT.isMPIPart).AND.(PDM%ParticleInside(i))) THEN ! reentered into myproc domain
             !! map from halo cells back to myproc domain
             !DoneSideID(1) = -MPIGEO%BC(4,GlobSideID)
             !IF(TriNum.EQ.1) DoneSideID(2) = 2      ! WILL ONLY WORK IF POINT NR1 IS SAME ON halocells
             !IF(TriNum.EQ.2) DoneSideID(2) = 1
             !IF (SideToElem(S2E_ELEM_ID,DoneSideID(1)).NE.-1) THEN
               !Element = SideToElem(S2E_ELEM_ID,DoneSideID(1))
             !ELSE
               !Element = SideToElem(S2E_NB_ELEM_ID,DoneSideID(1))
             !END IF
           !ELSE 
             !DONE = .TRUE.
           !END IF
!#endif
         !!---- If it is not a boundary side, it is either an inner side or periodic
         !ELSE   ! Inner side, GlobSideID here is also .GT. nBCSides
           !IF(GEO%PeriodicElemSide(LocalSide,Element).EQ.0) THEN ! inner side
             !DoneSideID(1) = GlobSideID
             !IF(TriNum.EQ.1) DoneSideID(2) = 2
             !IF(TriNum.EQ.2) DoneSideID(2) = 1
             !DoneLastElem(:,2) = DoneLastElem(:,1)
             !DoneLastElem(1,1) = Element
             !DoneLastElem(2,1) = LocalSide
             !DoneLastElem(3,1) = TriNum
             !IF (SideToElem(S2E_NB_ELEM_ID,GlobSideID).EQ.Element) THEN
               !Element = SideToElem(S2E_ELEM_ID,GlobSideID)
             !ELSE
               !Element = SideToElem(S2E_NB_ELEM_ID,GlobSideID)
             !END IF
           !ELSE  ! periodic side
             !PV = abs(GEO%PeriodicElemSide(LocalSide,Element))
             !PVsign = INT(GEO%PeriodicElemSide(LocalSide,Element)/PV)
             !lastPartPos(i,1) = lastPartPos(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
             !lastPartPos(i,2) = lastPartPos(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
             !lastPartPos(i,3) = lastPartPos(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
             !PartState(i,1)   = PartState(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
             !PartState(i,2)   = PartState(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
             !PartState(i,3)   = PartState(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
             !partShiftVector(1:3,i) = partShiftVector(1:3,i) + GEO%PeriodicVectors(1:3,PV)*PVsign
             !DoneSideID(1) = GlobSideID
             !IF(TriNum.EQ.1) DoneSideID(2) = 2
             !IF(TriNum.EQ.2) DoneSideID(2) = 1
             !DoneLastElem(:,2) = DoneLastElem(:,1)
             !DoneLastElem(1,1) = Element
             !DoneLastElem(2,1) = LocalSide
             !DoneLastElem(3,1) = TriNum
             !IF (SideToElem(S2E_NB_ELEM_ID,GlobSideID).EQ.Element) THEN
               !Element = SideToElem(S2E_ELEM_ID,GlobSideID)
             !ELSE
               !Element = SideToElem(S2E_NB_ELEM_ID,GlobSideID)
             !END IF
           !END IF
         !END IF
      !END IF
   !END DO
 !RETURN
!END SUBROUTINE Boundary_PIC_Particle

!#ifdef MPI
!SUBROUTINE Boundary_Particle_halocells(i,GlobSideID,isMPIPart,TriNum,DoneLastElemHC)
!!DEC$ ATTRIBUTES FORCEINLINE :: Boundary_PIC_Particle
!USE MOD_Particle_Vars!, ONLY: 
!USE MOD_Mesh_Vars,     ONLY : nInnerSides, nBCSides, BC, nSides
!USE MOD_part_MPI_Vars, ONLY : MPIGEO, PMPIVAR 
!USE MOD_BoundaryTools, ONLY : ParticleThroughSideCheck3DFast_halocells
!USE MOD_BoundaryTools, ONLY : ParticleInsideQuad3D_halocells
!USE MOD_BoundaryTools, ONLY : PerfectReflection3D_halocells
!USE MOD_BoundaryTools, ONLY : DiffuseReflection3D_halocells 
!USE MOD_BoundaryTools, ONLY : ParticleThroughSideLastPosCheck_halocells
!USE MOD_part_MPI_Vars, ONLY : partShiftVector
!USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartOut,PartEkinOut,PartAnalyzeStep
!USE MOD_Particle_Analyze       ,ONLY: CalcEkinPart
!!--------------------------------------------------------------------------------------------------!
  !IMPLICIT NONE                                                                                    !
!!--------------------------------------------------------------------------------------------------!
!! argument list declaration                                                                        !
  !INTEGER, INTENT(IN)              :: i                                                            !
  !INTEGER, INTENT(INOUT)           :: GlobSideID,TriNum, DoneLastElemHC(1:3)
  !LOGICAL, INTENT(OUT)             :: isMPIPart
!!--------------------------------------------------------------------------------------------------!
!! Local variable declaration                                                                       !
  !INTEGER                          :: n, k, Element, LocalSide, tempcount !
  !INTEGER                          :: NrOfThroughSides,DoneSideID(2), SecondNrOfThroughSides
  !INTEGER                          :: TempSideID,iLocSide, ind2, ind
  !INTEGER                          :: haloSideID, DoneLastElem(1:3,1:2)
  !INTEGER                          :: PV, PVsign, LocSidesTemp(1:6),TriNumTemp(1:6)
  !LOGICAL                          :: ThroughSide, InElementCheck, DONE, INNERFLAG
  !REAL                             :: det(6,2), detM,ratio,minRatio
  !REAL, PARAMETER                  :: eps = 0                                                      !
  !REAL                             :: RanNum
!!--------------------------------------------------------------------------------------------------!
   !isMPIPart = .FALSE.
   !DONE = .FALSE.
   !! do side mapping from inner cells to halo cells:
   !haloSideID = MPIGEO%haloMPINbSide(GlobSideID-nInnerSides-nBCSides)
   !! get halo cells ElemID. Not quite sure if ELEM_ID or NB_ELEM_ID. 
   !! maybe check this later and delete IF-query here to save cpu time
   !IF (MPIGEO%SideToElem(S2E_ELEM_ID,haloSideID).NE.-1) THEN
     !Element = MPIGEO%SideToElem(S2E_ELEM_ID,haloSideID)
   !ELSE
     !Element = MPIGEO%SideToElem(S2E_NB_ELEM_ID,haloSideID)
   !END IF

   !DoneSideID(1) = haloSideID
   !IF(TriNum.EQ.1) DoneSideID(2) = 2
   !IF(TriNum.EQ.2) DoneSideID(2) = 1
   !DoneLastElem(1,:) = DoneLastElemHC(1)
   !DoneLastElem(2,:) = DoneLastElemHC(2)
   !DoneLastElem(3,:) = DoneLastElemHC(3)
!!!$tempcount = 0
   !DO WHILE (.NOT.DONE)
!!!$tempcount = tempcount + 1
!!!$IF((tempcount.GT.5).AND.(i.EQ.600951))THEN
!!!$WRITE(*,*) 'partHC',i
!!!$STOP
!!!$END IF



     !!---- Check whether particle is in element

     !CALL ParticleInsideQuad3D_halocells(i,Element,InElementCheck,det)

     !!---- If it is, set new ElementNumber = Element and LocalizeOn = .FALSE. ->DONE
     !IF (InElementCheck) THEN
       !PEM%Element(i) = Element
       !isMPIPart = .TRUE.
       !DONE = .TRUE.
     !ELSE !---- If it is not, check through which side it moved
       !NrOfThroughSides = 0
       !LocSidesTemp(:) = 0
       !TriNumTemp(:) = 0
       !DO iLocSide=1,6
         !TempSideID=MPIGEO%ElemToSide(E2S_SIDE_ID,iLocSide,Element)
         !DO TriNum = 1,2
           !IF (det(iLocSide,TriNum).le.-eps) THEN
             !IF((TempSideID.EQ.DoneSideID(1)).AND.(TriNum.EQ.DoneSideID(2))) CYCLE
             !ThroughSide = .FALSE.
             !CALL ParticleThroughSideCheck3DFast_halocells(i,iLocSide,Element,ThroughSide,TriNum)
             !IF (ThroughSide) THEN
               !NrOfThroughSides = NrOfThroughSides + 1
               !LocSidesTemp(NrOfThroughSides) = iLocSide
               !TriNumTemp(NrOfThroughSides) = TriNum
               !haloSideID = TempSideID
               !LocalSide = iLocSide
             !END IF
           !END IF
         !END DO
       !END DO
       !TriNum = TriNumTemp(1)
!!!$IF(i.EQ.600951)THEN
!!!$WRITE(*,*) 'elHC',Element
!!!$WRITE(*,*) PartState(i,1:3)
!!!$WRITE(*,*) lastPartPos(i,1:3)
!!!$WRITE(*,*) LocSidesTemp
!!!$WRITE(*,*) TriNumTemp
!!!$END IF

       !!--- if no side is found use the slower search method
       !IF ((NrOfThroughSides.eq.0).OR.(haloSideID.LE.0)) THEN
         !WRITE(*,*)'Particle',i,' lost in halocell!'
         !WRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
         !WRITE(*,*) 'Pos:     ', PartState(i,1:3)
         !WRITE(*,*) 'Velo:    ', PartState(i,3:6)
         !PDM%ParticleInside(i) = .false.
         !DONE=.TRUE.
         !EXIT
       !END IF
       !IF (NrOfThroughSides.GT.1) THEN
         !SecondNrOfThroughSides = 0
         !minRatio = 0
         !DO ind2 = 1, NrOfThroughSides
           !IF(.NOT.((DoneLastElem(1,2).EQ.Element).AND. &
                    !(DoneLastElem(2,2).EQ.LocSidesTemp(ind2)).AND. &
                    !(DoneLastElem(3,2).EQ.TriNumTemp(ind2)))) THEN
             !CALL ParticleThroughSideLastPosCheck_halocells(i,LocSidesTemp(ind2),Element,InElementCheck,TriNumTemp(Ind2),detM)
!!!$IF(i.EQ.600951)THEN
!!!$WRITE(*,*) 'HCdet', LocSidesTemp(ind2),TriNumTemp(ind2)
!!!$WRITE(*,*) detM
!!!$END IF
             !IF (InElementCheck) THEN
               !ratio = det(LocSidesTemp(ind2),TriNumTemp(ind2))/detM
               !IF (ratio.LT.minRatio) THEN ! ratio is always negative, i.e. maximum abs is wanted!
                 !minRatio = ratio
                 !SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                 !haloSideID = MPIGEO%ElemToSide(E2S_SIDE_ID,LocSidesTemp(ind2),Element)
                 !LocalSide = LocSidesTemp(ind2)
                 !TriNum = TriNumTemp(ind2)
               !END IF
             !END IF
           !END IF
         !END DO
         !IF (SecondNrOfThroughSides.EQ.0) THEN
           !WRITE(*,*) 'Error in Boundary_treatment_halocells: Particle',i,'went through no Sides on second check'
           !WRITE(*,*) 'LastPos: ', LastPartPos(i,1:3)
           !WRITE(*,*) 'Pos:     ', PartState(i,1:3)
           !WRITE(*,*) 'Velo:    ', PartState(i,3:6)
           !PDM%ParticleInside(i) = .false.
           !DONE=.TRUE.
           !EXIT
         !END IF
       !END IF
       !IF (MPIGEO%BC(1,haloSideID).NE.0) THEN
         !IF ((MPIGEO%BC(1,haloSideID).NE.-1).AND.(MPIGEO%BC(1,haloSideID).NE.424242)) THEN 
           !!---- Side is boundary side: Call corresponding boundary subroutine
           !IF (PartBound%Map(MPIGEO%BC(1,haloSideID)).EQ.PartBound%OpenBC) THEN
             !PDM%ParticleInside(i) = .FALSE.
             !DONE = .TRUE.
             !IF(CalcPartBalance) THEN
               !!IF(MOD(iter+1,PartAnalyzeStep).EQ.0)THEN ! caution if correct
                 !nPartOut(PartSpecies(i))=nPartOut(PartSpecies(i)) + 1
                 !PartEkinOut(PartSpecies(i))=PartEkinOut(PartSpecies(i))+CalcEkinPart(i)
               !!END IF ! iter+1
             !END IF ! CalcPartBalance
           !ELSE IF (PartBound%Map(MPIGEO%BC(1,haloSideID)).EQ.PartBound%ReflectiveBC) THEN
             !CALL RANDOM_NUMBER(RanNum)
             !IF(RanNum.GE.PartBound%MomentumACC(MPIGEO%BC(1,haloSideID))) THEN
!!.... Specular re-emission
               !CALL PerfectReflection3D_halocells(i,LocalSide,Element,TriNum,PartBound%WallVelo(1:3,MPIGEO%BC(1,haloSideID)))
               !DoneSideID(1) = haloSideID
               !DoneSideID(2) = TriNum
               !DoneLastElem(:,:) = 0
               !DoneLastElemHC(:) = 0
             !ELSE
!!.... Diffuse re-emission
               !CALL DiffuseReflection3D_halocells(i,LocalSide,Element,TriNum,         &
                                                !PartBound%WallTemp(MPIGEO%BC(1,haloSideID)), &
                                                !PartBound%TransACC(MPIGEO%BC(1,haloSideID)), &
                                                !PartBound%VibACC(MPIGEO%BC(1,haloSideID)),   &
                                                !PartBound%RotACC(MPIGEO%BC(1,haloSideID)),   &
                                                !PartBound%WallVelo(1:3,MPIGEO%BC(1,haloSideID))       )
               !DoneSideID(1) = haloSideID
               !DoneSideID(2) = TriNum
               !DoneLastElem(:,:) = 0
               !DoneLastElemHC(:) = 0
             !END IF
           !ELSE IF (PartBound%Map(MPIGEO%BC(1,haloSideID)).EQ.PartBound%SimpleAnodeBC) THEN
             !WRITE(*,*)'ERROR: SimpleAnodeBC not implemented yet';STOP
           !ELSE IF (PartBound%Map(MPIGEO%BC(1,haloSideID)).EQ.PartBound%SimpleCathodeBC) THEN
             !WRITE(*,*)'ERROR: SimpleCathodeBC not implemented yet';STOP
           !END IF
         !ELSE IF (MPIGEO%BC(1,haloSideID).EQ.-1) THEN ! MPI-Side (case 1 and 3)
           !! ---TS: if this side is additionally a periodic side, we first shift the particle
           !IF (GEO%nPeriodicVectors.NE.0) THEN
             !IF(MPIGEO%PeriodicElemSide(LocalSide,Element).NE.0) THEN
               !PV = abs(MPIGEO%PeriodicElemSide(LocalSide,Element))
               !PVsign = INT(MPIGEO%PeriodicElemSide(LocalSide,Element)/PV)
               !lastPartPos(i,1) = lastPartPos(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
               !lastPartPos(i,2) = lastPartPos(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
               !lastPartPos(i,3) = lastPartPos(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
               !PartState(i,1)   = PartState(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
               !PartState(i,2)   = PartState(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
               !PartState(i,3)   = PartState(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
               !partShiftVector(1:3,i) = partShiftVector(1:3,i) + GEO%PeriodicVectors(1:3,PV)*PVsign
             !END IF
           !END IF
           !IF (MPIGEO%BC(2,haloSideID).EQ.PMPIVAR%iProc) THEN ! case1: reentered into my domain
             !DONE = .TRUE.
             !GlobSideID = haloSideID
             !DoneLastElemHC(1) = Element
             !DoneLastElemHC(2) = LocalSide
             !DoneLastElemHC(3) = TriNum
             !RETURN
           !ELSE ! case 3: crossed MPI-bound in halo cells
             !DoneSideID(1) = -MPIGEO%BC(4,haloSideID)
             !IF(TriNum.EQ.1) DoneSideID(2) = 2
             !IF(TriNum.EQ.2) DoneSideID(2) = 1
             !DoneLastElem(:,2) = DoneLastElem(:,1)
             !DoneLastElem(1,1) = Element
             !DoneLastElem(2,1) = LocalSide
             !DoneLastElem(3,1) = TriNum
             !IF (MPIGEO%SideToElem(S2E_ELEM_ID,DoneSideID(1)).NE.-1) THEN
               !Element = MPIGEO%SideToElem(S2E_ELEM_ID,DoneSideID(1))
             !ELSE
               !Element = MPIGEO%SideToElem(S2E_NB_ELEM_ID,DoneSideID(1))
             !END IF
           !END IF
         !ELSE IF (MPIGEO%BC(1,haloSideID).EQ.424242) THEN  !case 2: left halo cells in wrong direction
           !WRITE(*,*) 'Particle left halo cell domain:'
           !WRITE(*,*) 'Part',i,'Proc:',PMPIVAR%iProc
           !WRITE(*,*) 'PartState',PartState(i,:)
           !WRITE(*,*) "Error: particle left halo cell domain. Increase halo radius of proc ",PMPIVAR%iProc
           !STOP
         !END IF
       !ELSE   ! Inner or periodic side
         !INNERFLAG = .FALSE.
         !IF(GEO%nPeriodicVectors.NE.0)THEN  !periodic is possible
           !IF(MPIGEO%PeriodicElemSide(LocalSide,Element).EQ.0) INNERFLAG = .TRUE.
         !ELSE
           !INNERFLAG = .TRUE.  ! only inner is possible
         !END IF
         !IF(INNERFLAG) THEN ! inner side
           !DoneSideID(1) = haloSideID
           !IF(TriNum.EQ.1) DoneSideID(2) = 2
           !IF(TriNum.EQ.2) DoneSideID(2) = 1
           !DoneLastElem(:,2) = DoneLastElem(:,1)
           !DoneLastElem(1,1) = Element
           !DoneLastElem(2,1) = LocalSide
           !DoneLastElem(3,1) = TriNum
           !IF (MPIGEO%SideToElem(S2E_NB_ELEM_ID,haloSideID).EQ.Element) THEN
             !Element = MPIGEO%SideToElem(S2E_ELEM_ID,haloSideID)
           !ELSE
             !Element = MPIGEO%SideToElem(S2E_NB_ELEM_ID,haloSideID)
           !END IF
         !ELSE  ! periodic side
           !DoneSideID(1) = haloSideID
           !IF(TriNum.EQ.1) DoneSideID(2) = 2
           !IF(TriNum.EQ.2) DoneSideID(2) = 1
           !DoneLastElem(:,2) = DoneLastElem(:,1)
           !DoneLastElem(1,1) = Element
           !DoneLastElem(2,1) = LocalSide
           !DoneLastElem(3,1) = TriNum
           !PV = abs(MPIGEO%PeriodicElemSide(LocalSide,Element))
           !PVsign = INT(MPIGEO%PeriodicElemSide(LocalSide,Element)/PV)
           !lastPartPos(i,1) = lastPartPos(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
           !lastPartPos(i,2) = lastPartPos(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
           !lastPartPos(i,3) = lastPartPos(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
           !PartState(i,1)   = PartState(i,1) - GEO%PeriodicVectors(1,PV)*PVsign
           !PartState(i,2)   = PartState(i,2) - GEO%PeriodicVectors(2,PV)*PVsign
           !PartState(i,3)   = PartState(i,3) - GEO%PeriodicVectors(3,PV)*PVsign
           !partShiftVector(1:3,i) = partShiftVector(1:3,i) + GEO%PeriodicVectors(1:3,PV)*PVsign
           !IF (MPIGEO%SideToElem(S2E_NB_ELEM_ID,haloSideID).EQ.Element) THEN
             !Element = MPIGEO%SideToElem(S2E_ELEM_ID,haloSideID)
           !ELSE
             !Element = MPIGEO%SideToElem(S2E_NB_ELEM_ID,haloSideID)
           !END IF
         !END IF
       !END IF ! NrOfThroughSides
     !END IF ! InElementCheck
   !END DO ! While NOT DONE
   !RETURN
 !END SUBROUTINE Boundary_Particle_halocells


  !SUBROUTINE MPIParticleSend(n,nshape,indices,shape_indices)                                       !
    !!----------------------------------------------------------------------------------------------!
    !USE MOD_DSMC_Vars,     ONLY : useDSMC, CollisMode, PartStateIntEn, DSMC
    !USE MOD_Particle_Vars, ONLY : PDM,PEM,GEO,PartState,Pt_temp,PartSpecies,LastPartPos,Species,usevMPF,PartMPF
    !USE MOD_part_MPI_Vars, ONLY : MPIGEO, tMPIMessage,PMPIVAR,PMPIExchange,casematrix,NbrOfCases,&
                                  !partShiftVector
    !USE MOD_PICDepo_Vars,  ONLY : DepositionType
    !USE MOD_Globals,       ONLY : myRank
    !use mpi
    !!USE MOD_LD_Vars,       ONLY: UseLD, PartStateBulkValues
    !!----------------------------------------------------------------------------------------------!
    !IMPLICIT NONE                                                                                  !
    !!----------------------------------------------------------------------------------------------!
    !!INCLUDE 'mpif.h'                                                                               !
    !!----------------------------------------------------------------------------------------------!
    !! Argument list declaration                                                                    !
    !INTEGER, INTENT(IN)              :: n, nshape                                                  !
    !INTEGER, INTENT(IN), OPTIONAL    :: indices(:)                                                 !
    !INTEGER, INTENT(IN), OPTIONAL    :: shape_indices(:)                                           !
    !! local variables                                                                              !
    !TYPE(tMPIMessage)                :: send_message(0:PMPIVAR%nProcs-1)                           !
    !INTEGER                          :: send_request(0:PMPIVAR%nProcs-1,1:2)                       !
    !INTEGER                          :: send_status_list(1:MPI_STATUS_SIZE,0:PMPIVAR%nProcs-1,1:2) !
    !INTEGER                          :: counter_phys(0:PMPIVAR%nProcs-1)                                !
    !INTEGER                          :: counter_shape(0:PMPIVAR%nProcs-1)                                !
    !INTEGER                          :: PartMPIprocs(1:PMPIVAR%nProcs),nProcs                      !
    !INTEGER                          :: IERROR, allocStat                                          !
    !INTEGER                          :: nbrOfVariablesPerParticle                                  !
    !INTEGER                          :: i,j,iProc,ind,pos,CellX,CellY,CellZ                        !
    !INTEGER                          :: MsgLength                                                  !
    !INTEGER                          :: nRecv,totalnSendParticles,totalnRecvParticles              !
    !INTEGER                          :: myRealKind                                                 !
    !REAL                             :: myRealTestValue                                            !
    !INTEGER                          :: iCase, CellX2,CellY2,CellZ2
    !REAL                             :: Vec1(1:3), Vec2(1:3), Vec3(1:3), ShiftedPart(1:3)
    !LOGICAL                          :: NOPE
    !!----------------------------------------------------------------------------------------------!
    !!--- allocate message buffers and request buffers
    !ALLOCATE(PMPIExchange%send_message(0:PMPIVAR%nProcs-1),&
             !PMPIExchange%send_request(0:PMPIVAR%nProcs-1,1:2),&
             !PMPIExchange%NbrArray(1:PMPIVAR%nProcs*2),&
             !PMPIExchange%nbrOfSendParticles(0:PMPIVAR%nProcs-1,1:2),STAT=allocStat)
    !IF (allocStat .NE. 0) THEN
      !WRITE(*,*)'ERROR in MPIParticleSend: cannot allocate PMPIExchange variables!'; STOP
    !END IF
    !!--- determine datatype length for variables to be sent
    !myRealKind = KIND(myRealTestValue)
    !IF (myRealKind.EQ.4) THEN
      !myRealKind = MPI_REAL
    !ELSE IF (myRealKind.EQ.8) THEN
      !myRealKind = MPI_DOUBLE_PRECISION
    !ELSE
      !myRealKind = MPI_REAL
    !END IF
    !!--- define number of variables per particle for different cases
    !IF (useDSMC.AND.(CollisMode.NE.1)) THEN
      !IF (usevMPF .AND. DSMC%ElectronicState) THEN
        !nbrOfVariablesPerParticle = 18
      !ELSE IF (usevMPF ) THEN
        !nbrOfVariablesPerParticle = 17
      !ELSE IF ( DSMC%ElectronicState ) THEN
        !nbrOfVariablesPerParticle = 17
      !ELSE
        !nbrOfVariablesPerParticle = 16
      !END IF
    !ELSE
      !IF (usevMPF) THEN
        !nbrOfVariablesPerParticle = 15
      !ELSE
        !nbrOfVariablesPerParticle = 14
      !END IF
    !END IF
    !!--- 5 more PartValues for LD
    !!IF(UseLD) nbrOfVariablesPerParticle = nbrOfVariablesPerParticle + 5
!#if ((PP_TimeDiscMethod!=1) && (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6))  /* RK3 and RK4 only */
    !nbrOfVariablesPerParticle = nbrOfVariablesPerParticle - 6
!#endif
    !!--- determine number of particles sent to different procs
    !PMPIExchange%nbrOfSendParticles(:,:)=0
    !! --------------------------------------------------------------------
!!    IF (GEO%nPeriodicVectors.GT.0) THEN
!!       Vec1(1:3) = 0.
!!       Vec2(1:3) = 0.
!!       Vec3(1:3) = 0.
!!       IF (GEO%nPeriodicVectors.EQ.1) THEN
!!          Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!!       END IF
!!       IF (GEO%nPeriodicVectors.EQ.2) THEN
!!          Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!!          Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!!       END IF
!!       IF (GEO%nPeriodicVectors.EQ.3) THEN
!!          Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!!          Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!!          Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
!!       END IF
!!    END IF
    !DO ind=1,n + nshape
      !! --------------------------------------------------------------------
      !! Das hier ist noch unklar. Die Rolle vom Loadbalancing erscheint noch nicht so
      !! ganz einleuchtend? Ist damit der Prozess beim Restart gemeint? Wahrscheinlich 
      !! nur, wenn sich die anzahl der Prozesse ändert. Muss mit Jonathan abgeklärt werden. 
      !! ! ! Implement Load-Balancing strategy in case of nProcs changes during restart ! ! !
      !!IF (PMPIExchange%MPITags(ind).GT.0) THEN
      !!  !--- in case of load balancing, the particle needs to be sent to a single proc only
      !!  nProcs=1
      !!  PartMPIProcs(1)=PMPIExchange%MPIProcNbr(ind)
      !!ELSE
        !!--- else the procs that the particle has to be sent to are determined according to the BGM
        !IF(ind.LE.n) THEN ! process only pure MPI-parts
          !i=indices(ind)
        !ELSE              ! process only shape particles
          !i=shape_indices(ind-n)
        !END IF
        
        !! communicate the particle back to the own process because it will play
        !! a role for deposition (only do this for a pure MPI part)
        !IF ((Species(PartSpecies(i))%ChargeIC.NE.0).AND.&   ! case 3: mixed physical and shape parts
            !(DepositionType .EQ. 'shape_function')) THEN 
          !CellX = INT((PartState(i,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          !CellY = INT((PartState(i,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          !CellZ = INT((PartState(i,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          !NOPE = .FALSE.
          !IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
               !(CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
               !(CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
            !NOPE = .TRUE.
          !ELSE
            !IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
              !NOPE = .TRUE.
            !END IF
          !END IF
          !IF (NOPE) THEN
             !! it is possible that the particle has been moved over a periodic side 
             !IF (GEO%nPeriodicVectors.GT.0) THEN
               !ShiftedPart(1:3) = PartState(i,1:3) + partShiftVector(1:3,i)
               !CellX = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
               !CellY = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
               !CellZ = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
               !IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
                 !IPWRITE(*,*)'ERROR in MPIParticleSend: Particle outside BGM! Err2'
                 !WRITE(*,*)'i =',i,',ParticleInside =',PDM%ParticleInside(i)
                 !WRITE(*,'(3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMkmin,',minZ =',GEO%FIBGMlmin
                 !WRITE(*,'(3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
                 !WRITE(*,'(3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMkmax,',maxZ =',GEO%FIBGMlmax
                 !WRITE(*,'(3(A,ES13.5))')'PartX=',ShiftedPart(1),',PartY=',ShiftedPart(2),',PartZ=',ShiftedPart(3)
                 !STOP
               !END IF
             !ELSE
                !IPWRITE(*,*)'ERROR in MPIParticleSend: Particle outside BGM!'
                !WRITE(*,*)'i =',i,',ParticleInside =',PDM%ParticleInside(i)
                !WRITE(*,'(3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMkmin,',minZ =',GEO%FIBGMlmin
                !WRITE(*,'(3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
                !WRITE(*,'(3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMkmax,',maxZ =',GEO%FIBGMlmax
                !WRITE(*,'(3(A,ES13.5))')'PartX=',PartState(i,1),',PartY=',PartState(i,2),',PartZ=',PartState(i,3)
                !STOP
             !END IF
          !END IF
          !nProcs=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)
          !PartMPIProcs(1:nProcs)=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(2:nProcs+1)
          !! communicate to our own proc (only if has left our domain: ind<n)
          !IF (ind.LE.n) THEN  ! case 3 only (case 2 else)
            !nProcs=nProcs+1
            !PartMPIProcs(nProcs)=PMPIVAR%iProc
          !END IF
        !ELSE ! case 1: only physical particle (NBC)
          !! neutral particles or non-shape-function particle (NBC)
          !nProcs = 1
          !PartMPIProcs(nProcs) = MPIGEO%ElemMPIID(PEM%Element(i))
        !END IF
      !!END IF ! loadbalancing
      !DO j=1,nProcs
        !iProc=PartMPIProcs(j)
        !IF (ind.LE.n) THEN 
          !IF (iProc.EQ.MPIGEO%ElemMPIID(PEM%Element(i))) THEN ! pure MPI part 
            !PMPIExchange%nbrOfSendParticles(iProc,1)=PMPIExchange%nbrOfSendParticles(iProc,1)+1
          !ELSE             ! shape part (here also the particles are sent back to myproc)
            !PMPIExchange%nbrOfSendParticles(iProc,2)=PMPIExchange%nbrOfSendParticles(iProc,2)+1
          !END IF
        !ELSE               ! shape MPI part
          !PMPIExchange%nbrOfSendParticles(iProc,2)=PMPIExchange%nbrOfSendParticles(iProc,2)+1
        !END IF
      !END DO
    !END DO ! ind
    !!--- Don't send messages to myself unless shape_function!
    !IF ((SUM(PMPIExchange%nbrOfSendParticles(PMPIVAR%iProc,:)).NE.0).AND.&
        !(DepositionType .NE. 'shape_function')) THEN
       !WRITE(*,*) "ERROR: Process:",PMPIVAR%iProc,"sends particles to itself!!!"
       !STOP
    !END IF
    !!IF (DepositionType .NE. 'shape_function') PMPIExchange%nbrOfSendParticles(PMPIVAR%iProc)=0
    !!--- (non-blocking:) send number of particles to all procs
    !!    and allocate send messages
    !!--- build array
    !DO iProc=0,PMPIVAR%nProcs-1
      !PMPIExchange%NbrArray(iProc*2+1) = PMPIExchange%nbrOfSendParticles(iProc,1)
      !PMPIExchange%NbrArray(iProc*2+2) = PMPIExchange%nbrOfSendParticles(iProc,2)
    !END DO
    
    !DO iProc=0,PMPIVAR%nProcs-1
      !!--- Don't send messages to myself!
      !IF ((iProc.EQ.PMPIVAR%iProc).AND.(DepositionType.NE.'shape_function')) CYCLE
      !IF (.NOT.(PMPIVAR%MPINeighbor(iProc))) CYCLE
      !!--- MPI_ISEND lengths of lists of particles leaving local mesh
!!WRITE(*,*) 'send',PMPIVAR%iProc,iProc, PMPIExchange%NbrArray(iProc*2+1:iProc*2+2)
      !CALL MPI_ISEND(PMPIExchange%NbrArray(iProc*2+1:iProc*2+2), 2, MPI_INTEGER, iProc, 1001, PMPIVAR%COMM, &
                        !PMPIExchange%send_request(iProc,1), IERROR)
!!      WRITE(*,'(A,I2.2,A,I5,A,I2,A)')'Rank=',PMPIVAR%iProc,'; Sending ', &
!!                              PMPIExchange%nbrOfSendParticles(iProc),' particles to process ', iProc, ' ...'
      !!--- allocate send messages
      !IF (SUM(PMPIExchange%nbrOfSendParticles(iProc,:)).GT.0) THEN
        !MsgLength = PMPIExchange%nbrOfSendParticles(iProc,1) * nbrOfVariablesPerParticle & 
                   !+PMPIExchange%nbrOfSendParticles(iProc,2) * 7                           ! Pos,Vel,Charge for deposition
        !ALLOCATE( PMPIExchange%send_message(iProc)%content(1:MsgLength), STAT=allocStat )
        !IF (allocStat .NE. 0) THEN
          !WRITE(*,*)'ERROR in MPIParticleSend: cannot allocate PMPIExchange variables (proc ',iProc,')!'; STOP
        !END IF
      !END IF
    !END DO
    !!--- fill send messages
    !counter_phys(:)  = 0
    !counter_shape(:) = PMPIExchange%nbrOfSendParticles(:,1) &
                       !* nbrOfVariablesPerParticle
    !DO ind=1,n+nshape
      !IF(ind.LE.n) THEN ! process only pure MPI-parts
        !i=indices(ind)
      !ELSE              ! process only shapte particles
        !i=shape_indices(ind-n)
      !END IF
      !!IF (PMPIExchange%MPITags(ind).GT.0) THEN
      !!  nProcs=1
      !!  PartMPIProcs(1)=PMPIExchange%MPIProcNbr(ind)
      !!ELSE
        !IF ((Species(PartSpecies(i))%ChargeIC.NE.0).AND.&   ! case 3: mixed physical and shape parts
            !(DepositionType .EQ. 'shape_function')) THEN 
          !CellX = INT((PartState(i,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          !CellY = INT((PartState(i,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          !CellZ = INT((PartState(i,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          !NOPE = .FALSE.
          !IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
               !(CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
               !(CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
            !NOPE = .TRUE.
          !ELSE
            !IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
              !NOPE = .TRUE.
            !END IF
          !END IF
          !IF (NOPE) THEN
!!          IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
!!              (CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
!!              (CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
             !! Particle has to have moved over a periodic side 
             !ShiftedPart(1:3) = PartState(i,1:3) + partShiftVector(1:3,i)
             !CellX = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
             !CellY = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
             !CellZ = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
!!             DO iCase = 1, NbrOfCases
!!                IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
!!                    (casematrix(iCase,2).EQ.0) .AND. &
!!                    (casematrix(iCase,3).EQ.0)) CYCLE
!!                ShiftedPart(1:3) = PartState(i,1:3) + &
!!                     casematrix(iCase,1)*Vec1(1:3) + &
!!                     casematrix(iCase,2)*Vec2(1:3) + &
!!                     casematrix(iCase,3)*Vec3(1:3)
!!                CellX2 = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
!!                CellY2 = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
!!                CellZ2 = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
!!                IF ((CellX2.GT.GEO%FIBGMimax).OR.(CellX2.LT.GEO%FIBGMimin) .OR. &
!!                     (CellY2.GT.GEO%FIBGMkmax).OR.(CellY2.LT.GEO%FIBGMkmin) .OR. &
!!                     (CellZ2.GT.GEO%FIBGMlmax).OR.(CellZ2.LT.GEO%FIBGMlmin)) THEN
!!                ELSE
!!                   CellX = CellX2
!!                   CellY = CellY2
!!                   CellZ = CellZ2
!!                END IF
!!             END DO
          !END IF
          !nProcs=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)
          !PartMPIProcs(1:nProcs)=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(2:nProcs+1)
          !! communicate to our own proc (only if has left our domain: ind<n)
          !IF (ind.LE.n) THEN  ! case 3 only (case 2 else)
            !nProcs=nProcs+1
            !PartMPIProcs(nProcs)=PMPIVAR%iProc
          !END IF
        !ELSE ! case 1: only physical particle (NBC)
          !! neutral particles or non-shape-function particle (NBC)
          !nProcs = 1
          !PartMPIProcs(nProcs) = MPIGEO%ElemMPIID(PEM%Element(i))
        !END IF
      !!END IF
      !DO j=1,nProcs
        !iProc=PartMPIProcs(j)
        !!!--- Don't send messages to myself!
        !!IF ((iProc.EQ.PMPIVAR%iProc).AND.(DepositionType.NE.'shape_function')) CYCLE
        !IF (ind.LE.n) THEN 
          !! at this point the particle can be "deleted". Don't delete PartState since it is
          !! required for deposition. It is not required to set all quantities to zero
          !! since the ParticleInside-Flag is enough information.
          !PDM%ParticleInside(i) = .FALSE.  ! THIS IS IMPORTANT HERE - IT CANNOT BE SET IN BOUNDARY_PIC_PARTICLE!
          !IF (iProc.EQ.MPIGEO%ElemMPIID(PEM%Element(i))) THEN ! pure MPI part 
            !pos = counter_phys(iProc)
            !PMPIExchange%send_message(iProc)%content(pos+1:pos+6) = PartState(i,1:6)   ! to get the positions in the list...
            !PMPIExchange%send_message(iProc)%content(pos+7) = REAL(PartSpecies(i))
!#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* RK3 and RK4 only */
            !PMPIExchange%send_message(iProc)%content(pos+8:pos+13) = Pt_temp(i,1:6)
            !PMPIExchange%send_message(iProc)%content(pos+14) = MPIGEO%NativeElemID(PEM%Element(i))
            !!IF(.NOT.UseLD) THEN   
              !IF (useDSMC.AND.(CollisMode.NE.1)) THEN
                !IF (usevMPF .AND. DSMC%ElectronicState) THEN
                  !PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)    
                  !PMPIExchange%send_message(iProc)%content(pos+17) = PartMPF(i)
                  !PMPIExchange%send_message(iProc)%content(pos+18) = PartStateIntEn(i, 3)
                !ELSE IF (usevMPF) THEN
                  !PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)    
                  !PMPIExchange%send_message(iProc)%content(pos+17) = PartMPF(i)
                !ELSE IF ( DSMC%ElectronicState ) THEN
                  !PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)    
                  !PMPIExchange%send_message(iProc)%content(pos+17) = PartStateIntEn(i, 3)
                !ELSE
                  !PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)
                !END IF
              !ELSE
                !IF (usevMPF) PMPIExchange%send_message(iProc)%content(pos+15) = PartMPF(i)
              !END IF
            !!ELSE ! UseLD == true      =>      useDSMC == true
            !!  IF (CollisMode.NE.1) THEN
            !!    IF (usevMPF .AND. DSMC%ElectronicState) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)    
            !!      PMPIExchange%send_message(iProc)%content(pos+17) = PartMPF(i)
            !!      PMPIExchange%send_message(iProc)%content(pos+18) = PartStateIntEn(i, 3)
            !!      PMPIExchange%send_message(iProc)%content(pos+19:pos+23) = PartStateBulkValues(i,1:5)
            !!    ELSE IF (usevMPF) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)    
            !!      PMPIExchange%send_message(iProc)%content(pos+17) = PartMPF(i)
            !!      PMPIExchange%send_message(iProc)%content(pos+18:pos+22) = PartStateBulkValues(i,1:5)
            !!    ELSE IF ( DSMC%ElectronicState ) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)    
            !!      PMPIExchange%send_message(iProc)%content(pos+17) = PartStateIntEn(i, 3)
            !!      PMPIExchange%send_message(iProc)%content(pos+18:pos+22) = PartStateBulkValues(i,1:5)
            !!    ELSE
            !!      PMPIExchange%send_message(iProc)%content(pos+15) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+16) = PartStateIntEn(i, 2)
            !!      PMPIExchange%send_message(iProc)%content(pos+17:pos+21) = PartStateBulkValues(i,1:5)
            !!    END IF
            !!  ELSE
            !!    IF (usevMPF) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+15) = PartMPF(i)
            !!      PMPIExchange%send_message(iProc)%content(pos+16:pos+20) = PartStateBulkValues(i,1:5)
            !!    ELSE
            !!      PMPIExchange%send_message(iProc)%content(pos+15:pos+19) = PartStateBulkValues(i,1:5)
            !!    END IF
            !!  END IF
            !!END IF
!#else 
            !PMPIExchange%send_message(iProc)%content(pos+8) = MPIGEO%NativeElemID(PEM%Element(i))
            !!IF(.NOT.UseLD) THEN      
              !IF (useDSMC.AND.(CollisMode.NE.1)) THEN
                !IF (usevMPF .AND. DSMC%ElectronicState) THEN
                  !PMPIExchange%send_message(iProc)%content(pos+ 9) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)    
                  !PMPIExchange%send_message(iProc)%content(pos+11) = PartMPF(i)
                  !PMPIExchange%send_message(iProc)%content(pos+12) = PartStateIntEn(i, 3)
                !ELSE IF (usevMPF) THEN
                  !PMPIExchange%send_message(iProc)%content(pos+ 9) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)    
                  !PMPIExchange%send_message(iProc)%content(pos+11) = PartMPF(i)
                !ELSE IF ( DSMC%ElectronicState ) THEN
                  !PMPIExchange%send_message(iProc)%content(pos+ 9) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)    
                  !PMPIExchange%send_message(iProc)%content(pos+11) = PartStateIntEn(i, 3)
                !ELSE
                  !PMPIExchange%send_message(iProc)%content(pos+ 9) = PartStateIntEn(i, 1)
                  !PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)
                !END IF
              !ELSE
                !IF (usevMPF) PMPIExchange%send_message(iProc)%content(pos+9) = PartMPF(i)
              !END IF
            !!ELSE ! UseLD == true      =>      useDSMC == true
            !!  IF (CollisMode.NE.1) THEN
            !!    IF (usevMPF .AND. DSMC%ElectronicState) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+9) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)    
            !!      PMPIExchange%send_message(iProc)%content(pos+11) = PartMPF(i)
            !!      PMPIExchange%send_message(iProc)%content(pos+12) = PartStateIntEn(i, 3)
            !!      PMPIExchange%send_message(iProc)%content(pos+13:pos+17) = PartStateBulkValues(i,1:5)
            !!    ELSE IF (usevMPF) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+9) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)    
            !!      PMPIExchange%send_message(iProc)%content(pos+11) = PartMPF(i)
            !!      PMPIExchange%send_message(iProc)%content(pos+12:pos+16) = PartStateBulkValues(i,1:5)
            !!    ELSE IF ( DSMC%ElectronicState ) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+9) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)    
            !!      PMPIExchange%send_message(iProc)%content(pos+11) = PartStateIntEn(i, 3)
            !!      PMPIExchange%send_message(iProc)%content(pos+12:pos+16) = PartStateBulkValues(i,1:5)
            !!    ELSE
            !!      PMPIExchange%send_message(iProc)%content(pos+9) = PartStateIntEn(i, 1)
            !!      PMPIExchange%send_message(iProc)%content(pos+10) = PartStateIntEn(i, 2)
            !!      PMPIExchange%send_message(iProc)%content(pos+11:pos+15) = PartStateBulkValues(i,1:5)
            !!    END IF
            !!  ELSE
            !!    IF (usevMPF) THEN
            !!      PMPIExchange%send_message(iProc)%content(pos+9) = PartMPF(i)
            !!      PMPIExchange%send_message(iProc)%content(pos+10:pos+14) = PartStateBulkValues(i,1:5)
            !!    ELSE
            !!      PMPIExchange%send_message(iProc)%content(pos+9:pos+13) = PartStateBulkValues(i,1:5)
            !!    END IF
            !!  END IF
            !!END IF
!#endif 
            !counter_phys(iProc) = counter_phys(iProc) + nbrOfVariablesPerParticle   
          !ELSE             ! shape part (here also the particles are sent back to myproc)
            !pos = counter_shape(iProc)
            !PMPIExchange%send_message(iProc)%content(pos+1:pos+6) = PartState(i,1:6)   ! to get the positions in the list...
            !PMPIExchange%send_message(iProc)%content(pos+7) = REAL(PartSpecies(i))
            !counter_shape(iProc) = counter_shape(iProc) + 7                
          !END IF
        !ELSE               ! shape MPI part
          !pos = counter_shape(iProc)
          !PMPIExchange%send_message(iProc)%content(pos+1:pos+6) = PartState(i,1:6)   ! to get the positions in the list...
          !PMPIExchange%send_message(iProc)%content(pos+7) = REAL(PartSpecies(i))
          !counter_shape(iProc) = counter_shape(iProc) + 7                
        !END IF
      !END DO  ! nProcs
    !END DO ! ind

    !!--- (non-blocking:) send messages to all procs receiving particles from myself
    !DO iProc = 0,PMPIVAR%nProcs-1
      !IF (SUM(PMPIExchange%nbrOfSendParticles(iProc,:)) .GT. 0) THEN
        !!--- MPI_ISEND particle data of particles leaving local mesh
        !MsgLength = PMPIExchange%nbrOfSendParticles(iProc,1) * nbrOfVariablesPerParticle & 
                   !+PMPIExchange%nbrOfSendParticles(iProc,2) * 7                           ! Pos,Vel,Charge for deposition
        !CALL MPI_ISEND(PMPIExchange%send_message(iProc)%content,MsgLength,myRealKind,iProc,1002,PMPIVAR%COMM, &
                          !PMPIExchange%send_request(iProc,2), IERROR)
!!        WRITE(*,'(A,I2.2,A,I5,A,I2,A)')'Rank=',PMPIVAR%iProc,'; Sending ', &
!!                                PMPIExchange%nbrOfSendParticles(iProc),' particles to process ', iProc, '!'
      !END IF
    !END DO
  !END SUBROUTINE MPIParticleSend

  !SUBROUTINE MPIParticleRecv(nRecv)                                                                !
    !!----------------------------------------------------------------------------------------------!
    !USE MOD_DSMC_Vars,     ONLY : useDSMC, CollisMode, PartStateIntEn, DSMC
    !USE MOD_Particle_Vars, ONLY : PDM,PEM,GEO,PartState,Pt_temp,PartSpecies,LastPartPos,usevMPF,PartMPF
    !USE MOD_part_MPI_Vars, ONLY : tMPIMessage,PMPIVAR
    !USE MOD_part_MPI_Vars, ONLY : extPartState,extPartSpecies,extPartsAllocated,NbrOfAllocatedExtParts,NbrOfextParticles
    !USE MOD_PICDepo_Vars,  ONLY : DepositionType
    !use mpi
    !!USE MOD_LD_Vars,       ONLY: UseLD, PartStateBulkValues
    !!----------------------------------------------------------------------------------------------!
    !IMPLICIT NONE                                                                                  !
    !!----------------------------------------------------------------------------------------------!
    !!INCLUDE 'mpif.h'                                                                               !
    !!----------------------------------------------------------------------------------------------!
    !! Argument list declaration                                                                    !
    !INTEGER, INTENT(OUT)             :: nRecv                                                      !
    !! local variables                                                                              !
    !INTEGER                          :: nbrOfRecvParticles(0:PMPIVAR%nProcs-1,1:2)                 !
    !TYPE(tMPIMessage)                :: recv_message(0:PMPIVAR%nProcs-1)                           !
    !INTEGER                          :: recv_request(0:PMPIVAR%nProcs-1,1:2)                       !
    !INTEGER                          :: recv_status_list(1:MPI_STATUS_SIZE,0:PMPIVAR%nProcs-1,1:2) !
    !INTEGER                          :: nbrOfVariablesPerParticle                                  !
    !INTEGER                          :: IERROR, allocStat                                          !
    !INTEGER                          :: i,iProc,pos,counter1,counter2                              !
    !INTEGER                          :: dummy1,dummy2                                              !
    !INTEGER                          :: ParticleIndexNbr                                           !
    !INTEGER                          :: MsgLength                                                  !
    !INTEGER                          :: myRealKind                                                 !
    !REAL                             :: myRealTestValue                                            !
    !INTEGER                          :: NbrOfParticlesToReceive(1:2)                               !
    !!----------------------------------------------------------------------------------------------!
    !!--- determine datatype length for variables to be sent
    !myRealKind = KIND(myRealTestValue)
    !IF (myRealKind.EQ.4) THEN
      !myRealKind = MPI_REAL
    !ELSE IF (myRealKind.EQ.8) THEN
      !myRealKind = MPI_DOUBLE_PRECISION
    !ELSE
      !myRealKind = MPI_REAL
    !END IF
    !!--- define number of variables per particle for different cases
    !IF (useDSMC.AND.(CollisMode.NE.1)) THEN
      !IF (usevMPF .AND. DSMC%ElectronicState) THEN
        !nbrOfVariablesPerParticle = 18
      !ELSE IF (usevMPF .OR. DSMC%ElectronicState ) THEN
        !nbrOfVariablesPerParticle = 17
      !ELSE
        !nbrOfVariablesPerParticle = 16
      !END IF
    !ELSE
      !IF (usevMPF) THEN
        !nbrOfVariablesPerParticle = 15
      !ELSE
        !nbrOfVariablesPerParticle = 14
      !END IF
    !END IF
    !!--- 5 more PartValues for LD
    !!IF(UseLD) nbrOfVariablesPerParticle = nbrOfVariablesPerParticle + 5
!#if ((PP_TimeDiscMethod!=1) && (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6))  /* RK3 and RK4 only */
    !nbrOfVariablesPerParticle = nbrOfVariablesPerParticle - 6
!#endif
    !nbrOfRecvParticles(:,:) = 0
    !!--- (non-blocking:) post receive of number of particles from all procs
    !DO iProc=0,PMPIVAR%nProcs-1
      !NbrOfParticlesToReceive(:) = 0
      !!--- Don't send messages to myself unless shape_function!
      !IF ((iProc.EQ.PMPIVAR%iProc).AND.(DepositionType.NE.'shape_function')) CYCLE
      !IF (.NOT.(PMPIVAR%MPINeighbor(iProc))) CYCLE
      !CALL MPI_IRECV(NbrOfParticlesToReceive(:), 2, MPI_INTEGER, iProc, 1001, PMPIVAR%COMM, &
                        !recv_request(iProc,1), IERROR)
      !CALL MPI_WAIT(recv_request(iProc,1),recv_status_list(:,iProc,1),IERROR)
!!IF(PMPIVAR%iProc.EQ.iProc)THEN
!!WRITE(*,*) 'recv',PMPIVAR%iProc,iProc, NbrOfParticlesToReceive
!!END IF
      !nbrOfRecvParticles(iProc,1:2) = NbrOfParticlesToReceive(1:2)
    !END DO

    !!--- To be migrated / placed elsewhere?:
    !nRecv=0
    !counter1=0; counter2=0
!!    DO iProc=0,PMPIVAR%nProcs-1
!!      !--- Don't send messages to myself!
!!      IF ((iProc.EQ.PMPIVAR%iProc).AND.(DepositionType.NE.'shape_function')) CYCLE
!!      !--- wait for number of recv particles from proc iProc to be communicated
!!      CALL MPI_WAIT(recv_request(iProc,1),recv_status_list(:,iProc,1),IERROR)
!!!      WRITE(*,'(A,I2.2,A,I5,A,I2,A)')'Rank=',PMPIVAR%iProc,'; Receiving ', &
!!!                              nbrOfRecvParticles(iProc),' particles from process ', iProc, ' ...'
!!    END DO
    !DO iProc=0,PMPIVAR%nProcs-1
      !!--- if any particles are transmitted from proc iProc, ...
      !IF (SUM(nbrOfRecvParticles(iProc,:)).GT.0) THEN
        !!--- allocate recv buffer for particles communicated from proc iProc
        !MsgLength = nbrOfRecvParticles(iProc,1) * nbrOfVariablesPerParticle & 
                   !+nbrOfRecvParticles(iProc,2) * 7                           ! Pos,Vel,Charge for deposition
        !ALLOCATE( recv_message(iProc)%content(1:MsgLength), STAT=allocStat )
        !IF (allocStat .NE. 0) THEN
          !WRITE(*,*)'ERROR in MPIParticleRecv: cannot allocate recv_message (proc ',iProc,')!'; STOP
        !END IF
        !!--- allocate external particles for shape_function
        !IF ((DepositionType .EQ. 'shape_function') .AND. (extPartsAllocated .EQV. .FALSE.)) THEN
          !NbrOfAllocatedextParts=SUM(NbrOfRecvParticles(:,2))
          !ALLOCATE(extPartState(1:NbrOfAllocatedextParts,1:6), &
                   !extPartSpecies(1:NbrOfAllocatedextParts),   &
                   !STAT=allocStat)
          !IF (allocStat .NE. 0) THEN
            !WRITE(*,*)'ERROR in MPIParticleRecv: cannot allocate extPart... variables!',PMPIVAR%iProc; STOP
          !END IF
          !extPartsAllocated = .TRUE.
        !END IF
        !!--- (non-blocking:) post receive for particles to be expected from neighboring proc iProc
        !CALL MPI_IRECV(recv_message(iProc)%content,MsgLength,myRealKind,iProc,1002,PMPIVAR%COMM, &
                          !recv_request(iProc,2), IERROR)
        !!--- wait until particle data are received
        !!    I think it doesn't really matter that this MPI_WAIT directly follows its MPI_IRECV
        !!    since this message (MPITAG 2) is sent directly after the previous one (MPITAG 1)
        !!    for which we have already waited, i.e. which has already arrived.
        !!    Therefore it is probably best to just wait for the next message at this point.
        !CALL MPI_WAIT(recv_request(iProc,2),recv_status_list(:,iProc,2),IERROR)
!!        WRITE(*,'(A,I2.2,A,I5,A,I2,A)')'Rank=',PMPIVAR%iProc,'; Received ', &
!!                                nbrOfRecvParticles(iProc),' particles from process ', iProc, '!'
        !!--- check if particles are inside my subdomain
        !!    if so, insert particle data into respective arrays
        !pos = 0
        !DO i = 1,SUM(nbrOfRecvParticles(iProc,1:2))
          !! Only store full information if pure MPI part and if not a message 
          !! to myproc (these parts are automatically shape parts)
          !IF(i.LE.nbrOfRecvParticles(iProc,1)) THEN
            !nRecv=nRecv+1
            !ParticleIndexNbr = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
            !IF (ParticleIndexNbr .ne. 0) THEN
              !PartState(ParticleIndexNbr,1:6)   = recv_message(iProc)%content(pos+1:pos+6)! to get the positions in the list...
              !PartSpecies(ParticleIndexNbr)     = INT(recv_message(iProc)%content(pos+7))
!#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* RK3 and RK4 only */
              !Pt_temp(ParticleIndexNbr,1:6)     = recv_message(iProc)%content(pos+8:pos+13)
              !PEM%Element(ParticleIndexNbr)     = recv_message(iProc)%content(pos+14)
              !!IF(.NOT.UseLD) THEN
                !IF (useDSMC.AND.(CollisMode.NE.1)) THEN
                  !IF (usevMPF .AND. DSMC%ElectronicState) THEN
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+17)
                    !PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+18)
                  !ELSE IF ( usevMPF) THEN
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+17)
                  !ELSE IF ( DSMC%ElectronicState ) THEN
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+17)
                  !ELSE
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                  !END IF
                !ELSE
                  !IF (usevMPF) PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+15)
                !END IF
              !!ELSE ! UseLD == true      =>      useDSMC == true
                !!IF (CollisMode.NE.1) THEN
                  !!IF (usevMPF .AND. DSMC%ElectronicState) THEN
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !!PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+17)
                    !!PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+18)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+19:pos+23)
                  !!ELSE IF ( usevMPF) THEN
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !!PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+17)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+18:pos+22)
                  !!ELSE IF ( DSMC%ElectronicState ) THEN
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !!PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+17)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+18:pos+22)
                  !!ELSE
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+15)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+16)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+17:pos+21)
                  !!END IF
                !!ELSE
                  !!IF (usevMPF) THEN
                    !!PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+15)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+16:pos+20)
                  !!ELSE
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+15:pos+19)                    
                  !!END IF
                !!END IF
              !!END IF
!#else 
              !PEM%Element(ParticleIndexNbr)     = recv_message(iProc)%content(pos+8)
              !!IF(.NOT.UseLD) THEN
                !IF (useDSMC.AND.(CollisMode.NE.1)) THEN
                  !IF (usevMPF .AND. DSMC%ElectronicState) THEN
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+11)
                    !PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+12)
                  !ELSE IF ( usevMPF ) THEN
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+11)
                  !ELSE IF ( DSMC%ElectronicState) THEN
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+11)
                  !ELSE
                    !PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                  !END IF
                !ELSE
                  !IF (usevMPF) PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+9)
                !END IF
              !!ELSE ! UseLD == true      =>      useDSMC == true
                !!IF (CollisMode.NE.1) THEN
                  !!IF (usevMPF .AND. DSMC%ElectronicState) THEN
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !!PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+11)
                    !!PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+12)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+13:pos+17)
                  !!ELSE IF ( usevMPF) THEN
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !!PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+11)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+12:pos+16)
                  !!ELSE IF ( DSMC%ElectronicState ) THEN
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !!PartStateIntEn(ParticleIndexNbr, 3) = recv_message(iProc)%content(pos+11)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+12:pos+16)
                  !!ELSE
                    !!PartStateIntEn(ParticleIndexNbr, 1) = recv_message(iProc)%content(pos+9)
                    !!PartStateIntEn(ParticleIndexNbr, 2) = recv_message(iProc)%content(pos+10)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+11:pos+15)
                  !!END IF
                !!ELSE
                  !!IF (usevMPF) THEN
                    !!PartMPF(ParticleIndexNbr) = recv_message(iProc)%content(pos+9)
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+10:pos+14)
                  !!ELSE
                    !!PartStateBulkValues(ParticleIndexNbr,1:5) = recv_message(iProc)%content(pos+9:pos+13)                    
                  !!END IF
                !!END IF
              !!END IF
!#endif
              !! Set Flag for received parts in order to localize them later
              !PEM%lastElement(ParticleIndexNbr) = -888 
              !PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
!!              IF (recv_message(iProc)%content(pos+11).GT.0) THEN  ! load balancing, no need for doing new localization
!!                PEM%Element(ParticleIndexNbr)     = INT(recv_message(iProc)%content(pos+11))
!!                PEM%lastElement(ParticleIndexNbr) = INT(recv_message(iProc)%content(pos+11))
!!              ELSE                                               ! "normal" particle send through MPI side, localization needed!
!!                PEM%lastElement(ParticleIndexNbr) = -888 ! default
!!                !--- sum up only particles that actually are located in my subdomain
!!                nRecv=nRecv-1                                                    ! overwrite this entry
                 !! Shape-Function stuff: copied to below
!!              END IF
            !ELSE
              !WRITE(*,*)'ERROR in ParticleExchange_parallel:'
              !WRITE(*,'(A,I6,A)')'corrupted list: PIC%nextFreePosition(', nRecv, ') = 0!'
              !STOP
            !END IF
            !pos = pos + nbrOfVariablesPerParticle                                        ! since we're sending a 1-DIM array
          !ELSE ! shape part
            !NbrOfextParticles=NbrOfextParticles+1                                        ! set to zero in ParticleBoundary
            !ExtPartState(NbrOfextParticles,1:6) = recv_message(iProc)%content(pos+1:pos+6)
            !ExtPartSpecies(NbrOfextParticles)   = INT(recv_message(iProc)%content(pos+7))
            !pos = pos + 7
          !END IF ! pure MPI parts
        !END DO
        !!--- be nice: deallocate the receive buffer
        !DEALLOCATE( recv_message(iProc)%content )
      !END IF
    !END DO
    !PDM%ParticleVecLength = PDM%ParticleVecLength + SUM(nbrOfRecvParticles(:,1))
    !PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + SUM(nbrOfRecvParticles(:,1))
   !!--- error checking
!!   CALL MPI_ALLREDUCE(counter1,dummy1,1,MPI_INTEGER,MPI_SUM,PMPIVAR%COMM,IERROR)
!!   CALL MPI_ALLREDUCE(counter2,dummy2,1,MPI_INTEGER,MPI_SUM,PMPIVAR%COMM,IERROR)
!!   IF ((PMPIVAR%iProc.EQ.0)) THEN
!!     WRITE(*,'(A,I0,A,I0)')'tag method: ',dummy2,', search method: ',dummy1
!!   END IF
  !END SUBROUTINE MPIParticleRecv
!#endif

 !FUNCTION beta(z,w)                                                                                                
   !USE nr
   !IMPLICIT NONE
   !REAL beta, w, z                                                                                                  
   !beta = exp(gammln(z)+gammln(w)-gammln(z+w))                                                                    
 !END FUNCTION beta 

END MODULE MOD_part_boundary
