#include "boltzplatz.h"

MODULE MOD_ParticleInit
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE InitParticleGeometry
  MODULE PROCEDURE InitParticleGeometry
END INTERFACE

INTERFACE InitElemVolumes
  MODULE PROCEDURE InitElemVolumes
END INTERFACE


PUBLIC::InitParticles,InitParticleGeometry,InitElemVolumes
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
! Subroutine for particle initialization 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_Mesh_Vars,          ONLY : nElems,nNodes
USE MOD_Mesh_Vars,          ONLY : Elems,offsetElem,ElemToSide
USE MOD_Particle_Vars,      ONLY : GEO
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, iLocSide, iNode, jNode, GlobSideID,iCoord
INTEGER           :: nStart, NodeNum
INTEGER           :: ALLOCSTAT
INTEGER           :: NodeMap(4,6)
REAL              :: A(3,3),detcon,Elem_SP(3)
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE GEOMETRY INFORMATION...'
NodeMap(:,1)=(/1,4,3,2/)
NodeMap(:,2)=(/1,2,6,5/)
NodeMap(:,3)=(/2,3,7,6/)
NodeMap(:,4)=(/3,4,8,7/)
NodeMap(:,5)=(/1,5,8,4/)
NodeMap(:,6)=(/5,6,7,8/)
ALLOCATE(GEO%ElemToNodeID(1:8,1:nElems),       &
         GEO%ElemSideNodeID(1:4,1:6,1:nElems), &
         GEO%NodeCoords(1:3,1:nNodes),         &
         GEO%ConcaveElemSide(1:6,1:nElems),    &
         GEO%PeriodicElemSide(1:6,1:nElems),    STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(__STAMP__&
 ,'ERROR in InitParticleGeometry: Cannot allocate GEO%... stuff!')
END IF
GEO%ElemToNodeID(:,:)=0
GEO%ElemSideNodeID(:,:,:)=0
GEO%NodeCoords(:,:)=0.
GEO%PeriodicElemSide(:,:)=0
GEO%ConcaveElemSide(:,:)=.FALSE.
iNode=0
DO iElem=1,nElems
  DO jNode=1,8
    Elems(iElem+offsetElem)%ep%node(jNode)%np%ind=0
  END DO
END DO
DO iElem=1,nElems
  !--- Save corners of sides
  DO jNode=1,8
    IF (Elems(iElem+offsetElem)%ep%node(jNode)%np%ind.EQ.0) THEN
      iNode=iNode+1
      Elems(iElem+offsetElem)%ep%node(jNode)%np%ind=iNode
      GEO%NodeCoords(1:3,iNode)=Elems(iElem+offsetElem)%ep%node(jNode)%np%x(1:3)
    END IF
    GEO%ElemToNodeID(jNode,iElem)=Elems(iElem+offsetElem)%ep%node(jNode)%np%ind
  END DO
END DO
DO iElem=1,nElems
  DO iLocSide=1,6
    nStart=MAX(0,ElemToSide(E2S_FLIP,iLocSide,iElem)-1)
    GEO%ElemSideNodeID(1:4,iLocSide,iElem)=(/Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart  ,4)+1,iLocSide))%np%ind,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+1,4)+1,iLocSide))%np%ind,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+2,4)+1,iLocSide))%np%ind,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+3,4)+1,iLocSide))%np%ind/)
  END DO
END DO
!--- Initialize Periodic Side Info
DO iElem=1,nElems
  DO iLocSide=1,6
    GlobSideID = ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    IF ((Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%BCindex.GT.0)           .AND. &
        (ASSOCIATED(Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%connection))) THEN
      GEO%PeriodicElemSide(iLocSide,iElem) = -1
    END IF
  END DO
END DO
!--- Save whether Side is concave or convex
DO iElem = 1,nElems
  DO iLocSide = 1,6
    !--- Check whether the bilinear side is concave
    !--- Node Number 4 and triangle 1-2-3
    DO NodeNum = 1,3               ! for all 3 nodes of triangle
      A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,iElem)) &
                   - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,iLocSide,iElem))
    END DO
    !--- concave if detcon < 0:
    detcon = ((A(2,1) * A(3,2) - A(3,1) * A(2,2)) * A(1,3) +     &
              (A(3,1) * A(1,2) - A(1,1) * A(3,2)) * A(2,3) +     &
              (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * A(3,3))
    IF (detcon.LT.0) GEO%ConcaveElemSide(iLocSide,iElem)=.TRUE.
  END DO
END DO
! for supersampeledsurfaces see particle_surfases, sub: GetSuperSampledSurfaces
!--- check for elements with intersecting sides (e.g. very flat elements)
CALL WeirdElementCheck()

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GEOMETRY INFORMATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleGeometry

SUBROUTINE WeirdElementCheck()
!===================================================================================================================================
! Calculate whether element edges intersect other sides
! If this is the case it means that part of the element is turned inside-out
! which results in a warning so the user can decide whether it is a problem that 
! necessitates a new mesh. 
! Fixing the problem would involve defining the bilinear edge between nodes 2 and 4
! (instead of 1 and 3). This information would need to be stored and used throughout
! the particle treatment. Additionally, since the edge would need to be changed 
! for both neighboring elements, it is possible that both element might have the problem
! hence no solution exists.
! tl;dr: Hard/maybe impossible to fix, hence only a warning is given so the user can decide
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_Mesh_Vars,          ONLY : nElems
USE MOD_Particle_Vars,      ONLY : GEO, WeirdElems
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, iLocSide, kLocSide, iNode, WeirdElemNbrs(1:nElems)
REAL              :: vec(1:3), Node(1:3,1:4),det(1:3)
LOGICAL           :: WEIRD, TRICHECK, TRIABSCHECK
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' CHECKING FOR WEIRD ELEMENTS...'

WeirdElems = 0
DO iElem = 1, nElems ! go through all elements
  WEIRD = .FALSE.
  DO iLocSide = 1,5  ! go through local sides
    IF (.not.WEIRD) THEN  ! if one is found there is no need to continue
      IF (GEO%ConcaveElemSide(iLocSide,iElem)) THEN  ! only concave elements need to be checked
        ! build vector from node 1 to node 3
        vec(:) = GEO%NodeCoords(:,GEO%ElemSideNodeID(3,iLocSide,iElem)) &
               - GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,iElem))
        ! check all other sides
        DO kLocSide = iLocSide + 1, 6
          IF (GEO%ConcaveElemSide(kLocSide,iElem)) THEN  ! only concave elements need to be checked
            ! build 4 vectors from point 1 of edge to 4 nodes of kLocSide
            DO iNode = 1,4
              Node(:,iNode) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,iElem)) &
                            - GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,kLocSide,iElem))
            END DO
            ! Compute whether any of the triangle intersects with the vector vec:
            ! If all three volumes built by the vector vec and the vectors Node
            ! are either positive or negative then there is an intersection

            ! Triangle 1 (Nodes 1,2,3)
            ! Only check this if neither point of vec is part of the triangle.
            ! If points of vec correspont to point 1 or 3 or triangle then both
            ! triangles can be skipped (triabscheck = true), else point 4 needs to be checked
            ! separately for triangle 2 (see below)
            TRICHECK = .FALSE.
            TRIABSCHECK = .FALSE.
            DO iNode = 1,3
              det(:) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,iElem)) &
                     - GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,kLocSide,iElem))
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
              det(:) = GEO%NodeCoords(:,GEO%ElemSideNodeID(3,iLocSide,iElem)) &
                     - GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,kLocSide,iElem))
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
            END DO
            IF (.not.TRICHECK) THEN
              det(1) = ((Node(2,1) * Node(3,2) - Node(3,1) * Node(2,2)) * vec(1)  + &
                        (Node(3,1) * Node(1,2) - Node(1,1) * Node(3,2)) * vec(2)  + & 
                        (Node(1,1) * Node(2,2) - Node(2,1) * Node(1,2)) * vec(3))
              det(2) = ((Node(2,2) * Node(3,3) - Node(3,2) * Node(2,3)) * vec(1)  + &
                        (Node(3,2) * Node(1,3) - Node(1,2) * Node(3,3)) * vec(2)  + & 
                        (Node(1,2) * Node(2,3) - Node(2,2) * Node(1,3)) * vec(3))
              det(3) = ((Node(2,3) * Node(3,1) - Node(3,3) * Node(2,1)) * vec(1)  + &
                        (Node(3,3) * Node(1,1) - Node(1,3) * Node(3,1)) * vec(2)  + & 
                        (Node(1,3) * Node(2,1) - Node(2,3) * Node(1,1)) * vec(3))
              IF ((det(1).LT.0).AND.(det(2).LT.0).AND.(det(3).LT.0)) WEIRD = .TRUE.
              IF ((det(1).GT.0).AND.(det(2).GT.0).AND.(det(3).GT.0)) WEIRD = .TRUE.
            END IF

            ! Triangle 2 (Nodes 1,3,4)
            TRICHECK = .FALSE.
            IF (.not.TRIABSCHECK) THEN
              ! Node 4 needs to be checked separately (see above)
              det(:) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,iElem)) &
                     - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,kLocSide,iElem))
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              det(:) = GEO%NodeCoords(:,GEO%ElemSideNodeID(3,iLocSide,iElem)) &
                     - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,kLocSide,iElem))
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              IF (.not.TRICHECK) THEN
                det(1) = ((Node(2,1) * Node(3,3) - Node(3,1) * Node(2,3)) * vec(1)  + &
                          (Node(3,1) * Node(1,3) - Node(1,1) * Node(3,3)) * vec(2)  + & 
                          (Node(1,1) * Node(2,3) - Node(2,1) * Node(1,3)) * vec(3))
                det(2) = ((Node(2,3) * Node(3,4) - Node(3,3) * Node(2,4)) * vec(1)  + &
                          (Node(3,3) * Node(1,4) - Node(1,3) * Node(3,4)) * vec(2)  + & 
                          (Node(1,3) * Node(2,4) - Node(2,3) * Node(1,4)) * vec(3))
                det(3) = ((Node(2,4) * Node(3,1) - Node(3,4) * Node(2,1)) * vec(1)  + &
                          (Node(3,4) * Node(1,1) - Node(1,4) * Node(3,1)) * vec(2)  + & 
                          (Node(1,4) * Node(2,1) - Node(2,4) * Node(1,1)) * vec(3))
                IF ((det(1).LT.0).AND.(det(2).LT.0).AND.(det(3).LT.0)) WEIRD = .TRUE.
                IF ((det(1).GT.0).AND.(det(2).GT.0).AND.(det(3).GT.0)) WEIRD = .TRUE.
              END IF
            END IF
          END IF
        END DO
      END IF
    END IF
  END DO
  IF (WEIRD) THEN 
    WeirdElems = WeirdElems + 1
    WeirdElemNbrs(WeirdElems) = iElem
  END IF
END DO
              
SWRITE(UNIT_stdOut,'(A)')' CHECKING FOR WEIRD ELEMENTS DONE!'
SWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
SWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
DO iElem = 1,WeirdElems
  SWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
!  DO iLocSide = 1,6
!    DO iNode = 1,4
!      SWRITE(UNIT_stdOut,*) GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,iLocSide,WeirdElemNbrs(iElem)))
!    END DO
!  END DO
!  STOP
END DO
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WeirdElementCheck

SUBROUTINE InitElemVolumes()
!===================================================================================================================================
! Calculate Element volumes for later use in particle routines
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_Mesh_Vars,          ONLY : nElems,NGeo,sJ
USE MOD_Interpolation_Vars, ONLY : wGP
USE MOD_Particle_Vars,      ONLY : GEO, usevMPF
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem
INTEGER           :: i,j,k
INTEGER           :: ALLOCSTAT
REAL              :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE GEOMETRY INFORMATION (Element Volumes)...'
ALLOCATE(GEO%Volume(nElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR in InitParticleGeometry: Cannot allocate GEO%Volume!')
END IF
usevMPF = GETLOGICAL('Part-vMPF','.FALSE.')
IF(usevMPF) THEN
  ALLOCATE(GEO%DeltaEvMPF(nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    ,'ERROR in InitParticleGeometry: Cannot allocate GEO%DeltaEvMPF!')
  END IF
  GEO%DeltaEvMPF(:) = 0.0
END IF
DO iElem=1,nElems
  !--- Calculate and save volume of element iElem
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  GEO%Volume(iElem) = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    GEO%Volume(iElem) = GEO%Volume(iElem) + wGP(i)*wGP(j)*wGP(k)*J_N(1,i,j,k)
  END DO; END DO; END DO
END DO
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GEOMETRY INFORMATION (Element Volumes) DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitElemVolumes

SUBROUTINE InitParticles()
!===================================================================================================================================
! Glue Subroutine for particle initialization 
!===================================================================================================================================
! MODULES
USE MOD_Globals!,       ONLY: MPIRoot,UNIT_STDOUT
USE MOD_ReadInTools
USE MOD_Particle_Vars, ONLY: ParticlesInitIsDone, WriteMacroValues, nSpecies
USE MOD_part_emission, ONLY: InitializeParticleEmission
USE MOD_DSMC_Init,     ONLY: InitDSMC
!USE MOD_LD_Init,       ONLY: InitLD
!USE MOD_LD_Vars,       ONLY: useLD
USE MOD_DSMC_Vars,     ONLY: useDSMC, DSMC, SampDSMC
USE MOD_Mesh_Vars,     ONLY : nElems
USE MOD_InitializeBackgroundField
USE MOD_PICInterpolation_Vars, ONLY: useBGField
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef MPI
#endif
!===================================================================================================================================
IF(ParticlesInitIsDone)THEN
   SWRITE(*,*) "InitParticles already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES ...'

CALL InitializeVariables()
IF(useBGField) CALL InitializeBackgroundField()
CALL InitializeParticleEmission()

IF (useDSMC) THEN
  CALL  InitDSMC()
  !IF (useLD) CALL InitLD
ELSE IF (WriteMacroValues) THEN
  DSMC%SampNum = 0
  ALLOCATE(SampDSMC(nElems,nSpecies))
  SampDSMC(1:nElems,1:nSpecies)%PartV(1)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV(2)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV(3)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(1) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(2) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(3) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartNum   = 0
  SampDSMC(1:nElems,1:nSpecies)%SimPartNum   = 0
  SampDSMC(1:nElems,1:nSpecies)%ERot      = 0
  SampDSMC(1:nElems,1:nSpecies)%EVib      = 0
END IF

ParticlesInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticles


SUBROUTINE InitializeVariables()
!===================================================================================================================================
! Initialize the variables first 
!===================================================================================================================================
! MODULES
USE MOD_Globals!, ONLY:MPIRoot,UNIT_STDOUT,myRank,nProcessors
USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_Particle_Vars!, ONLY: 
USE MOD_Mesh_Vars,             ONLY: nElems, BoundaryName,BoundaryType, nBCs
USE MOD_Restart_Vars,          ONLY: DoRestart
USE MOD_DSMC_Vars,             ONLY: useDSMC
USE MOD_Particle_Output_Vars,  ONLY: WriteFieldsToVTK, OutputMesh
USE MOD_part_MPFtools,         ONLY: DefinePolyVec, DefineSplitVec
USE MOD_PICInterpolation_Vars, ONLY: InterpolationType
USE MOD_PICInterpolation,      ONLY: InitializeInterpolation
USE MOD_PICInit,               ONLY: InitPIC
#ifdef MPI
USE MOD_part_MPI_Vars,         ONLY: PMPIVAR
USE MOD_Particle_MPI_Vars,     ONLY: SafetyFactor,halo_eps_velo
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iSpec, iInit, iPartBound, iSeed
  INTEGER               :: SeedSize, iPBC, iBC, iSwaps
  INTEGER               :: ALLOCSTAT
  CHARACTER(32)         :: hilf , hilf2
  CHARACTER(200)        :: tmpString
  LOGICAL               :: TrueRandom                                                  !
  INTEGER,ALLOCATABLE   :: iSeeds(:)
  REAL                  :: iRan, aVec, bVec   ! random numbers for random vectors
  REAL                  :: lineVector(3), v_drift_line, A_ins
  INTEGER               :: iVec, MaxNbrOfSpeciesSwaps
#ifdef MPI
#endif
!===================================================================================================================================
!#ifdef MPI
!   PMPIVAR%COMM   = MPI_COMM_WORLD
!   PMPIVAR%iProc  = myRank
!   PMPIVAR%nProcs = nProcessors
!   CALL MPI_COMM_GROUP(PMPIVAR%COMM,PMPIVAR%GROUP,IERROR)
!   PMPIVAR%GROUPWORLD=PMPIVAR%GROUP
!!IPWRITE(*,*)'INIT: GROUPWORLD',PMPIVAR%GROUPWORLD
!#endif

! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* RK3 and RK4 only */
print*, "SFSDRWE#"
ALLOCATE(Pt_temp(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT)  
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
#endif 


! predefine random vectors

NumRanVec = GETINT('Particles-NumberOfRandomVectors','100000')
IF ((usevMPF).OR.(useDSMC)) THEN
  ALLOCATE(RandomVec(NumRanVec, 3))
  RandomVec = 0
  DO iVec = 1, NumRanVec  ! calculation of NumRanVec different Vectors
    CALL RANDOM_NUMBER(iRan)
    bVec              = 1 - 2*iRan
    aVec              = SQRT(1 - bVec**2)
    RandomVec(iVec,1) = bVec
    CALL RANDOM_NUMBER(iRan)
    bVec              = Pi *2 * iRan
    RandomVec(iVec,2) = aVec * COS(bVec)
    RandomVec(iVec,3) = aVec * SIN(bVec)
  END DO
END IF

ALLOCATE(PartState(1:PDM%maxParticleNumber,1:6)       , &
         LastPartPos(1:PDM%maxParticleNumber,1:3)     , &
         Pt(1:PDM%maxParticleNumber,1:3)              , &
         PartSpecies(1:PDM%maxParticleNumber)         , &
         PDM%ParticleInside(1:PDM%maxParticleNumber)  , &
         PDM%nextFreePosition(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
PDM%ParticleInside(1:PDM%maxParticleNumber) = .FALSE.
LastPartPos(1:PDM%maxParticleNumber,1:3)    = 0.
IF (.NOT.DoRestart) THEN
  PartState          = 0.
!  LastPartPos        = 0.
  Pt                 = 0.
  PartSpecies        = 0
!  PDM%ParticleInside = .FALSE.                ! Initialize with no particles inside (will be filled in SetParticlePosition)
!  Pt_temp            = 0. ! gets initialized in RK4-stepper 
END IF

nSpecies = GETINT('Part-nSpecies','1')


! init varibale MPF per particle
IF (usevMPF) THEN
  enableParticleMerge = GETLOGICAL('Part-vMPFPartMerge','.FALSE.')
  IF (enableParticleMerge) THEN
    vMPFMergePolyOrder = GETINT('Part-vMPFMergePolOrder','2')
    vMPFMergeCellSplitOrder = GETINT('Part-vMPFCellSplitOrder','15')
    vMPFMergeParticleTarget = GETINT('Part-vMPFMergeParticleTarget','0')
    IF (vMPFMergeParticleTarget.EQ.0) WRITE(*,*) 'vMPFMergeParticleTarget equals zero: no merging is performed!'
    vMPFSplitParticleTarget = GETINT('Part-vMPFSplitParticleTarget','0')
    IF (vMPFSplitParticleTarget.EQ.0) WRITE(*,*) 'vMPFSplitParticleTarget equals zero: no split is performed!'
    vMPFMergeParticleIter = GETINT('Part-vMPFMergeParticleIter','100')
    vMPF_velocityDistribution = GETSTR('Part-vMPFvelocityDistribution','OVDR')
    ALLOCATE(vMPF_SpecNumElem(1:nElems,1:nSpecies))
  END IF
  ALLOCATE(PartMPF(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
  END IF
END IF

!WriteOutputMesh (=vtk mesh at start, seperate mesh for each proc
OutputMesh = GETLOGICAL('Part-WriteOutputMesh','.FALSE.')
           
! output of macroscopic values
WriteMacroValues = GETLOGICAL('Part-WriteMacroValues','.FALSE.')
MacroValSamplIterNum = GETINT('Part-IterationForMacroVal','1')
!ParticlePushMethod = GETSTR('Part-ParticlePushMethod','boris_leap_frog_scheme')
WriteFieldsToVTK = GETLOGICAL('Part-WriteFieldsToVTK','.FALSE.')

!!!! Logicals for Constant Pressure in Cells
! are particles to be ADDED to cells in order to reach constant pressure? Default YES
PartPressAddParts = GETLOGICAL('Part-ConstPressAddParts','.TRUE.')
! are particles to be REMOVED from cells in order to reach constant pressure? Default NO
PartPressRemParts = GETLOGICAL('Part-ConstPressRemParts','.FALSE.')

! Read particle species data
!nSpecies = CNTSTR('Part-Species-SpaceIC')

IF (nSpecies.LE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF
PartPressureCell = .FALSE.
ALLOCATE(Species(1:nSpecies))

DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(hilf)//'-nInits','0')
  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits)) 
    DO iInit = 0, Species(iSpec)%NumberOfInits
      ! set help characters
    IF(iInit.EQ.0)THEN
      hilf2=TRIM(hilf)
    ELSE ! iInit >0
      WRITE(UNIT=hilf2,FMT='(I2)') iInit
      hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
    END IF ! iInit
    ! get species values // only once
    IF(iInit.EQ.0)THEN
      !General Species Values
      Species(iSpec)%ChargeIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-ChargeIC','0.')
      Species(iSpec)%MassIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-MassIC','0.')
      Species(iSpec)%MacroParticleFactor   = GETREAL('Part-Species'//TRIM(hilf2)//'-MacroParticleFactor','1.')
    END IF ! iInit
    ! get emmission and init data
    Species(iSpec)%Init(iInit)%UseForInit           = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForInit','.TRUE.')
    Species(iSpec)%Init(iInit)%UseForEmission       = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForEmission','.TRUE.')
    Species(iSpec)%Init(iInit)%SpaceIC               = GETSTR('Part-Species'//TRIM(hilf2)//'-SpaceIC','cuboid')
    Species(iSpec)%Init(iInit)%velocityDistribution  = GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution','constant')
    Species(iSpec)%Init(iInit)%initialParticleNumber = GETINT('Part-Species'//TRIM(hilf2)//'-initialParticleNumber','0')
    Species(iSpec)%Init(iInit)%RadiusIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusIC','1.')
    Species(iSpec)%Init(iInit)%Radius2IC             = GETREAL('Part-Species'//TRIM(hilf2)//'-Radius2IC','0.')
    Species(iSpec)%Init(iInit)%RadiusICGyro          = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusICGyro','1.')
    Species(iSpec)%Init(iInit)%NormalIC              = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-NormalIC',3,'0. , 0. , 1.')
    Species(iSpec)%Init(iInit)%BasePointIC           = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BasePointIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector1IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector1IC',3,'1. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector2IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector2IC',3,'0. , 1. , 0.')
    Species(iSpec)%Init(iInit)%CuboidHeightIC        = GETREAL('Part-Species'//TRIM(hilf2)//'-CuboidHeightIC','1.')
    Species(iSpec)%Init(iInit)%CylinderHeightIC      = GETREAL('Part-Species'//TRIM(hilf2)//'-CylinderHeightIC','1.')
    Species(iSpec)%Init(iInit)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC','0.')
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%Amplitude             = GETREAL('Part-Species'//TRIM(hilf2)//'-Amplitude','0.01')
    Species(iSpec)%Init(iInit)%WaveNumber            = GETREAL('Part-Species'//TRIM(hilf2)//'-WaveNumber','2.')
    Species(iSpec)%Init(iInit)%maxParticleNumberX    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-x','0')
    Species(iSpec)%Init(iInit)%maxParticleNumberY    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-y','0')
    Species(iSpec)%Init(iInit)%maxParticleNumberZ    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-z','0')
    Species(iSpec)%Init(iInit)%Alpha                 = GETREAL('Part-Species'//TRIM(hilf2)//'-Alpha','0.')
    Species(iSpec)%Init(iInit)%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC','0.')
    Species(iSpec)%Init(iInit)%ConstantPressure      = GETREAL('Part-Species'//TRIM(hilf2)//'-ConstantPressure','0.')
    Species(iSpec)%Init(iInit)%ConstPressureRelaxFac = GETREAL('Part-Species'//TRIM(hilf2)//'-ConstPressureRelaxFac','1.')
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity','0.')
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(hilf2)//'-ParticleEmissionType','2')
    Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL('Part-Species'//TRIM(hilf2)//'-ParticleEmission','0.')
    Species(iSpec)%Init(iInit)%NSigma                = GETREAL('Part-Species'//TRIM(hilf2)//'-NSigma','10.')
    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'weibel') THEN
      Species(iSpec)%Init(iInit)%WeibelVeloPar       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPar','0')
      Species(iSpec)%Init(iInit)%WeibelVeloPer       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPer','0')
    END IF
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ. 'OneD-twostreaminstabilty') THEN
      Species(iSpec)%Init(iInit)%OneDTwoStreamVelo   = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamVelo','0')
      Species(iSpec)%Init(iInit)%OneDTwoStreamTransRatio = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamTransRatio','0')
    END IF

    !-- various checks/calculations after read-in of Species(i)%Init(iInit)%-data
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
      IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1).AND.(Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2)) &
        CALL abort(__STAMP__&
        ,' Wrong emission-type for virtual Pre-Inserting region!')
      IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell_lpn') &
        CALL abort(__STAMP__&
        ,' Only maxwell_lpn is implemened as velocity-distribution for virtual Pre-Inserting region!')
      IF (Species(iSpec)%Init(iInit)%initialParticleNumber.GT.0) THEN
        CALL abort(__STAMP__&
          ,' initialParticleNumber does not work for virtual Pre-Inserting. Use additional Init!')
      ELSE
        Species(iSpec)%Init(iInit)%UseForInit = .FALSE. !Do not make particle_emission stuff during initialization
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled, as both VPI and"
        SWRITE(*,*) "initialParticleNumber=0 detected for Species, Init ", iSpec, iInit
      END IF
      Species(iSpec)%Init(iInit)%VirtPreInsert = .TRUE.
      SWRITE(*,*) "Virtual Pre-Inserting is used for Species, Init ", iSpec, iInit
      IF (Species(iSpec)%Init(iInit)%PartDensity .EQ. 0.) THEN
        SWRITE(*,*) "WARNING: If VPI-BC is open, a backflow might not be compensated"
        SWRITE(*,*) "         (use PartDensity instead of ParticleEmission)!"
      END IF
    ELSE
      Species(iSpec)%Init(iInit)%VirtPreInsert = .FALSE.
    END IF
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         ((Species(iSpec)%Init(iInit)%ParticleEmission-INT(Species(iSpec)%Init(iInit)%ParticleEmission)).NE.0)) THEN
       CALL abort(__STAMP__&
       ,' If ParticleEmissionType = 2 (parts per iteration), ParticleEmission has to be an integer number')
    END IF
    IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 4).OR. &
        (Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6)) PartPressureCell = .TRUE.
    !--- normalize VeloVecIC and NormalIC (and BaseVector 1 & 2 IC for cylinder) for Inits
    IF (.NOT. ALL(Species(iSpec)%Init(iInit)%VeloVecIC(:).eq.0.)) THEN
      Species(iSpec)%Init(iInit)%VeloVecIC = Species(iSpec)%Init(iInit)%VeloVecIC            / &
        SQRT(Species(iSpec)%Init(iInit)%VeloVecIC(1)*Species(iSpec)%Init(iInit)%VeloVecIC(1) + &
        Species(iSpec)%Init(iInit)%VeloVecIC(2)*Species(iSpec)%Init(iInit)%VeloVecIC(2)      + &
        Species(iSpec)%Init(iInit)%VeloVecIC(3)*Species(iSpec)%Init(iInit)%VeloVecIC(3))
    END IF
    Species(iSpec)%Init(iInit)%NormalIC = Species(iSpec)%Init(iInit)%NormalIC /                 &
      SQRT(Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1) + &
      Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2) + &
      Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')&
        .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
        Species(iSpec)%Init(iInit)%BaseVector1IC =&
                  Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector1IC /     &
        SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)*Species(iSpec)%Init(iInit)%BaseVector1IC(1) + &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2)*Species(iSpec)%Init(iInit)%BaseVector1IC(2) + &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3)*Species(iSpec)%Init(iInit)%BaseVector1IC(3))
        Species(iSpec)%Init(iInit)%BaseVector2IC =&
                   Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector2IC /    &
        SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)*Species(iSpec)%Init(iInit)%BaseVector2IC(1) + &
        Species(iSpec)%Init(iInit)%BaseVector2IC(2)*Species(iSpec)%Init(iInit)%BaseVector2IC(2)      + &
        Species(iSpec)%Init(iInit)%BaseVector2IC(3)*Species(iSpec)%Init(iInit)%BaseVector2IC(3))
    END IF

    !--- stuff for calculating ParticleEmission from PartDensity
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.).AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'LD_insert')) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) &
        CALL abort(__STAMP__&
        , 'Only emission-type 1 is supported for PartDensity without LD!')
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
        IF  (((TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') ) THEN
          IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
            CALL abort(__STAMP__&
            ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
          ELSE
            !---calculation of Base-Area and corresponding component of VeloVecIC
            lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
              Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
            lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
              Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
            lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
              Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
            IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
               CALL abort(__STAMP__&
               ,'BaseVectors are parallel!')
            ELSE
              A_ins = SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + lineVector(3) * lineVector(3))
              lineVector = lineVector / A_ins
              v_drift_line = Species(iSpec)%Init(iInit)%VeloIC * &
                ( Species(iSpec)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(iSpec)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
                + Species(iSpec)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
              IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .EQ. 'cylinder' ) &
                A_ins = Pi * (Species(iSpec)%Init(iInit)%RadiusIC**2-Species(iSpec)%Init(iInit)%Radius2IC**2)
              !---calculation of particle flow (macroparticles/s) through boundary
              Species(iSpec)%Init(iInit)%ParticleEmission = Species(iSpec)%Init(iInit)%PartDensity &
                                                          / Species(iSpec)%MacroParticleFactor * v_drift_line * A_ins
              END IF
          END IF
        ELSE
          CALL abort(__STAMP__&
          ,'Only const. or maxwell(_lpn) is supported as velocityDistr. for PartDensity without LD!')
        END IF
      ELSE IF (Species(iSpec)%Init(iInit)%VirtPreInsert) THEN
        IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
               CALL abort(__STAMP__&
          ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
        ELSE
          SWRITE(*,*) "PartDensity is used for VPI of Species, Init ", iSpec, iInit
        END IF
      ELSE
        CALL abort(__STAMP__&
        ,'PartDensity without LD is only supported for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
      END IF
    END IF

    IF(iInit.EQ.0)THEN
      !!!for new case: check if to be included!!!
      IF((( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0)&
        .AND.(Species(iSpec)%Init(iInit)%ParticleEmission.EQ.0.) )  &
        .AND.(Species(iSpec)%Init(iInit)%PartDensity.EQ.0.) )       &
        .AND.(Species(iSpec)%Init(iInit)%ConstantPressure.EQ.0.)    &
        .AND.(Species(iSpec)%NumberOfInits.GT.0))       THEN 
        Species(iSpec)%StartnumberOfInits = 1 ! only new style paramaters defined (Part-Species(i)-Init(iInit)-***)
      ELSE
        Species(iSpec)%StartnumberOfInits = 0 ! old style parameters has been defined for inits/emissions (Part-Species(i)-***)
      END IF
      SWRITE(*,*) "StartnumberOfInits of Species ", iSpec, " = ", Species(iSpec)%StartnumberOfInits
    END IF ! iInit .EQ.0

  END DO ! iInit
END DO ! iSpec 

! Which Lorentz boost method should be used?
PartLorentzType = GETINT('Part-LorentzType','3')

! Read parameter for FastInitBackgroundMesh (FIBGM)
GEO%FIBGMdeltas(1:3)              = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
GEO%FactorFIBGM(1:3)              = GETREALARRAy('Part-FactorFIBGM',3,'1. , 1. , 1.')
GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

! Read in boundary parameters
nPartBound = GETINT('Part-nBounds','1.')
IF (nPartBound.LE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR: nPartBound .LE. 0:', nPartBound)
END IF
ALLOCATE(PartBound%SourceBoundName(1:nPartBound))
ALLOCATE(PartBound%TargetBoundCond(1:nPartBound))
ALLOCATE(PartBound%MomentumACC(1:nPartBound))
ALLOCATE(PartBound%WallTemp(1:nPartBound))
ALLOCATE(PartBound%TransACC(1:nPartBound))
ALLOCATE(PartBound%VibACC(1:nPartBound))
ALLOCATE(PartBound%RotACC(1:nPartBound))
ALLOCATE(PartBound%WallVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientCondition(1:nPartBound))
ALLOCATE(PartBound%AmbientTemp(1:nPartBound))
ALLOCATE(PartBound%AmbientMeanPartMass(1:nPartBound))
ALLOCATE(PartBound%AmbientBeta(1:nPartBound))
ALLOCATE(PartBound%AmbientVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientDens(1:nPartBound))
ALLOCATE(PartBound%AmbientDynamicVisc(1:nPartBound))
ALLOCATE(PartBound%AmbientThermalCond(1:nPartBound))

ALLOCATE(PartBound%Voltage(1:nPartBound))
ALLOCATE(PartBound%NbrOfSpeciesSwaps(1:nPartBound))
!--determine MaxNbrOfSpeciesSwaps for correct allocation
MaxNbrOfSpeciesSwaps=0
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I2)') iPartBound
  PartBound%NbrOfSpeciesSwaps(iPartBound)= GETINT('Part-Boundary'//TRIM(hilf)//'-NbrOfSpeciesSwaps','0')
  MaxNbrOfSpeciesSwaps=max(PartBound%NbrOfSpeciesSwaps(iPartBound),MaxNbrOfSpeciesSwaps)
END DO
IF (MaxNbrOfSpeciesSwaps.gt.0) THEN
  ALLOCATE(PartBound%ProbOfSpeciesSwaps(1:nPartBound))
  ALLOCATE(PartBound%SpeciesSwaps(1:2,1:MaxNbrOfSpeciesSwaps,1:nPartBound))
END IF
!--
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I2)') iPartBound
  tmpString = GETSTR('Part-Boundary'//TRIM(hilf)//'-Condition','open')
  SELECT CASE (TRIM(tmpString))
  CASE('open')
     PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
     PartBound%AmbientCondition(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-AmbientCondition','.FALSE.')
     IF(PartBound%AmbientCondition(iPartBound)) THEN
       PartBound%AmbientTemp(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientTemp','0')
       PartBound%AmbientMeanPartMass(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientMeanPartMass','0')
       PartBound%AmbientBeta(iPartBound) = &
       SQRT(PartBound%AmbientMeanPartMass(iPartBound)/(2*BoltzmannConst*PartBound%AmbientTemp(iPartBound)))
       PartBound%AmbientVelo(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-AmbientVelo',3,'0. , 0. , 0.')
       PartBound%AmbientDens(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientDens','0')
       PartBound%AmbientDynamicVisc(iPartBound)=&
           GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientDynamicVisc','1.72326582572253E-5') ! N2:T=288K
       PartBound%AmbientThermalCond(iPartBound)=&
           GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientThermalCond','2.42948500556027E-2') ! N2:T=288K
     END IF
     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
  CASE('reflective')
     PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
     PartBound%MomentumACC(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-MomentumACC','0')
     PartBound%WallTemp(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp','0')
     PartBound%TransACC(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-TransACC','0')
     PartBound%VibACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-VibACC','0')
     PartBound%RotACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotACC','0')
     PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-WallVelo',3,'0. , 0. , 0.')
     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
     IF (PartBound%NbrOfSpeciesSwaps(iPartBound).gt.0) THEN  
       !read Species to be changed at wall (in, out), out=0: delete
       PartBound%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-Boundary'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
       DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(iPartBound)
         WRITE(UNIT=hilf2,FMT='(I2)') iSwaps
         PartBound%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
             GETINTARRAY('Part-Boundary'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
       END DO
     END IF
  CASE('periodic')
     PartBound%TargetBoundCond(iPartBound) = PartBound%PeriodicBC
  CASE('simple_anode')
     PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleAnodeBC
  CASE('simple_cathode')
     PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleCathodeBC
  CASE DEFAULT
     SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
     CALL abort(__STAMP__&
         ,'Particle Boundary Condition does not exist')
  END SELECT
  PartBound%SourceBoundName(iPartBound) = GETSTR('Part-Boundary'//TRIM(hilf)//'-SourceName')
END DO
DEALLOCATE(PartBound%AmbientMeanPartMass)
DEALLOCATE(PartBound%AmbientTemp)
! Set mapping from field boundary to particle boundary index
ALLOCATE(PartBound%Map(1:nBCs))
PartBound%Map(:)=-10
DO iPBC=1,nPartBound
  DO iBC = 1, nBCs
    IF (BoundaryType(iBC,1).EQ.0) THEN
      PartBound%Map(iBC) = 1
      SWRITE(*,*)"PartBound",iPBC,"is internal bound, no mapping needed"
    END IF
    IF (TRIM(BoundaryName(iBC)).EQ.TRIM(PartBound%SourceBoundName(iPBC))) THEN
      PartBound%Map(iBC) = PartBound%TargetBoundCond(iPBC) 
      SWRITE(*,*)"Mapped PartBound",iPBC,"on FieldBound",BoundaryType(iBC,1),",i.e.:",TRIM(BoundaryName(iBC))
    END IF
  END DO
END DO
! Errorhandler for PartBound-Types that could not be mapped to the 
! FieldBound-Types.
DO iBC = 1,nBCs
  IF (PartBound%Map(iBC).EQ.-10) THEN
    CALL abort(__STAMP__&
    ,' PartBound%Map for Boundary is not set. iBC: :',iBC)
  END IF
END DO

ALLOCATE(PEM%Element(1:PDM%maxParticleNumber), PEM%lastElement(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(__STAMP__&
  ,' Cannot allocate PEM arrays!')
END IF
IF (useDSMC.OR.PartPressureCell) THEN
  ALLOCATE(PEM%pStart(1:nElems)                         , &
           PEM%pNumber(1:nElems)                        , &
           PEM%pEnd(1:nElems)                           , &
           PEM%pNext(1:PDM%maxParticleNumber)           , STAT=ALLOCSTAT) 
           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    , ' Cannot allocate DSMC PEM arrays!')
  END IF
END IF
IF (useDSMC) THEN
  ALLOCATE(PDM%PartInit(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    ,' Cannot allocate DSMC PEM arrays!')
  END IF
END IF

!--- Read Manual Time Step
useManualTimeStep = .FALSE.
ManualTimeStep = GETREAL('Particles-ManualTimeStep', '0.0')
IF (ManualTimeStep.GT.0.0) THEN
  useManualTimeStep=.True.
END IF

!--- initialize randomization (= random if one or more seeds are 0 or random is wanted)
nrSeeds = GETINT('Part-NumberOfRandomSeeds','0')
IF (nrSeeds.GT.0) THEN
   ALLOCATE(seeds(1:nrSeeds))
   DO iSeed = 1, nrSeeds
      WRITE(UNIT=hilf,FMT='(I2)') iSeed
      seeds(iSeed) = GETINT('Particles-RandomSeed'//TRIM(hilf),'0')
   END DO
END IF

CALL RANDOM_SEED(Size = SeedSize)                       ! Check for number of needed Seeds
TrueRandom = .FALSE.                             ! FALSE for defined random seed

IF (nrSeeds.GT.0) THEN
   IF (nrSeeds.NE.SeedSize) THEN
      IPWRITE(*,*) 'Error: Number of seeds for RNG must be ',SeedSize
      IPWRITE(*,*) 'Random RNG seeds are used'
      TrueRandom = .TRUE.
   END IF
   DO iSeed = 1, nrSeeds
      IF (Seeds(iSeed).EQ.0) THEN
         IPWRITE(*,*) 'Error: ',SeedSize,' seeds for RNG must be defined'
         IPWRITE(*,*) 'Random RNG seeds are used'
         TrueRandom = .TRUE.
      END IF
   END DO
ELSE
   TrueRandom = .TRUE.
END IF

IF (TrueRandom) THEN
   CALL RANDOM_SEED()
ELSE
#ifdef MPI
   Seeds(1:SeedSize) = Seeds(1:SeedSize)+PMPIVAR%iProc
#endif
   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
END IF

ALLOCATE(iseeds(SeedSize))
iseeds(:)=0
CALL RANDOM_SEED(GET = iseeds(1:SeedSize))
IPWRITE(*,*) 'Random seeds in PIC_init:'
DO iSeed = 1,SeedSize
   IPWRITE(*,*) iseeds(iSeed)
END DO
DEALLOCATE(iseeds)

DelayTime = GETREAL('Part-DelayTime','0.')
! init interpolation
CALL InitializeInterpolation() ! not any more required ! has to be called earliear
CALL InitPIC()
#ifdef MPI
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')
#endif /*MPI*/
!-- Finalizing InitializeVariables
CALL DomainUpdate()
IF(enableParticleMerge) THEN
 !IF (TRIM(InterpolationType).NE.'particle_position') CALL DefineElemT_inv()
 CALL DefinePolyVec(vMPFMergePolyOrder) 
 CALL DefineSplitVec(vMPFMergeCellSplitOrder)
END IF

END SUBROUTINE InitializeVariables

SUBROUTINE DomainUpdate()
!===================================================================================================================================
! Update of domain
!===================================================================================================================================
! MODULES
USE MOD_Globals!,               ONLY : UNIT_errOut,nProcessors
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,              ONLY : nElems, nNodes
USE MOD_PICDepo,                ONLY : InitializeDeposition
USE MOD_PreProc
USE MOD_part_boundary_periodic, ONLY : InitPeriodic
USE MOD_DSMC_Analyze,           ONLY : WriteOutputMesh
USE MOD_PICInit,                ONLY : InitPIC
USE MOD_ReadInTools
USE MOD_PICInterpolation,      ONLY: InitializeInterpolation
!USE MOD_part_pressure,          ONLY : ParticlePressureIni,ParticlePressureCellIni
USE MOD_Particle_Output_Vars,   ONLY : OutputMesh
#ifdef MPI
USE MOD_Equation_Vars,          ONLY : c
USE MOD_part_boundary_mpi,      ONLY : Initialize
USE MOD_CalcTimeStep,           ONLY : CalcTimeStep
USE MOD_PICDepo_Vars,           ONLY : DepositionType, r_sf
USE MOD_part_MPI_Vars,          ONLY : PMPIVAR, FIBGMCellPadding, SafetyFactor, halo_eps, NbrOfCases, casematrix
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,k,l,ElemID,iLocSide,iNode, iInit                                 !
INTEGER                          :: BGMimin,BGMimax,BGMkmin,BGMkmax,BGMlmin,BGMlmax             !
REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax                          !
INTEGER                          :: BGMCellXmax,BGMCellXmin                                     !
INTEGER                          :: BGMCellYmax,BGMCellYmin                                     !
INTEGER                          :: BGMCellZmax,BGMCellZmin                                     !
INTEGER                          :: ALLOCSTAT                                                   !
LOGICAL                          :: exitTrue
#ifdef MPI
INTEGER                          :: BGMCells, j, m, CurrentProc, Cell, Procs
INTEGER                          :: imin, imax, kmin, kmax, lmin, lmax                          !
INTEGER                          :: ii, kk, ll,j_offset
INTEGER                          :: nPaddingCellsX, nPaddingCellsY, nPaddingCellsZ
INTEGER                          :: nShapePaddingX, nShapePaddingY, nShapePaddingZ
INTEGER                          :: NbrOfBGMCells(0:PMPIVAR%nProcs-1)                           !
INTEGER                          :: Displacement(1:PMPIVAR%nProcs)                              !
INTEGER, ALLOCATABLE             :: BGMCellsArray(:),CellProcNum(:,:,:)                         !
INTEGER, ALLOCATABLE             :: GlobalBGMCellsArray(:), ReducedBGMArray(:)                  !
INTEGER                          :: ReducedNbrOfBGMCells(0:PMPIVAR%nProcs-1)                    !
INTEGER, ALLOCATABLE             :: CellProcList(:,:,:,:)
INTEGER                          :: tempproclist(0:PMPIVAR%nProcs-1)                            !
REAL                             :: deltaT, halo_eps_velo
INTEGER                          :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
INTEGER                          :: ind, Shift(1:3), iCase
! out commented
!INTEGER                          :: MPIGROUPMAP,MPIProcPIC                                      !
#endif
!===================================================================================================================================

! zeros
#ifdef MPI
ll=0
ii=0
kk=0
#endif

!#ifdef MPI
!   !--- If this MPI process does not contain particles, step out
!   IF (PMPIVAR%GROUP.EQ.MPI_GROUP_EMPTY) RETURN
!#endif
!--- calc min and max coordinates for mesh
xmin = 1.0E200
xmax = -1.0E200
ymin = 1.0E200
ymax = -1.0E200
zmin = 1.0E200
zmax = -1.0E200
DO iNode=1,nNodes
  xmin=MIN(xmin,GEO%NodeCoords(1,iNode))
  xmax=MAX(xmax,GEO%NodeCoords(1,iNode))
  ymin=MIN(ymin,GEO%NodeCoords(2,iNode))
  ymax=MAX(ymax,GEO%NodeCoords(2,iNode))
  zmin=MIN(zmin,GEO%NodeCoords(3,iNode))
  zmax=MAX(zmax,GEO%NodeCoords(3,iNode))
END DO
GEO%xmin = xmin
GEO%xmax = xmax
GEO%ymin = ymin
GEO%ymax = ymax
GEO%zmin = zmin
GEO%zmax = zmax

CALL InitPeriodic()
!CALL InitializeInterpolation()
CALL InitializeDeposition()
CALL InitPIC()



#ifdef MPI
  ! allocate and initialize MPINeighbor
  ALLOCATE(PMPIVAR%MPINeighbor(0:PMPIVAR%nProcs-1))
  PMPIVAR%MPINeighbor(:) = .FALSE.
#endif


   !--- calculate variables necessary for fast initial localization

#ifdef MPI
  CALL MPI_ALLREDUCE(GEO%xmin, GEO%xminglob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PMPIVAR%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%ymin, GEO%yminglob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PMPIVAR%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%zmin, GEO%zminglob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PMPIVAR%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%xmax, GEO%xmaxglob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PMPIVAR%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%ymax, GEO%ymaxglob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PMPIVAR%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%zmax, GEO%zmaxglob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PMPIVAR%COMM, IERROR)
#else
  GEO%xminglob=GEO%xmin
  GEO%yminglob=GEO%ymin
  GEO%zminglob=GEO%zmin
  GEO%xmaxglob=GEO%xmax
  GEO%ymaxglob=GEO%ymax
  GEO%zmaxglob=GEO%zmax
#endif   

IF (ASSOCIATED(GEO%FIBGM)) THEN
  DO i=GEO%FIBGMimin,GEO%FIBGMimax
    DO k=GEO%FIBGMkmin,GEO%FIBGMkmax
      DO l=GEO%FIBGMlmin,GEO%FIBGMlmax
        SDEALLOCATE(GEO%FIBGM(i,k,l)%Element)
#ifdef MPI
        SDEALLOCATE(GEO%FIBGM(i,k,l)%ShapeProcs)
        SDEALLOCATE(GEO%FIBGM(i,k,l)%PaddingProcs)
!           SDEALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs)
#endif
      END DO
    END DO
  END DO
  DEALLOCATE(GEO%FIBGM)
END IF
!--- compute number of background cells in each direction
BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)+1.00001)
BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)+0.99999)
BGMkmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)+1.00001)
BGMkmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)+0.99999)
BGMlmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)+1.00001)
BGMlmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)+0.99999)
!   IF(PMPIVAR%iProc.EQ.1) THEN
!     print*, "INT",INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
!     print*,"zmaxminglob", GEO%zminglob, GEO%zmax, GEO%zmin
!   END IF
!     STOP
   
#ifdef MPI
  !--- JN: For MPI communication, information also about the neighboring FIBGM cells is needed
  !--- AS: shouldn't we add up here the nPaddingCells? 
  !--- TS: What we need to do is increase the BGM area for shape_function ONLY
  !        Reason: if a particle moves outside the domain, there still needs to be a
  !                BGM with an associated ShapeProc at the particle position
  !        Particle may only move c*dt*Safetyfactor.

  IF (ManualTimeStep.EQ.0.0) THEN
    deltaT=CALCTIMESTEP()
  ELSE
    deltaT=ManualTimeStep
  END IF
  halo_eps_velo = GETREAL('Particles-HaloEpsVelo','0')
  IF (halo_eps_velo.EQ.0) halo_eps_velo = c
#if (PP_TimeDiscMethod==4 || PP_TimeDiscMethod==200 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==1000)
  IF (halo_eps_velo.EQ.c) THEN
     CALL abort(__STAMP__&
     , 'Halo Eps Velocity for MPI not defined')
  END IF
#endif
#if (PP_TimeDiscMethod==201)
  deltaT=CALCTIMESTEP()
  halo_eps = c*deltaT*SafetyFactor*3.8
#else
  halo_eps = halo_eps_velo*deltaT*SafetyFactor
#endif
  !SWRITE(*,*)'  - eps=',halo_eps
  IF (DepositionType.EQ.'shape_function') THEN
    BGMimax = INT((GEO%xmax+halo_eps-GEO%xminglob)/GEO%FIBGMdeltas(1)+1.00001)
    BGMimin = INT((GEO%xmin-halo_eps-GEO%xminglob)/GEO%FIBGMdeltas(1)+0.99999)
    BGMkmax = INT((GEO%ymax+halo_eps-GEO%yminglob)/GEO%FIBGMdeltas(2)+1.00001)
    BGMkmin = INT((GEO%ymin-halo_eps-GEO%yminglob)/GEO%FIBGMdeltas(2)+0.99999)
    BGMlmax = INT((GEO%zmax+halo_eps-GEO%zminglob)/GEO%FIBGMdeltas(3)+1.00001)
    BGMlmin = INT((GEO%zmin-halo_eps-GEO%zminglob)/GEO%FIBGMdeltas(3)+0.99999)
  END IF


!!   BGMimax = BGMimax + 1
!!   BGMimin = BGMimin - 1
!!   BGMkmax = BGMkmax + 1
!!   BGMkmin = BGMkmin - 1
!!   BGMlmax = BGMlmax + 1
!!   BGMlmin = BGMlmin - 1

#endif

!print*,"BGM-Indices:",PMPIVAR%iProc, BGMimin, BGMimax, BGMkmin, BGMkmax, BGMlmin, BGMlmax
GEO%FIBGMimax=BGMimax
GEO%FIBGMimin=BGMimin
GEO%FIBGMkmax=BGMkmax
GEO%FIBGMkmin=BGMkmin
GEO%FIBGMlmax=BGMlmax
GEO%FIBGMlmin=BGMlmin

ALLOCATE(GEO%FIBGM(BGMimin:BGMimax,BGMkmin:BGMkmax,BGMlmin:BGMlmax), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  WRITE(*,'(A,6(I0,A))')'Problem allocating GEO%FIBGM(',BGMimin,':',BGMimax,',', &
                                                        BGMkmin,':',BGMkmax,',', &
                                                        BGMlmin,':',BGMlmax,')'
 CALL abort(__STAMP__&
     , 'STOP' )
END IF

DO i = BGMimin,BGMimax
   DO k = BGMkmin,BGMkmax
      DO l = BGMlmin,BGMlmax
         GEO%FIBGM(i,k,l)%nElem = 0
      END DO
   END DO
END DO
!--- compute number of elements in each background cell
DO ElemID=1,nElems
   xmin = 1.0E200
   xmax = -1.0E200
   ymin = 1.0E200
   ymax = -1.0E200
   zmin = 1.0E200
   zmax = -1.0E200    
   DO iLocSide = 1,6
     DO iNode=1,4
       xmin=MIN(xmin,GEO%NodeCoords(1,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       xmax=MAX(xmax,GEO%NodeCoords(1,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       ymin=MIN(ymin,GEO%NodeCoords(2,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       ymax=MAX(ymax,GEO%NodeCoords(2,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       zmin=MIN(zmin,GEO%NodeCoords(3,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       zmax=MAX(zmax,GEO%NodeCoords(3,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
     END DO
   END DO
   !--- find minimum and maximum BGM cell for current element
   BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
   BGMCellXmax = MIN(BGMCellXmax,BGMimax)
   BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
   BGMCellXmin = MAX(BGMCellXmin,BGMimin)
   BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
   BGMCellYmax = MIN(BGMCellYmax,BGMkmax)
   BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
   BGMCellYmin = MAX(BGMCellYmin,BGMkmin)
   BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
   BGMCellZmax = MIN(BGMCellZmax,BGMlmax)
   BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
   BGMCellZmin = MAX(BGMCellZmin,BGMlmin)      
   DO i = BGMCellXmin,BGMCellXmax
      DO k = BGMCellYmin,BGMCellYmax
         DO l = BGMCellZmin,BGMCellZmax
            GEO%FIBGM(i,k,l)%nElem = GEO%FIBGM(i,k,l)%nElem + 1
         END DO
      END DO
   END DO
END DO
!--- allocate mapping variable and clean number for mapping (below)
DO i = BGMimin,BGMimax
   DO k = BGMkmin,BGMkmax
      DO l = BGMlmin,BGMlmax
         ALLOCATE(GEO%FIBGM(i,k,l)%Element(1:GEO%FIBGM(i,k,l)%nElem))
         GEO%FIBGM(i,k,l)%nElem = 0
      END DO
   END DO
END DO
!--- map elements to background cells
DO ElemID=1,nElems
   xmin = 1.0E200
   xmax = -1.0E200
   ymin = 1.0E200
   ymax = -1.0E200
   zmin = 1.0E200
   zmax = -1.0E200    
   DO iLocSide = 1,6
     DO iNode=1,4
       xmin=MIN(xmin,GEO%NodeCoords(1,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       xmax=MAX(xmax,GEO%NodeCoords(1,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       ymin=MIN(ymin,GEO%NodeCoords(2,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       ymax=MAX(ymax,GEO%NodeCoords(2,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       zmin=MIN(zmin,GEO%NodeCoords(3,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
       zmax=MAX(zmax,GEO%NodeCoords(3,GEO%ElemSideNodeID(iNode,iLocSide,ElemID)))
     END DO
   END DO
   BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
   BGMCellXmax = MIN(BGMCellXmax,BGMimax)
   BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
   BGMCellXmin = MAX(BGMCellXmin,BGMimin)
   BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
   BGMCellYmax = MIN(BGMCellYmax,BGMkmax)
   BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
   BGMCellYmin = MAX(BGMCellYmin,BGMkmin)
   BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
   BGMCellZmax = MIN(BGMCellZmax,BGMlmax)
   BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
   BGMCellZmin = MAX(BGMCellZmin,BGMlmin)     
   DO i = BGMCellXmin,BGMCellXmax
      DO k = BGMCellYmin,BGMCellYmax
         DO l = BGMCellZmin,BGMCellZmax
            GEO%FIBGM(i,k,l)%nElem = GEO%FIBGM(i,k,l)%nElem + 1    
            GEO%FIBGM(i,k,l)%Element(GEO%FIBGM(i,k,l)%nElem) = ElemID
         END DO
      END DO
   END DO
END DO

#ifdef MPI
!--- MPI stuff for background mesh (FastinitBGM)
BGMCells=0 
ALLOCATE(BGMCellsArray(1:(BGMimax-BGMimin+1)*(BGMkmax-BGMkmin+1)*(BGMlmax-BGMlmin+1)*3))
DO i=BGMimin, BGMimax  !Count BGMCells with Elements inside and save their indices in BGMCellsArray
  DO k=BGMkmin, BGMkmax
    DO l=BGMlmin, BGMlmax
      IF (GEO%FIBGM(i,k,l)%nElem .GT. 0) THEN
        BGMCellsArray(BGMCells*3+1)= i
        BGMCellsArray(BGMCells*3+2)= k
        BGMCellsArray(BGMCells*3+3)= l
        BGMCells=BGMCells+1
      END IF
    END DO !i
  END DO !k
END DO !l

!Communicate number of BGMCells
CALL MPI_ALLGATHER(BGMCells, 1, MPI_INTEGER, NbrOfBGMCells(0:PMPIVAR%nProcs-1), 1, MPI_INTEGER, PMPIVAR%COMM, IERROR) 
ALLOCATE(GlobalBGMCellsArray(1:SUM(NbrOfBGMCells)*3))
Displacement(1)=0
DO i=2, PMPIVAR%nProcs
  Displacement(i) = SUM(NbrOfBGMCells(0:i-2))*3
END DO
!Gather indices of every Procs' Cells
CALL MPI_ALLGATHERV(BGMCellsArray(1:BGMCells*3), BGMCells*3, MPI_INTEGER, GlobalBGMCellsArray, &    
                   & NbrOfBGMCells(0:PMPIVAR%nProcs-1)*3, Displacement, MPI_INTEGER, PMPIVAR%COMM, IERROR)

!--- JN: first: count required array size for ReducedBGMArray
!--- TS: Define padding stencil (max of halo and shape padding)
!        Reason: This padding is used to build the ReducedBGM, so any information 
!                outside this region is lost 
FIBGMCellPadding(1:3) = INT(halo_eps/GEO%FIBGMdeltas(1:3))+1
nShapePaddingX = 0
nShapePaddingY = 0
nShapePaddingZ = 0
IF (DepositionType.EQ.'shape_function') THEN
  nShapePaddingX = INT(r_sf/GEO%FIBGMdeltas(1)+0.9999999)
  nShapePaddingY = INT(r_sf/GEO%FIBGMdeltas(2)+0.9999999)
  nShapePaddingZ = INT(r_sf/GEO%FIBGMdeltas(3)+0.9999999)
! 0.999999 in order to prevent stencil to get too big in case of r_sf==c_int*deltas
!  -> worst case: last 0.000001 gets cut off -> insignificant
END IF
nPaddingCellsX = MAX(nShapePaddingX,FIBGMCellPadding(1))
nPaddingCellsY = MAX(nShapePaddingY,FIBGMCellPadding(2))
nPaddingCellsZ = MAX(nShapePaddingZ,FIBGMCellPadding(3))

j=0
CurrentProc=0
DO i=1, SUM(NbrOfBGMCells)*3, 3
  IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PMPIVAR%nProcs-1) THEN
    CurrentProc=CurrentProc+1
  END IF
  IF  (.NOT.(GlobalBGMCellsArray(i) .LT. BGMimin-nPaddingCellsX .OR. GlobalBGMCellsArray(i).GT. BGMimax+nPaddingCellsX &
      & .OR. GlobalBGMCellsArray(i+1) .LT. BGMkmin-nPaddingCellsY .OR. GlobalBGMCellsArray(i+1) .GT. BGMkmax+nPaddingCellsY &
      & .OR. GlobalBGMCellsArray(i+2) .LT. BGMlmin-nPaddingCellsZ .OR. GlobalBGMCellsArray(i+2) .GT. BGMlmax+nPaddingCellsZ &
      & .OR. CurrentProc .EQ. PMPIVAR%iProc)) THEN
    j=j+3
  END IF
END DO !i

! Periodic: ReducedBGMArray needs to include cells on the other side of periodic vectors
Vec1(1:3) = 0
Vec2(1:3) = 0
Vec3(1:3) = 0
IF (GEO%nPeriodicVectors.GT.0) THEN
  ! build case matrix
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    DO ind = 1,3
      Vec1(ind) = INT(GEO%PeriodicVectors(ind,1)/GEO%FIBGMdeltas(ind)+0.1)
    END DO
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    DO ind = 1,3
      Vec1(ind) = INT(GEO%PeriodicVectors(ind,1)/GEO%FIBGMdeltas(ind)+0.1)
      Vec2(ind) = INT(GEO%PeriodicVectors(ind,2)/GEO%FIBGMdeltas(ind)+0.1)
    END DO
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    DO ind = 1,3
      Vec1(ind) = INT(GEO%PeriodicVectors(ind,1)/GEO%FIBGMdeltas(ind)+0.1)
      Vec2(ind) = INT(GEO%PeriodicVectors(ind,2)/GEO%FIBGMdeltas(ind)+0.1)
      Vec3(ind) = INT(GEO%PeriodicVectors(ind,3)/GEO%FIBGMdeltas(ind)+0.1)
    END DO
  END IF
  CurrentProc=0
  DO i=1, SUM(NbrOfBGMCells)*3, 3
    DO iCase = 1, NbrOfCases
      IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
          (casematrix(iCase,2).EQ.0) .AND. &
          (casematrix(iCase,3).EQ.0)) CYCLE
      Shift(1:3) = casematrix(iCase,1)*Vec1(1:3) + &
                   casematrix(iCase,2)*Vec2(1:3) + &
                   casematrix(iCase,3)*Vec3(1:3)
      IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PMPIVAR%nProcs-1) THEN
        CurrentProc=CurrentProc+1
      END IF
      IF  (.NOT.(GlobalBGMCellsArray(i)  +Shift(1) .LT. BGMimin-nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i)  +Shift(1) .GT. BGMimax+nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i+1)+Shift(2) .LT. BGMkmin-nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+1)+Shift(2) .GT. BGMkmax+nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+2)+Shift(3) .LT. BGMlmin-nPaddingCellsZ &
           .OR.  GlobalBGMCellsArray(i+2)+Shift(3) .GT. BGMlmax+nPaddingCellsZ &
           .OR. CurrentProc .EQ. PMPIVAR%iProc)) THEN
        j=j+3
      END IF
    END DO !iCase
  END DO !i
END IF !nPeriodic>0

ALLOCATE(ReducedBGMArray(1:j))
!Reduce GlobalBGMCellsArray: erase cells far away from iprocs domain
!--- JN: ReducedBGMArray contains data only from other MPI procs!

IF (GEO%nPeriodicVectors.GT.0) THEN  !Periodic (can't be done below because ReducedBGMArray is sorted by proc)
  j=1
  CurrentProc=0
  ReducedBGMArray=0
  ReducedNbrOfBGMCells=0
  DO i=1, SUM(NbrOfBGMCells)*3, 3
    DO iCase = 1, NbrOfCases         ! This time INCLUDING non-moved
      Shift(1:3) = casematrix(iCase,1)*Vec1(1:3) + &
                   casematrix(iCase,2)*Vec2(1:3) + &
                   casematrix(iCase,3)*Vec3(1:3)
      IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PMPIVAR%nProcs-1) THEN
        CurrentProc=CurrentProc+1
      END IF
      IF  (.NOT.(GlobalBGMCellsArray(i)   +Shift(1) .LT. BGMimin-nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i)   +Shift(1) .GT. BGMimax+nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i+1) +Shift(2) .LT. BGMkmin-nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+1) +Shift(2) .GT. BGMkmax+nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+2) +Shift(3) .LT. BGMlmin-nPaddingCellsZ &
           .OR.  GlobalBGMCellsArray(i+2) +Shift(3) .GT. BGMlmax+nPaddingCellsZ &
           .OR.  CurrentProc .EQ. PMPIVAR%iProc)) THEN
        ReducedBGMArray(j)=GlobalBGMCellsArray(i)     +Shift(1)
        ReducedBGMArray(j+1)=GlobalBGMCellsArray(i+1) +Shift(2)
        ReducedBGMArray(j+2)=GlobalBGMCellsArray(i+2) +Shift(3)
        j=j+3
        ReducedNbrOfBGMCells(CurrentProc)=ReducedNbrOfBGMCells(CurrentProc)+1
      END IF
    END DO ! iCase
  END DO !i
ELSE ! non periodic case
  j=1
  CurrentProc=0
  ReducedBGMArray=0
  ReducedNbrOfBGMCells=0
  DO i=1, SUM(NbrOfBGMCells)*3, 3
    IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PMPIVAR%nProcs-1) THEN
      CurrentProc=CurrentProc+1
    END IF
    IF  (.NOT.(GlobalBGMCellsArray(i) .LT. BGMimin-nPaddingCellsX .OR. GlobalBGMCellsArray(i).GT. BGMimax+nPaddingCellsX &
         & .OR. GlobalBGMCellsArray(i+1) .LT. BGMkmin-nPaddingCellsY .OR. GlobalBGMCellsArray(i+1) .GT. BGMkmax+nPaddingCellsY &
         & .OR. GlobalBGMCellsArray(i+2) .LT. BGMlmin-nPaddingCellsZ .OR. GlobalBGMCellsArray(i+2) .GT. BGMlmax+nPaddingCellsZ &
         & .OR. CurrentProc .EQ. PMPIVAR%iProc)) THEN
      ReducedBGMArray(j)=GlobalBGMCellsArray(i)
      ReducedBGMArray(j+1)=GlobalBGMCellsArray(i+1)
      ReducedBGMArray(j+2)=GlobalBGMCellsArray(i+2)
      j=j+3
      ReducedNbrOfBGMCells(CurrentProc)=ReducedNbrOfBGMCells(CurrentProc)+1
    END IF
  END DO !i
END IF !periodic


!--- JN: Determine required size of CellProcList array (hope this works, everytime I try to again understand this
!        shape function parallelization stuff, I get confused...)
!--- JN: But therefore we first have to refill BGMCellsArray to not only contain
!        cells with PIC%FastInitBGM%nElem.GT.0 but also those adjacent to them!
!--- TS: Actually, not the adjacent cell needs to be considered but a shape_proc stencil
!        Usually, the shape function radius is chosen to be the size of one BGM, but this 
!        is not necessarily always true. Hence new shape_proc padding:
BGMCells=0 
DO i=BGMimin, BGMimax  !Count BGMCells with Elements inside or adjacent and save their indices in BGMCellsArray
  DO k=BGMkmin, BGMkmax
    DO l=BGMlmin, BGMlmax
      iMin=MAX(i-nShapePaddingX,BGMimin); iMax=MIN(i+nShapePaddingX,BGMimax)
      kMin=MAX(k-nShapePaddingY,BGMkmin); kMax=MIN(k+nShapePaddingY,BGMkmax)
      lMin=MAX(l-nShapePaddingZ,BGMlmin); lMax=MIN(l+nShapePaddingZ,BGMlmax)
      IF (SUM(GEO%FIBGM(iMin:iMax,kMin:kMax,lMin:lMax)%nElem) .GT. 0) THEN
        BGMCellsArray(BGMCells*3+1)= i
        BGMCellsArray(BGMCells*3+2)= k
        BGMCellsArray(BGMCells*3+3)= l
        BGMCells=BGMCells+1
      END IF
    END DO !i
  END DO !k
END DO !l

! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved 
! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
ALLOCATE(CellProcNum(BGMimin-nShapePaddingX:BGMimax+nShapePaddingX, &
                     BGMkmin-nShapePaddingY:BGMkmax+nShapePaddingY, &
                     BGMlmin-nShapePaddingZ:BGMlmax+nShapePaddingZ))
CellProcNum = 0
Procs = 0 ! = maximum number of procs in one BGM cell
DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
  IF((ReducedBGMArray(j).GE.BGMimin-nShapePaddingX).AND.(ReducedBGMArray(j).LE.BGMimax+nShapePaddingX))THEN
    IF((ReducedBGMArray(j+1).GE.BGMkmin-nShapePaddingY).AND.(ReducedBGMArray(j+1).LE.BGMkmax+nShapePaddingY))THEN
      IF((ReducedBGMArray(j+2).GE.BGMlmin-nShapePaddingZ).AND.(ReducedBGMArray(j+2).LE.BGMlmax+nShapePaddingZ))THEN !inside
        CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
        Procs = MAX(Procs, CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)))
      END IF
    END IF
  END IF
END DO
! allocate the temporary array
ALLOCATE(CellProcList(BGMimin-nShapePaddingX:BGMimax+nShapePaddingX, &
                      BGMkmin-nShapePaddingY:BGMkmax+nShapePaddingY, &
                      BGMlmin-nShapePaddingZ:BGMlmax+nShapePaddingZ, &
                      1:Procs))
CellProcList = -1

! fill array with proc numbers

CellProcNum = 0
j_offset = 0
DO CurrentProc = 0,PMPIVAR%nProcs-1
  DO j = 1+j_offset, ReducedNbrOfBGMCells(CurrentProc)*3-2+j_offset,3
    IF((ReducedBGMArray(j).GE.BGMimin-nShapePaddingX).AND.(ReducedBGMArray(j).LE.BGMimax+nShapePaddingX))THEN
      IF((ReducedBGMArray(j+1).GE.BGMkmin-nShapePaddingY).AND.(ReducedBGMArray(j+1).LE.BGMkmax+nShapePaddingY))THEN
        IF((ReducedBGMArray(j+2).GE.BGMlmin-nShapePaddingZ).AND.(ReducedBGMArray(j+2).LE.BGMlmax+nShapePaddingZ))THEN
          CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
          CellProcList(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2), &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2))) = CurrentProc
        END IF
      END IF
    END IF
  END DO
  j_offset = j_offset + ReducedNbrOfBGMCells(CurrentProc)*3
END DO
! fill real array
DO Cell=0, BGMCells-1
  TempProcList=0
  DO i = BGMCellsArray(Cell*3+1)-nShapePaddingX, BGMCellsArray(Cell*3+1)+nShapePaddingX
    DO k = BGMCellsArray(Cell*3+2)-nShapePaddingY, BGMCellsArray(Cell*3+2)+nShapePaddingY
      DO l = BGMCellsArray(Cell*3+3)-nShapePaddingZ, BGMCellsArray(Cell*3+3)+nShapePaddingZ
        DO m = 1,CellProcNum(i,k,l)
          TempProcList(CellProcList(i,k,l,m))=1       ! every proc that is within the stencil gets a 1
        END DO ! m
        ll = l
      END DO !l
      kk = k
    END DO !k
    ii = i
  END DO !i
  Procs=SUM(TempProcList)
  IF (Procs.NE.0) THEN
    ALLOCATE(GEO%FIBGM(ii-nShapePaddingX,kk-nShapePaddingY,ll-nShapePaddingZ)%ShapeProcs(1:Procs+1))
    GEO%FIBGM(ii-nShapePaddingX,kk-nShapePaddingY,ll-nShapePaddingZ)%ShapeProcs(1) = Procs
    j=2
    DO m=0,PMPIVAR%nProcs-1
      IF (TempProcList(m) .EQ. 1) THEN
        PMPIVAR%MPINeighbor(m) = .true.
        GEO%FIBGM(ii-nShapePaddingX,kk-nShapePaddingY,ll-nShapePaddingZ)%ShapeProcs(j)=m
        j=j+1
      END IF
    END DO !m
  END IF
END DO !Cell

   !Compare own BGMCells and their Neighbors with ReducedBGMArray and save other Processes in BGM-Cells
   !--- JN: ReducedBGMArray contains data only from other MPI procs!
   !--- JN: BGMCellsArray contains in index triplets (i,k,l) all BGM cells containing elements from the local MPI proc
   !        plus the index triplets of BGM cells adjacent to cells containing elements from the local MPI proc

!   !--- JN: First identify only procs that share the exact same BGM cell as I (SharedProcs)
!   Procs = 0
!   CellProcList=-1
!   DO Cell=0, BGMCells-1
!     TempProcList=0
!     i = BGMCellsArray(Cell*3+1)
!     k = BGMCellsArray(Cell*3+2)
!     l = BGMCellsArray(Cell*3+3)
!     IF (GEO%FIBGM(i,k,l)%nElem.EQ.0) CYCLE
!     CurrentProc=0
!     m=2
!     DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
!       !--- JN: Slide CurrentProc to the MPI Proc that the currently checked BGMCell belongs to
!       DO WHILE (j .GT. SUM(ReducedNbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PMPIVAR%nProcs-1)
!         CurrentProc=CurrentProc+1
!       END DO
!       IF (i .EQ. ReducedBGMArray(j) .AND. k .EQ. ReducedBGMArray(j+1) .AND. l .EQ. ReducedBGMArray(j+2)) THEN
!         IF (m .GT. MaxShapeProcs) THEN
!           CALL abort(__STAMP__,&
!                                'ERROR in Boundary_PIC.f90: Cellproclist can contain only MaxShapeProcs=',MaxShapeProcs,999.)
!         END IF
!         CellProcList(i,k,l,m)=CurrentProc
!         m=m+1
!         TempProcList(CurrentProc)=1
!       END IF
!     END DO !j
!     Procs=SUM(TempProcList)
!     ALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs(1:Procs+1)) 
!     GEO%FIBGM(i,k,l)%SharedProcs(1) = Procs
!     j=2
!     DO m=0,PMPIVAR%nProcs-1
!       IF (TempProcList(m) .EQ. 1) THEN
!         GEO%FIBGM(i,k,l)%SharedProcs(j)=m
!         j=j+1
!       END IF
!     END DO !m
!   END DO !Cell


! ----------------------------------------------------------------!
!--- AS: Do it again for Paddingcells
DEALLOCATE(CellProcList)
DEALLOCATE(CellProcNum)
!--- JN: Determine required size of CellProcList array (hope this works, everytime I try to again understand this
!        shape function parallelization stuff, I get confused...)
!--- JN: But therefore we first have to refill BGMCellsArray to not only contain
!        cells with PIC%FastInitBGM%nElem.GT.0 but also those adjacent and the paddingcells to them!
BGMCells=0
DO i=BGMimin, BGMimax  !Count BGMCells with Elements inside or adjacent and save their indices in BGMCellsArray
  DO k=BGMkmin, BGMkmax
    DO l=BGMlmin, BGMlmax
      iMin=MAX(i-nPaddingCellsX,BGMimin); iMax=MIN(i+nPaddingCellsX,BGMimax)
      kMin=MAX(k-nPaddingCellsY,BGMkmin); kMax=MIN(k+nPaddingCellsY,BGMkmax)
      lMin=MAX(l-nPaddingCellsZ,BGMlmin); lMax=MIN(l+nPaddingCellsZ,BGMlmax)
      IF (SUM(GEO%FIBGM(iMin:iMax,kMin:kMax,lMin:lMax)%nElem) .GT. 0) THEN
        BGMCellsArray(BGMCells*3+1)= i
        BGMCellsArray(BGMCells*3+2)= k
        BGMCellsArray(BGMCells*3+3)= l
        BGMCells=BGMCells+1
      END IF
    END DO !i
  END DO !k
END DO !l

! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved 
! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
ALLOCATE(CellProcNum(BGMimin-nPaddingCellsX:BGMimax+nPaddingCellsX, &
                     BGMkmin-nPaddingCellsY:BGMkmax+nPaddingCellsY, &
                     BGMlmin-nPaddingCellsZ:BGMlmax+nPaddingCellsZ))
CellProcNum = 0
Procs = 0
DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
   IF((ReducedBGMArray(j).GE.BGMimin-nPaddingCellsX).AND.(ReducedBGMArray(j).LE.BGMimax+nPaddingCellsX))THEN
     IF((ReducedBGMArray(j+1).GE.BGMkmin-nPaddingCellsY).AND.(ReducedBGMArray(j+1).LE.BGMkmax+nPaddingCellsY))THEN
       IF((ReducedBGMArray(j+2).GE.BGMlmin-nPaddingCellsZ).AND.(ReducedBGMArray(j+2).LE.BGMlmax+nPaddingCellsZ))THEN
        CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
        Procs = MAX(Procs, CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)))
       END IF
     END IF
   END IF
END DO
! allocate the temporary array
ALLOCATE(CellProcList(BGMimin-nPaddingCellsX:BGMimax+nPaddingCellsX, &
                      BGMkmin-nPaddingCellsY:BGMkmax+nPaddingCellsY, &
                      BGMlmin-nPaddingCellsZ:BGMlmax+nPaddingCellsZ, &
                      1:Procs))
CellProcList = -1

! fill array with proc numbers

CellProcNum = 0
j_offset = 0
DO CurrentProc = 0,PMPIVAR%nProcs-1
  DO j = 1+j_offset, j_offset+ReducedNbrOfBGMCells(CurrentProc)*3-2,3
    CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
    CellProcList(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2), &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2))) = CurrentProc
  END DO
  j_offset = j_offset + ReducedNbrOfBGMCells(CurrentProc)*3
END DO

! fill real array
DO Cell=0, BGMCells-1
  TempProcList=0
  DO i = BGMCellsArray(Cell*3+1)-nPaddingCellsX, BGMCellsArray(Cell*3+1)+nPaddingCellsX
    DO k = BGMCellsArray(Cell*3+2)-nPaddingCellsY, BGMCellsArray(Cell*3+2)+nPaddingCellsY
      DO l = BGMCellsArray(Cell*3+3)-nPaddingCellsZ, BGMCellsArray(Cell*3+3)+nPaddingCellsZ
        DO m = 1,CellProcNum(i,k,l)
          TempProcList(CellProcList(i,k,l,m))=1       ! every proc that is within the stencil gets a 1
        END DO ! m
        ll = l
      END DO !l
      kk = k
    END DO !k
    ii = i
  END DO !i
  Procs=SUM(TempProcList)
  IF (Procs.NE.0) THEN
    ALLOCATE(GEO%FIBGM(ii-nPaddingCellsX,kk-nPaddingCellsY,ll-nPaddingCellsZ)%PaddingProcs(1:Procs+1))
    GEO%FIBGM(ii-nPaddingCellsX,kk-nPaddingCellsY,ll-nPaddingCellsZ)%PaddingProcs(1) = Procs
    j=2
    DO m=0,PMPIVAR%nProcs-1
      IF (TempProcList(m) .EQ. 1) THEN
        GEO%FIBGM(ii-nPaddingCellsX,kk-nPaddingCellsY,ll-nPaddingCellsZ)%PaddingProcs(j)=m
        j=j+1
      END IF
    END DO !m
  END IF
END DO !Cell
DEALLOCATE(ReducedBGMArray, BGMCellsArray, CellProcList, GlobalBGMCellsArray, CellProcNum)
CALL Initialize()  ! Initialize parallel environment for particle exchange between MPI domains
#endif

!exitTrue=.false.
!DO i = 1,nSpecies
!  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
!    IF((Species(i)%Init(iInit)%ParticleEmissionType .EQ. 3).OR.(Species(i)%Init(iInit)%ParticleEmissionType .EQ. 5)) THEN
!      CALL ParticlePressureIni()
!      exitTrue=.true.
!      EXIT
!    END IF
!  END DO
!  IF (exitTrue) EXIT
!END DO
!
!exitTrue=.false.
!DO i = 1,nSpecies
!  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
!    IF ((Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4).OR.(Species(i)%Init(iInit)%ParticleEmissionType .EQ. 6)) THEN
!      CALL ParticlePressureCellIni()
!      exitTrue=.true.
!      EXIT
!    END IF
!  END DO
!  IF (exitTrue) EXIT
!END DO

IF (OutputMesh) CALL WriteOutputMesh()

END SUBROUTINE DomainUpdate

END MODULE MOD_ParticleInit
