#include "boltzplatz.h"

MODULE MOD_Particle_Mesh
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

INTERFACE InitParticleMesh
  MODULE PROCEDURE InitParticleMesh
END INTERFACE

INTERFACE FinalizeParticleMesh
  MODULE PROCEDURE FinalizeParticleMesh
END INTERFACE

INTERFACE InitFIBGM
  MODULE PROCEDURE InitFIBGM
END INTERFACE

INTERFACE SingleParticleToExactElement
  MODULE PROCEDURE SingleParticleToExactElement
END INTERFACE

INTERFACE SingleParticleToExactElementNoMap
  MODULE PROCEDURE SingleParticleToExactElementNoMap
END INTERFACE

INTERFACE InitElemVolumes
  MODULE PROCEDURE InitElemVolumes
END INTERFACE

PUBLIC::InitElemVolumes
PUBLIC::InitParticleMesh,FinalizeParticleMesh, InitFIBGM, SingleParticleToExactElement, SingleParticleToExactElementNoMap
!===================================================================================================================================
!
CONTAINS

SUBROUTINE InitParticleMesh()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars
USE MOD_Mesh_Vars,              ONLY:Elems,nElems,nSides,SideToElem,ElemToSide,offsetElem
!USE MOD_Particle_Surfaces_Vars, ONLY:neighborElemID,neighborLocSideID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ALLOCSTAT
INTEGER           :: iElem, ilocSide,SideID,flip,iSide
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH ...'
IF(ParticleMeshInitIsDone) CALL abort(__STAMP__&
     , ' Particle-Mesh is already initialized.')
! allocate and duplicate partElemToside
nTotalSides=nSides
nTotalElems=nElems
ALLOCATE(PartElemToSide(1:2,1:6,1:nTotalSides)    &
        ,PartSideToElem(1:5,1:nTotalSides)        &
        ,PartNeighborElemID(1:6,1:nTotalElems)    &
        ,PartNeighborLocSideID(1:6,1:nTotalElems) &
        ,STAT=ALLOCSTAT                      )
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
 ,'  Cannot allocate particle mesh vars!')

!--- Initialize Periodic Side Info
ALLOCATE(SidePeriodicType(1:nSides)) 
SidePeriodicType=0
DO iElem=1,nElems
  DO iLocSide=1,6
    SideID = ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    IF ((Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%BCindex.GT.0)           .AND. &
        (ASSOCIATED(Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%connection))) THEN
      SidePeriodicType(SideID) = -1
    END IF
  END DO
END DO

! copy
PartElemToSide=ElemToSide
PartSideToelem=SideToElem

! get conection
DO iElem=1,PP_nElems
  DO ilocSide=1,6
    flip = ElemToSide(E2S_FLIP,ilocSide,iElem)
    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    IF(flip.EQ.0)THEN
      ! SideID of slave
      PartneighborlocSideID(ilocSide,iElem)=SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      PartneighborElemID   (ilocSide,iElem)=SideToElem(S2E_NB_ELEM_ID,SideID)
    ELSE
      ! SideID of master
      PartneighborlocSideID(ilocSide,iElem)=SideToElem(S2E_LOC_SIDE_ID,SideID)
      PartneighborElemID   (ilocSide,iElem)=SideToElem(S2E_ELEM_ID,SideID)
    END IF
  END DO ! ilocSide
END DO ! Elem
!PartNeighborElemID=NeighborElemID
!PartNeighborLocSideID=NeighborLocSideID

ParticleMeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE FinalizeParticleMesh()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(PartElemToSide)
SDEALLOCATE(PartSideToElem)
SDEALLOCATE(PartNeighborElemID)
SDEALLOCATE(PartNeighborLocSideID)
SDEALLOCATE(PartBCSideList)
ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh


SUBROUTINE SingleParticleToExactElement(iPart,doHalo)                                                         
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,          ONLY:PartState,PEM,PDM,PartPosRef
USE MOD_TimeDisc_Vars,          ONLY:dt
USE MOD_Equation_Vars,          ONLY:c_inv,c
USE MOD_Particle_Mesh_Vars,     ONLY:Geo
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol,OneMepsilon,epsilonOne,SuperSampledNodes,NPartCurved,ElemBaryNGeo,doRefMapping
USE MOD_Particle_Surfaces_Vars, ONLY:ClipHit
USE MOD_Mesh_Vars,              ONLY:ElemToSide,XCL_NGeo
USE MOD_Eval_xyz,               ONLY:eval_xyz_elemcheck
USE MOD_Utils,                  ONLY:BubbleSortID
USE MOD_PICDepo_Vars,           ONLY:DepositionType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
LOGICAL,INTENT(IN)                :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                           :: ilocSide,SideID,i
LOGICAL                           :: InElementCheck,ParticleFound                                
REAL                              :: xi(1:3),vBary(1:3)
REAL,ALLOCATABLE                  :: Distance(:)
INTEGER,ALLOCATABLE               :: ListDistance(:)
REAL,PARAMETER                    :: eps=1e-8 ! same value as in eval_xyz_elem
REAL,PARAMETER                    :: eps2=1e-3
REAL                              :: epsOne,OneMeps
!===================================================================================================================================

epsOne=1.0+eps
OneMeps=1.0-eps
ParticleFound = .FALSE.
IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
     (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
     (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
   PDM%ParticleInside(iPart) = .FALSE.
   RETURN
END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
CellX = MIN(GEO%FIBGMimax,CellX)                             
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MIN(GEO%FIBGMjmax,CellY) 
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MIN(GEO%FIBGMkmax,CellZ)
!   print*,'cell indices',CellX,CellY,CellZ
!   print*,'number of cells in bgm',GEO%FIBGM(CellX,CellY,CellZ)%nElem
!   read*

!--- check all cells associated with this beckground mesh cell
nBGMElems=GEO%FIBGM(CellX,CellY,CellZ)%nElem
ALLOCATE( Distance(1:nBGMElems) &
        , ListDistance(1:nBGMElems) )


! get closest element barycenter
Distance=1.
ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  Distance(iBGMElem)=(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
                    +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
                    +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) 
  Distance(iBGMElem)=SQRT(Distance(iBGMElem))
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

!print*,'earlier',Distance,ListDistance
CALL BubbleSortID(Distance,ListDistance,nBGMElems)
!print*,'after',Distance,ListDistance

! loop through sorted list and start by closest element  
DO iBGMElem=1,nBGMElems
  ElemID=ListDistance(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  CALL Eval_xyz_elemcheck(PartState(iPart,1:3),xi,ElemID)
  IF(ALL(ABS(Xi).LT.ClipHit)) THEN ! particle inside
    InElementCheck=.TRUE.
  ELSE IF(ANY(ABS(Xi).GT.ClipHit))THEN ! particle outside
  !  print*,'ici'
    InElementCheck=.FALSE.
  ELSE ! particle at face,edge or node, check most possible point
    ! alter particle position
    ! 1) compute vector to cell centre
    vBary=ElemBaryNGeo(1:3,ElemID)-PartState(iPart,1:3)
    ! 2) move particle pos along vector
    PartState(iPart,1:3) = PartState(iPart,1:3)+eps2*VBary(1:3)
    CALL Eval_xyz_elemcheck(PartState(iPart,1:3),xi,ElemID)
    !print*,xi
    IF(ALL(ABS(Xi).LT.ClipHit)) THEN ! particle inside
      InElementCheck=.TRUE.
    ELSE
!      IPWRITE(*,*) ' PartPos', PartState(iPart,1:3)
!      IPWRITE(*,*) ' xi',      XI(1:3)
      !SWRITE(*,*) ' Particle not located!'
      !SWRITE(*,*) ' PartPos', PartState(iPart,1:3)
      InElementCheck=.FALSE.
    END IF
  END IF
  IF (InElementCheck) THEN !  !     print*,Element
 ! read*
    PEM%Element(iPart) = ElemID
    IF(DoRefMapping) PartPosRef(1:3,iPart) = Xi
    ParticleFound = .TRUE.
    EXIT
  END IF
END DO ! iBGMElem

! particle not found
IF (.NOT.ParticleFound) THEN
  PDM%ParticleInside(iPart) = .FALSE.
END IF
! deallocate lists
DEALLOCATE( Distance,ListDistance)
!read*
END SUBROUTINE SingleParticleToExactElement


SUBROUTINE SingleParticleToExactElementNoMap(iPart,doHALO,debug) 
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,              ONLY:NGeo
USE MOD_Particle_Vars,          ONLY:PartState,PEM,PDM,LastPartPos
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide
USE MOD_Particle_Mesh_Vars,     ONLY:Geo
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol,OneMepsilon,epsilonOne,ElemBaryNGeo,BezierControlPoints3D,SideType
USE MOD_Utils,                  ONLY:BubbleSortID
USE MOD_Particle_Intersection,  ONLY:ComputePlanarInterSectionBezier,ComputeBilinearIntersectionSuperSampled2
USE MOD_Particle_Intersection,  ONLY:ComputeBezierIntersection
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
LOGICAL,INTENT(IN)                :: doHalo
LOGICAL,INTENT(IN),OPTIONAL       :: debug
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                           :: ilocSide,SideID,flip
LOGICAL                           :: ParticleFound,isHit
REAL                              :: vBary(1:3),lengthPartTrajectory,tmpPos(3),xNodes(1:3,1:4),tmpLastPartPos(1:3)
REAL,ALLOCATABLE                  :: Distance(:)
INTEGER,ALLOCATABLE               :: ListDistance(:)
REAL,PARAMETER                    :: eps=1e-8 ! same value as in eval_xyz_elem
REAL,PARAMETER                    :: eps2=1e-3
REAL                              :: epsOne,OneMeps,locAlpha(6),xi,eta
REAL                              :: PartTrajectory(1:3)
!===================================================================================================================================

epsOne=1.0+eps
OneMeps=1.0-eps
ParticleFound = .FALSE.
IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
     (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
     (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
   PDM%ParticleInside(iPart) = .FALSE.
   RETURN
END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
CellX = MIN(GEO%FIBGMimax,CellX)                             
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MIN(GEO%FIBGMjmax,CellY) 
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MIN(GEO%FIBGMkmax,CellZ)

!--- check all cells associated with this beckground mesh cell
nBGMElems=GEO%FIBGM(CellX,CellY,CellZ)%nElem
ALLOCATE( Distance(1:nBGMElems) &
        , ListDistance(1:nBGMElems) )

! get closest element barycenter
Distance=1.
ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  Distance(iBGMElem)=(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
                    +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
                    +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) 
  Distance(iBGMElem)=SQRT(Distance(iBGMElem))
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

CALL BubbleSortID(Distance,ListDistance,nBGMElems)

! loop through sorted list and start by closest element  
tmpPos=PartState(iPart,1:3)
tmpLastPartPos(1:3)=LastPartPos(iPart,1:3)
LastPartPos(iPart,1:3)=PartState(iPart,1:3)
DO iBGMElem=1,nBGMElems
  ElemID=ListDistance(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  PartState(iPart,1:3)=ElemBaryNGeo(:,ElemID)
  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                           +PartTrajectory(2)*PartTrajectory(2) &
                           +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  lengthPartTrajectory=lengthPartTrajectory+epsilontol
  locAlpha=1.
  DO ilocSide=1,6
    !SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
    flip  = PartElemToSide(E2S_FLIP,ilocSide,ElemID)
    SELECT CASE(SideType(SideID))
    CASE(PLANAR)
      CALL ComputePlanarIntersectionBezier(ishit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                              ,xi                 &
                                                                              ,eta             ,iPart,flip,SideID)
                                                                              !,eta             ,iPart,ilocSide,SideID,ElemID)
    CASE(BILINEAR)
      xNodes(1:3,1)=BezierControlPoints3D(1:3,0   ,0   ,SideID)
      xNodes(1:3,2)=BezierControlPoints3D(1:3,NGeo,0   ,SideID)
      xNodes(1:3,3)=BezierControlPoints3D(1:3,NGeo,NGeo,SideID)
      xNodes(1:3,4)=BezierControlPoints3D(1:3,0   ,NGeo,SideID)
      CALL ComputeBiLinearIntersectionSuperSampled2(ishit,xNodes &
                                                          ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                        ,xi                       &
                                                                                        ,eta                ,iPart,flip,SideID)
 


!      CALL ComputeBiLinearIntersectionSuperSampled2(ishit,[BezierControlPoints3D(1:3,0   ,0   ,SideID)  &
!                                                          ,BezierControlPoints3D(1:3,NGeo,0   ,SideID)  &
!                                                          ,BezierControlPoints3D(1:3,NGeo,NGeo,SideID)  &
!                                                          ,BezierControlPoints3D(1:3,0   ,NGeo,SideID)] &
!                                                          ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                                        ,xi                       &
                                                                                        !,eta                ,iPart,flip,SideID)
    CASE(CURVED)
      CALL ComputeBezierIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                              ,xi                 &
                                                                              ,eta                ,iPart,SideID)
    END SELECT
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      IF((ABS(xi).GT.1.0).OR.(ABS(eta).GT.1.0)) locAlpha(ilocSide)=-1.0
    END IF
  END DO ! ilocSide
  !IF((present(debug)).and.(ipart.eq.39)) print*,'blabla',locAlpha
  IF(ALMOSTEQUAL(MAXVAL(locAlpha(:)),-1.0))THEN
    ! no intersection found and particle is in final element
    PartState(iPart,1:3)=tmpPos
    LastPartPos(iPart,1:3)=tmpLastPartPos
    PEM%Element(iPart) = ElemID
    ParticleFound=.TRUE.
    EXIT
  END IF
END DO ! iBGMElem

! particle not found
IF (.NOT.ParticleFound) THEN
  PDM%ParticleInside(iPart) = .FALSE.
END IF
! deallocate lists
DEALLOCATE( Distance,ListDistance)
!read*
END SUBROUTINE SingleParticleToExactElementNoMap


SUBROUTINE InitFIBGM()
!===================================================================================================================================
! Build Fast-Init-Background-Mesh.
! The BGM is a cartesian mesh for easier locating of particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,                        ONLY:GetRealArray
USE MOD_Particle_Surfaces,                  ONLY:GetSideType,GetBCSideType,BuildElementBasis
USE MOD_Particle_Surfaces_Vars,             ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO,nTotalElems
USE MOD_Particle_Surfaces_Vars,             ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis,ElemRadiusNGeo!,DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI,                       ONLY:InitHALOMesh
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================


!! Read parameter for FastInitBackgroundMesh (FIBGM)
GEO%FIBGMdeltas(1:3)              = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
GEO%FactorFIBGM(1:3)              = GETREALARRAY('Part-FactorFIBGM',3,'1. , 1. , 1.')
GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

CALL GetFIBGM()
#ifdef MPI
!CALL Initialize()  ! Initialize parallel environment for particle exchange between MPI domains
CALL InitHaloMesh()
! HALO mesh and region build. Unfortunately, the local FIBGM has to be extended to include the HALO elements :(
! rebuild is a local operation
CALL AddHALOCellsToFIBGM()
#endif /*MPI*/

IF(DoRefMapping) THEN
  CALL GetBCSideType()
ELSE
  CALL GetSideType()
END IF
ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:nTotalElems) &
        ,slenXiEtaZetaBasis(1:6,1:nTotalElems) &
        ,ElemRadiusNGeo(1:nTotalElems)         &
        ,ElemBaryNGeo(1:3,1:nTotalElems)       )
CALL BuildElementBasis()

END SUBROUTINE InitFIBGM


SUBROUTINE GetFIBGM()
!===================================================================================================================================
! build local FIBGM mesh for process local FIBGM mesh including HALO region
! mode 1: build local BGM and interconnections with other processes
! mode 2: rebuild BGM including HALO region
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D
USE MOD_Mesh_Vars,                          ONLY:XCL_NGeo
USE MOD_Mesh_Vars,                          ONLY:nSides,NGeo
USE MOD_Partilce_Periodic_BC,               ONLY:InitPeriodicBC
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO,nTotalElems,nTotalSides
USE MOD_PICDepo,                            ONLY:InitializeDeposition
USE MOD_Particle_Surfaces_Vars,             ONLY:DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI,                       ONLY:InitHALOMesh
USE MOD_Equation_Vars,                      ONLY:c
USE MOD_Particle_Mesh_Vars,                 ONLY:FIBGMCellPadding,PartElemToSide,PartSideToElem
USE MOD_PICDepo_Vars,                       ONLY:DepositionType, r_sf
USE MOD_Particle_MPI_Vars,                  ONLY:PartMPI,SafetyFactor,halo_eps_velo,halo_eps,halo_eps2
USE MOD_Particle_MPI_Vars,                  ONLY:NbrOfCases,casematrix
USE MOD_Particle_Vars,                      ONLY:manualtimestep
USE MOD_CalcTimeStep,                       ONLY:CalcTimeStep
#endif /*MPI*/

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN)    :: mode
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                  :: localXmin,localXmax,localymin,localymax,localzmin,localzmax
INTEGER                          :: BGMimin,BGMimax,BGMjmin,BGMjmax,BGMkmin,BGMkmax
REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                          :: iSide,iBGM,jBGM,kBGM,SideID,iElem,ilocSide
INTEGER                          :: BGMCellXmax,BGMCellXmin
INTEGER                          :: BGMCellYmax,BGMCellYmin
INTEGER                          :: BGMCellZmax,BGMCellZmin
INTEGER                          :: ALLOCSTAT
INTEGER                          :: iSpec,iProc
#ifdef MPI
REAL                             :: deltaT
INTEGER                          :: ii,jj,kk,i,j,k
INTEGER                          :: BGMCells,  m, CurrentProc, Cell, Procs
INTEGER                          :: imin, imax, kmin, kmax, jmin, jmax
INTEGER                          :: nPaddingCellsX, nPaddingCellsY, nPaddingCellsZ
INTEGER                          :: nShapePaddingX, nShapePaddingY, nShapePaddingZ
INTEGER                          :: NbrOfBGMCells(0:PartMPI%nProcs-1)
INTEGER                          :: Displacement(1:PartMPI%nProcs)
INTEGER, ALLOCATABLE             :: BGMCellsArray(:),CellProcNum(:,:,:)
INTEGER, ALLOCATABLE             :: GlobalBGMCellsArray(:), ReducedBGMArray(:)
INTEGER                          :: ReducedNbrOfBGMCells(0:PartMPI%nProcs-1)
INTEGER, ALLOCATABLE             :: CellProcList(:,:,:,:)
INTEGER                          :: tempproclist(0:PartMPI%nProcs-1)
INTEGER                          :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
INTEGER                          :: ind, Shift(1:3), iCase
INTEGER                          :: j_offset
#endif /*MPI*/
!===================================================================================================================================

! zeros
#ifdef MPI
ii=0
jj=0
kk=0
#endif /*MPI*/

!#ifdef MPI
!   !--- If this MPI process does not contain particles, step out
!   IF (PMPIVAR%GROUP.EQ.MPI_GROUP_EMPTY) RETURN
!#endif
!--- calc min and max coordinates for mesh
xmin = HUGE(1.0)
xmax =-HUGE(1.0)
ymin = HUGE(1.0)
ymax =-HUGE(1.0)
zmin = HUGE(1.0)
zmax =-HUGE(1.0)

! serch for min,max of BezierControlPoints, e.g. the convec hull of the domain
! more accurate, XCL_NGeo
DO iElem=1,nTotalElems
  xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
  xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
  ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
  ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
  zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
  zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
END DO ! iElem

!  DO iSide=1,nTotalSides
!    xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,iSide)))
!    xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,iSide)))
!    ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,iSide)))
!    ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,iSide)))
!    zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,iSide)))
!    zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,iSide)))
!  END DO ! iSide

GEO%xmin=xmin
GEO%xmax=xmax
GEO%ymin=ymin
GEO%ymax=ymax
GEO%zmin=zmin
GEO%zmax=zmax

#ifdef MPI
  ! allocate and initialize MPINeighbor
  ALLOCATE(PartMPI%isMPINeighbor(0:PartMPI%nProcs-1))
  PartMPI%isMPINeighbor(:) = .FALSE.
  PartMPI%nMPINeighbors=0

! get global min, max
  CALL MPI_ALLREDUCE(GEO%xmin, GEO%xminglob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PartMPI%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%ymin, GEO%yminglob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PartMPI%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%zmin, GEO%zminglob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PartMPI%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%xmax, GEO%xmaxglob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PartMPI%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%ymax, GEO%ymaxglob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PartMPI%COMM, IERROR)
  CALL MPI_ALLREDUCE(GEO%zmax, GEO%zmaxglob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PartMPI%COMM, IERROR)
#else
  GEO%xminglob=GEO%xmin
  GEO%yminglob=GEO%ymin
  GEO%zminglob=GEO%zmin
  GEO%xmaxglob=GEO%xmax
  GEO%ymaxglob=GEO%ymax
  GEO%zmaxglob=GEO%zmax
#endif   

  CALL InitPeriodicBC()
  IF(DoRefMapping) CALL ReshapeBezierSides()
  !CALL InitializeInterpolation() ! not any more required ! has to be called earliear
  CALL InitializeDeposition()     ! has to remain here, because domain can have changed
  !CALL InitPIC()                 ! does not depend on domain

! deallocate stuff // required for dynamic load balance
IF (ALLOCATED(GEO%FIBGM)) THEN
  DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
    DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
      DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%ShapeProcs)
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)
!           SDEALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs)
      END DO
    END DO
  END DO
  DEALLOCATE(GEO%FIBGM)
END IF

!--- compute number of background cells in each direction
BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)+1.00001)
BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)+0.99999)
BGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)+1.00001)
BGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)+0.99999)
BGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)+1.00001)
BGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)+0.99999)

#ifdef MPI
  !--- JN: For MPI communication, information also about the neighboring FIBGM cells is needed
  !--- AS: shouldn't we add up here the nPaddingCells? 
  !--- TS: What we need to do is increase the BGM area for shape_function ONLY
  !        Reason: if a particle moves outside the domain, there still needs to be a
  !                BGM with an associated ShapeProc at the particle position
  !        Particle may only move c*dt*Safetyfactor.
  !--- PO: modified for curved and shape-function influence
  !        c*dt*SafetyFactor+r_cutoff
  IF (ManualTimeStep.EQ.0.0) THEN
    deltaT=CALCTIMESTEP()
  ELSE
    deltaT=ManualTimeStep
  END IF
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
  halo_eps = halo_eps_velo*deltaT*SafetyFactor ! for RK too large
#endif
  halo_eps2=halo_eps*halo_eps
  IF (DepositionType.EQ.'shape_function') THEN
    BGMimax = INT((GEO%xmax+halo_eps-GEO%xminglob)/GEO%FIBGMdeltas(1)+1.00001)
    BGMimin = INT((GEO%xmin-halo_eps-GEO%xminglob)/GEO%FIBGMdeltas(1)+0.99999)
    BGMjmax = INT((GEO%ymax+halo_eps-GEO%yminglob)/GEO%FIBGMdeltas(2)+1.00001)
    BGMjmin = INT((GEO%ymin-halo_eps-GEO%yminglob)/GEO%FIBGMdeltas(2)+0.99999)
    BGMkmax = INT((GEO%zmax+halo_eps-GEO%zminglob)/GEO%FIBGMdeltas(3)+1.00001)
    BGMkmin = INT((GEO%zmin-halo_eps-GEO%zminglob)/GEO%FIBGMdeltas(3)+0.99999)
  END IF
#endif

!print*,"BGM-Indices:",PartMPI%iProc, BGMimin, BGMimax, BGMjmin, BGMjmax, BGMkmin, BGMkmax
GEO%FIBGMimax=BGMimax
GEO%FIBGMimin=BGMimin
GEO%FIBGMjmax=BGMjmax
GEO%FIBGMjmin=BGMjmin
GEO%FIBGMkmax=BGMkmax
GEO%FIBGMkmin=BGMkmin

! allocate space for BGM
ALLOCATE(GEO%FIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  WRITE(*,'(A,6(I0,A))')'Problem allocating GEO%FIBGM(',BGMimin,':',BGMimax,',', &
                                                        BGMjmin,':',BGMjmax,',', &
                                                        BGMkmin,':',BGMkmax,')'
#ifdef MPI
  iProc=PartMPI%MyRank
#else
  iProc=0
#endif /*MPI*/
  CALL abort(__STAMP__&
      , 'Problem allocating GEO%FIBGM!' )
END IF

! null number of element per BGM cell
DO iBGM = BGMimin,BGMimax
   DO jBGM = BGMjmin,BGMjmax
      DO kBGM = BGMkmin,BGMkmax
         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
      END DO
   END DO
END DO

!--- compute number of elements in each background cell
DO iElem=1,nTotalElems
  xmin = HUGE(1.0)
  xmax =-HUGE(1.0)
  ymin = HUGE(1.0)
  ymax =-HUGE(1.0)
  zmin = HUGE(1.0)
  zmax =-HUGE(1.0)

  ! use XCL_NGeo of each element :)
  xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
  xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
  ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
  ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
  zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
  zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  !! get min,max of BezierControlPoints of Element
  !DO iLocSide = 1,6
  !  SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, iElem)
  !  xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,SideID)))
  !  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,SideID)))
  !  ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,SideID)))
  !  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,SideID)))
  !  zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,SideID)))
  !  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,SideID)))
  !END DO ! ilocSide
  !--- find minimum and maximum BGM cell for current element
  BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmax = MIN(BGMCellXmax,BGMimax)
  BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmin = MAX(BGMCellXmin,BGMimin)
  BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmax = MIN(BGMCellYmax,BGMjmax)
  BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmin = MAX(BGMCellYmin,BGMjmin)
  BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmax = MIN(BGMCellZmax,BGMkmax)
  BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmin = MAX(BGMCellZmin,BGMkmin)      
  ! add ecurrent element to number of BGM-elems
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

!--- allocate mapping variable and clean number for mapping (below)
DO iBGM = BGMimin,BGMimax
  DO jBGM = BGMjmin,BGMjmax
    DO kBGM = BGMkmin,BGMkmax
      ALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%FIBGM(iBGM,jBGM,kBGM)%nElem))
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!--- map elements to background cells
DO iElem=1,nTotalElems
  xmin = HUGE(1.0)
  xmax =-HUGE(1.0)
  ymin = HUGE(1.0)
  ymax =-HUGE(1.0)
  zmin = HUGE(1.0)
  zmax =-HUGE(1.0)

  ! use XCL_NGeo of each element :)
  xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
  xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
  ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
  ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
  zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
  zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  !! get min,max of BezierControlPoints of Element
  !DO iLocSide = 1,6
  !  SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, iElem)
  !  xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,SideID)))
  !  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,SideID)))
  !  ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,SideID)))
  !  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,SideID)))
  !  zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,SideID)))
  !  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,SideID)))
  !END DO ! ilocSide

  ! same as above
  BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmax = MIN(BGMCellXmax,BGMimax)
  BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmin = MAX(BGMCellXmin,BGMimin)
  BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmax = MIN(BGMCellYmax,BGMjmax)
  BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmin = MAX(BGMCellYmin,BGMjmin)
  BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmax = MIN(BGMCellZmax,BGMkmax)
  BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmin = MAX(BGMCellZmin,BGMkmin)     
  ! add current Element to BGM-Elem
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1    
        GEO%FIBGM(iBGM,jBGM,kBGM)%Element(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem


!IF(mode.EQ.2) RETURN
#ifdef MPI
!--- MPI stuff for background mesh (FastinitBGM)
BGMCells=0 
ALLOCATE(BGMCellsArray(1:(BGMimax-BGMimin+1)*(BGMjmax-BGMjmin+1)*(BGMkmax-BGMkmin+1)*3))
DO iBGM=BGMimin, BGMimax  !Count BGMCells with Elements inside and save their indices in BGMCellsArray
  DO jBGM=BGMjmin, BGMjmax
    DO kBGM=BGMkmin, BGMkmax
      IF (GEO%FIBGM(iBGM,jBGM,kBGM)%nElem .GT. 0) THEN
        !print*,"1",BGMCells*3+1,iBGM
        !print*,"2",BGMCells*3+2,jBGM
        !print*,"3",BGMCells*3+3,kBGM
        BGMCellsArray(BGMCells*3+1)= iBGM
        BGMCellsArray(BGMCells*3+2)= jBGM
        BGMCellsArray(BGMCells*3+3)= kBGM
        BGMCells=BGMCells+1
      END IF
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!Communicate number of BGMCells
CALL MPI_ALLGATHER(BGMCells, 1, MPI_INTEGER, NbrOfBGMCells(0:PartMPI%nProcs-1), 1, MPI_INTEGER, PartMPI%COMM, IERROR) 
ALLOCATE(GlobalBGMCellsArray(1:SUM(NbrOfBGMCells)*3))
Displacement(1)=0
DO i=2, PartMPI%nProcs
  Displacement(i) = SUM(NbrOfBGMCells(0:i-2))*3
END DO
!print*,'displacement',displacement
!Gather indices of every Procs' Cells
CALL MPI_ALLGATHERV(BGMCellsArray(1:BGMCells*3), BGMCells*3, MPI_INTEGER, GlobalBGMCellsArray, &    
                   & NbrOfBGMCells(0:PartMPI%nProcs-1)*3, Displacement, MPI_INTEGER, PartMPI%COMM, IERROR)

!stop
!--- JN: first: count required array size for ReducedBGMArray
!--- TS: Define padding stencil (max of halo and shape padding)
!        Reason: This padding is used to build the ReducedBGM, so any information 
!                outside this region is lost 
FIBGMCellPadding(1:3) = INT(halo_eps/GEO%FIBGMdeltas(1:3))+1
! halo region already included in BGM
!FIBGMCellPadding(1:3) = 0
nShapePaddingX = 0
nShapePaddingY = 0
nShapePaddingZ = 0
IF (DepositionType.EQ.'shape_function') THEN
  nShapePaddingX = INT(r_sf/GEO%FIBGMdeltas(1)+0.9999999)
  nShapePaddingY = INT(r_sf/GEO%FIBGMdeltas(2)+0.9999999)
  nShapePaddingZ = INT(r_sf/GEO%FIBGMdeltas(3)+0.9999999)
 ! IF(mode.EQ.2) THEN
 !   IF((nShapePaddingX.EQ.0)    &
 !     .OR.(nShapePaddingY.EQ.0) &
 !     .OR.(nShapePaddingZ.EQ.0))THEN 
 !       CALL abort(__STAMP__,&
 !         'Error in stencil calculation for FIBGM and shape function')
 !   END IF
 ! END IF
! 0.999999 in order to prevent stencil to get too big in case of r_sf==c_int*deltas
!  -> worst case: last 0.000001 gets cut off -> insignificant
END IF
nPaddingCellsX = MAX(nShapePaddingX,FIBGMCellPadding(1))
nPaddingCellsY = MAX(nShapePaddingY,FIBGMCellPadding(2))
nPaddingCellsZ = MAX(nShapePaddingZ,FIBGMCellPadding(3))

j=0
CurrentProc=0
DO i=1, SUM(NbrOfBGMCells)*3, 3
  IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
    CurrentProc=CurrentProc+1
  END IF
  IF  (.NOT.(GlobalBGMCellsArray(i) .LT. BGMimin-nPaddingCellsX .OR. GlobalBGMCellsArray(i).GT. BGMimax+nPaddingCellsX &
      & .OR. GlobalBGMCellsArray(i+1) .LT. BGMjmin-nPaddingCellsY .OR. GlobalBGMCellsArray(i+1) .GT. BGMjmax+nPaddingCellsY &
      & .OR. GlobalBGMCellsArray(i+2) .LT. BGMkmin-nPaddingCellsZ .OR. GlobalBGMCellsArray(i+2) .GT. BGMkmax+nPaddingCellsZ &
      & .OR. CurrentProc .EQ. PartMPI%Myrank)) THEN
    j=j+3
  END IF
END DO !i

! Periodic: ReducedBGMArray needs to include cells on the other side of periodic vectors
! --- PO: CAUTION: changes throuogh curved
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
      IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
        CurrentProc=CurrentProc+1
      END IF
      IF  (.NOT.(GlobalBGMCellsArray(i)  +Shift(1) .LT. BGMimin-nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i)  +Shift(1) .GT. BGMimax+nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i+1)+Shift(2) .LT. BGMjmin-nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+1)+Shift(2) .GT. BGMjmax+nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+2)+Shift(3) .LT. BGMkmin-nPaddingCellsZ &
           .OR.  GlobalBGMCellsArray(i+2)+Shift(3) .GT. BGMkmax+nPaddingCellsZ &
           .OR. CurrentProc .EQ. PartMPI%MyRank)) THEN
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
      IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
        CurrentProc=CurrentProc+1
      END IF
      IF  (.NOT.(GlobalBGMCellsArray(i)   +Shift(1) .LT. BGMimin-nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i)   +Shift(1) .GT. BGMimax+nPaddingCellsX &
           .OR.  GlobalBGMCellsArray(i+1) +Shift(2) .LT. BGMjmin-nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+1) +Shift(2) .GT. BGMjmax+nPaddingCellsY &
           .OR.  GlobalBGMCellsArray(i+2) +Shift(3) .LT. BGMkmin-nPaddingCellsZ &
           .OR.  GlobalBGMCellsArray(i+2) +Shift(3) .GT. BGMkmax+nPaddingCellsZ &
           .OR.  CurrentProc .EQ. PartMPI%MyRank)) THEN
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
    IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
      CurrentProc=CurrentProc+1
    END IF
    IF  (.NOT.(GlobalBGMCellsArray(i)   .LT. BGMimin-nPaddingCellsX .OR. GlobalBGMCellsArray(i).GT.    BGMimax+nPaddingCellsX &
        & .OR. GlobalBGMCellsArray(i+1) .LT. BGMjmin-nPaddingCellsY .OR. GlobalBGMCellsArray(i+1) .GT. BGMjmax+nPaddingCellsY &
        & .OR. GlobalBGMCellsArray(i+2) .LT. BGMkmin-nPaddingCellsZ .OR. GlobalBGMCellsArray(i+2) .GT. BGMkmax+nPaddingCellsZ &
         & .OR. CurrentProc .EQ. PartMPI%MyRank)) THEN
      ReducedBGMArray(j  )=GlobalBGMCellsArray(i  )
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
DO iBGM=BGMimin, BGMimax  !Count BGMCells with Elements inside or adjacent and save their indices in BGMCellsArray
  DO jBGM=BGMjmin, BGMjmax
    DO kBGM=BGMkmin, BGMkmax
      iMin=MAX(iBGM-nShapePaddingX,BGMimin); iMax=MIN(iBGM+nShapePaddingX,BGMimax)
      jMin=MAX(jBGM-nShapePaddingY,BGMjmin); jMax=MIN(jBGM+nShapePaddingY,BGMjmax)
      kMin=MAX(kBGM-nShapePaddingZ,BGMkmin); kMax=MIN(kBGM+nShapePaddingZ,BGMkmax)
      IF (SUM(GEO%FIBGM(iMin:iMax,jMin:jMax,kMin:kMax)%nElem) .GT. 0) THEN
        ! debug here changed i,j,k to ibgm,jbgm,kbgm
        BGMCellsArray(BGMCells*3+1)= iBGM
        BGMCellsArray(BGMCells*3+2)= jBGM
        BGMCellsArray(BGMCells*3+3)= kBGM
        BGMCells=BGMCells+1
      END IF
    END DO !iBGM
  END DO !jBGM
END DO !kBGM
!print*,'BGMCellsArray',BGMCellsArray

! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved 
! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
!print*,'BGMimin',BGMimin,BGMimax
!print*,'nShapePaddingz',nShapePaddingx
!
!print*,'BGMjmin',BGMjmin,BGMjmax
!print*,'nShapePaddingy',nShapePaddingy
!
!print*,'BGMkmin',BGMkmin,BGMkmax
!print*,'nShapePaddingZ',nShapePaddingz
ALLOCATE(CellProcNum(BGMimin-nShapePaddingX:BGMimax+nShapePaddingX, &
                     BGMjmin-nShapePaddingY:BGMjmax+nShapePaddingY, &
                     BGMkmin-nShapePaddingZ:BGMkmax+nShapePaddingZ))
CellProcNum = 0
Procs = 0 ! = maximum number of procs in one BGM cell
DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
  IF((ReducedBGMArray(j).GE.BGMimin-nShapePaddingX).AND.(ReducedBGMArray(j).LE.BGMimax+nShapePaddingX))THEN
    IF((ReducedBGMArray(j+1).GE.BGMjmin-nShapePaddingY).AND.(ReducedBGMArray(j+1).LE.BGMjmax+nShapePaddingY))THEN
      IF((ReducedBGMArray(j+2).GE.BGMkmin-nShapePaddingZ).AND.(ReducedBGMArray(j+2).LE.BGMkmax+nShapePaddingZ))THEN !inside
        CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
        Procs = MAX(Procs, CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)))
      END IF
    END IF
  END IF
END DO
! allocate the temporary array
ALLOCATE(CellProcList(BGMimin-nShapePaddingX:BGMimax+nShapePaddingX, &
                      BGMjmin-nShapePaddingY:BGMjmax+nShapePaddingY, &
                      BGMkmin-nShapePaddingZ:BGMkmax+nShapePaddingZ, &
                      1:Procs))
CellProcList = -1

! fill array with proc numbers

CellProcNum = 0
j_offset = 0
DO CurrentProc = 0,PartMPI%nProcs-1
  DO j = 1+j_offset, ReducedNbrOfBGMCells(CurrentProc)*3-2+j_offset,3
    IF((ReducedBGMArray(j).GE.BGMimin-nShapePaddingX).AND.(ReducedBGMArray(j).LE.BGMimax+nShapePaddingX))THEN
      IF((ReducedBGMArray(j+1).GE.BGMjmin-nShapePaddingY).AND.(ReducedBGMArray(j+1).LE.BGMjmax+nShapePaddingY))THEN
        IF((ReducedBGMArray(j+2).GE.BGMkmin-nShapePaddingZ).AND.(ReducedBGMArray(j+2).LE.BGMkmax+nShapePaddingZ))THEN
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
  DO iBGM = BGMCellsArray(Cell*3+1)-nShapePaddingX, BGMCellsArray(Cell*3+1)+nShapePaddingX
    DO jBGM = BGMCellsArray(Cell*3+2)-nShapePaddingY, BGMCellsArray(Cell*3+2)+nShapePaddingY
      DO kBGM = BGMCellsArray(Cell*3+3)-nShapePaddingZ, BGMCellsArray(Cell*3+3)+nShapePaddingZ
        DO m = 1,CellProcNum(iBGM,jBGM,kBGM)
          TempProcList(CellProcList(iBGM,jBGM,kBGM,m))=1       ! every proc that is within the stencil gets a 1
        END DO ! m
        kk = kBGM
      END DO !kBGM
      jj = jBGM
    END DO !jBGM
    ii = iBGM
  END DO !iBGM
  Procs=SUM(TempProcList)
  IF (Procs.NE.0) THEN
    ALLOCATE(GEO%FIBGM(ii-nShapePaddingX,jj-nShapePaddingY,kk-nShapePaddingZ)%ShapeProcs(1:Procs+1))
    GEO%FIBGM(ii-nShapePaddingX,jj-nShapePaddingY,kk-nShapePaddingZ)%ShapeProcs(1) = Procs
    j=2
    DO m=0,PartMPI%nProcs-1
      IF (TempProcList(m) .EQ. 1) THEN
        IF(.NOT.PartMPI%isMPINeighbor(m))THEN
          !IF(mode.EQ.2)THEN
          !  IPWRITE(*,*) ' Warning, something wrong with halo region'
          !  CALL abort(__STAMP__&
          !      , ' Something wrong with Halo region' )
          !END IF
          PartMPI%isMPINeighbor(m) = .true.
          PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
        END IF
        GEO%FIBGM(ii-nShapePaddingX,jj-nShapePaddingY,kk-nShapePaddingZ)%ShapeProcs(j)=m
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
DO iBGM=BGMimin, BGMimax  !Count BGMCells with Elements inside or adjacent and save their indices in BGMCellsArray
  DO jBGM=BGMjmin, BGMjmax
    DO kBGM=BGMkmin, BGMkmax
      iMin=MAX(iBGM-nPaddingCellsX,BGMimin); iMax=MIN(iBGM+nPaddingCellsX,BGMimax)
      jMin=MAX(jBGM-nPaddingCellsY,BGMjmin); jMax=MIN(jBGM+nPaddingCellsY,BGMjmax)
      kMin=MAX(kBGM-nPaddingCellsZ,BGMkmin); kMax=MIN(kBGM+nPaddingCellsZ,BGMkmax)
      IF (SUM(GEO%FIBGM(iMin:iMax,jMin:jMax,kMin:kMax)%nElem) .GT. 0) THEN
        BGMCellsArray(BGMCells*3+1)= iBGM
        BGMCellsArray(BGMCells*3+2)= jBGM
        BGMCellsArray(BGMCells*3+3)= kBGM
        BGMCells=BGMCells+1
      END IF
    END DO !iBGM
  END DO !jBGM
END DO !kBGM

! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved 
! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
ALLOCATE(CellProcNum(BGMimin-nPaddingCellsX:BGMimax+nPaddingCellsX, &
                     BGMjmin-nPaddingCellsY:BGMjmax+nPaddingCellsY, &
                     BGMkmin-nPaddingCellsZ:BGMkmax+nPaddingCellsZ))
CellProcNum = 0
Procs = 0
DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
   IF((ReducedBGMArray(j).GE.BGMimin-nPaddingCellsX).AND.(ReducedBGMArray(j).LE.BGMimax+nPaddingCellsX))THEN
     IF((ReducedBGMArray(j+1).GE.BGMjmin-nPaddingCellsY).AND.(ReducedBGMArray(j+1).LE.BGMjmax+nPaddingCellsY))THEN
       IF((ReducedBGMArray(j+2).GE.BGMkmin-nPaddingCellsZ).AND.(ReducedBGMArray(j+2).LE.BGMkmax+nPaddingCellsZ))THEN
        CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
        Procs = MAX(Procs, CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)))
       END IF
     END IF
   END IF
END DO
! allocate the temporary array
ALLOCATE(CellProcList(BGMimin-nPaddingCellsX:BGMimax+nPaddingCellsX, &
                      BGMjmin-nPaddingCellsY:BGMjmax+nPaddingCellsY, &
                      BGMkmin-nPaddingCellsZ:BGMkmax+nPaddingCellsZ, &
                      1:Procs))
CellProcList = -1

! fill array with proc numbers

CellProcNum = 0
j_offset = 0
DO CurrentProc = 0,PartMPI%nProcs-1
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
  DO iBGM = BGMCellsArray(Cell*3+1)-nPaddingCellsX, BGMCellsArray(Cell*3+1)+nPaddingCellsX
    DO jBGM = BGMCellsArray(Cell*3+2)-nPaddingCellsY, BGMCellsArray(Cell*3+2)+nPaddingCellsY
      DO kBGM = BGMCellsArray(Cell*3+3)-nPaddingCellsZ, BGMCellsArray(Cell*3+3)+nPaddingCellsZ
        DO m = 1,CellProcNum(iBGM,jBGM,kBGM)
          TempProcList(CellProcList(iBGM,jBGM,kBGM,m))=1       ! every proc that is within the stencil gets a 1
        END DO ! m
        kk = kBGM
      END DO !l
      jj = jBGM
    END DO !k
    ii = iBGM
  END DO !i
  Procs=SUM(TempProcList)
  IF (Procs.NE.0) THEN
    ALLOCATE(GEO%FIBGM(ii-nPaddingCellsX,jj-nPaddingCellsY,kk-nPaddingCellsZ)%PaddingProcs(1:Procs+1))
    GEO%FIBGM(ii-nPaddingCellsX,jj-nPaddingCellsY,kk-nPaddingCellsZ)%PaddingProcs(1) = Procs
    j=2
    DO m=0,PartMPI%nProcs-1
      IF (TempProcList(m) .EQ. 1) THEN
        GEO%FIBGM(ii-nPaddingCellsX,jj-nPaddingCellsY,kk-nPaddingCellsZ)%PaddingProcs(j)=m
        j=j+1
      END IF
    END DO !m
  END IF
END DO !Cell
DEALLOCATE(ReducedBGMArray, BGMCellsArray, CellProcList, GlobalBGMCellsArray, CellProcNum)
#endif /*MPI*/

END SUBROUTINE GetFIBGM


#ifdef MPI
SUBROUTINE AddHALOCellsToFIBGM()
!===================================================================================================================================
! remap all elements including halo-elements into FIBGM
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D
USE MOD_Mesh_Vars,                          ONLY:XCL_NGeo
USE MOD_Mesh_Vars,                          ONLY:nSides,NGeo
USE MOD_Partilce_Periodic_BC,               ONLY:InitPeriodicBC
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO,nTotalElems,nTotalSides
USE MOD_PICDepo,                            ONLY:InitializeDeposition
USE MOD_Particle_Surfaces_Vars,             ONLY:DoRefMapping
USE MOD_Particle_MPI,                       ONLY:InitHALOMesh
USE MOD_Equation_Vars,                      ONLY:c
USE MOD_Particle_Mesh_Vars,                 ONLY:FIBGMCellPadding,PartElemToSide,PartSideToElem
USE MOD_PICDepo_Vars,                       ONLY:DepositionType, r_sf
USE MOD_Particle_MPI_Vars,                  ONLY:PartMPI,SafetyFactor,halo_eps_velo,halo_eps,halo_eps2
USE MOD_Particle_MPI_Vars,                  ONLY:NbrOfCases,casematrix
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN)    :: mode
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: BGMimin,BGMimax,BGMjmin,BGMjmax,BGMkmin,BGMkmax
REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                          :: iSide,iBGM,jBGM,kBGM,SideID,iElem,ilocSide
INTEGER                          :: BGMCellXmax,BGMCellXmin
INTEGER                          :: BGMCellYmax,BGMCellYmin
INTEGER                          :: BGMCellZmax,BGMCellZmin
INTEGER                          :: ALLOCSTAT
INTEGER                          :: iSpec,iProc
INTEGER                          :: ii,jj,kk,i,j,k
INTEGER                          :: BGMCells,  m, CurrentProc, Cell, Procs
INTEGER                          :: imin, imax, kmin, kmax, jmin, jmax
INTEGER                          :: nPaddingCellsX, nPaddingCellsY, nPaddingCellsZ
INTEGER                          :: nShapePaddingX, nShapePaddingY, nShapePaddingZ
INTEGER                          :: NbrOfBGMCells(0:PartMPI%nProcs-1)
INTEGER                          :: Displacement(1:PartMPI%nProcs)
INTEGER, ALLOCATABLE             :: BGMCellsArray(:),CellProcNum(:,:,:)
INTEGER, ALLOCATABLE             :: GlobalBGMCellsArray(:), ReducedBGMArray(:)
INTEGER                          :: ReducedNbrOfBGMCells(0:PartMPI%nProcs-1)
INTEGER, ALLOCATABLE             :: CellProcList(:,:,:,:)
INTEGER                          :: tempproclist(0:PartMPI%nProcs-1)
INTEGER                          :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
INTEGER                          :: ind, Shift(1:3), iCase
INTEGER                          :: j_offset
LOGICAL, ALLOCATABLE             :: ElementFound(:)
!===================================================================================================================================

! simplify writting
BGMimax=GEO%FIBGMimax
BGMimin=GEO%FIBGMimin
BGMjmax=GEO%FIBGMjmax
BGMjmin=GEO%FIBGMjmin
BGMkmax=GEO%FIBGMkmax
BGMkmin=GEO%FIBGMkmin

! delete all elements form FIBGM && zero nElems per BGM-cell
DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
  DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
    DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
      SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO
  END DO
END DO

ALLOCATE( ElementFound(1:nTotalElems) )
ElementFound = .FALSE.

!--- compute number of elements in each background cell
DO iElem=1,nTotalElems
  xmin = HUGE(1.0)
  xmax =-HUGE(1.0)
  ymin = HUGE(1.0)
  ymax =-HUGE(1.0)
  zmin = HUGE(1.0)
  zmax =-HUGE(1.0)

  ! use XCL_NGeo of each element :)
  IF(DoRefMapping)THEN
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  ELSE
    ! get min,max of BezierControlPoints of Element
    DO iLocSide = 1,6
      SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, iElem)
      xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,SideID)))
      xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,SideID)))
      ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,SideID)))
      ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,SideID)))
      zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,SideID)))
      zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,SideID)))
    END DO ! ilocSide
  END IF
  !--- find minimum and maximum BGM cell for current element
  BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmax = MIN(BGMCellXmax,BGMimax)
  BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmin = MAX(BGMCellXmin,BGMimin)
  BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmax = MIN(BGMCellYmax,BGMjmax)
  BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmin = MAX(BGMCellYmin,BGMjmin)
  BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmax = MIN(BGMCellZmax,BGMkmax)
  BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmin = MAX(BGMCellZmin,BGMkmin)      
  ! add ecurrent element to number of BGM-elems
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
         ElementFound(iElem) = .TRUE.
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

!--- allocate mapping variable and clean number for mapping (below)
DO iBGM = BGMimin,BGMimax
  DO jBGM = BGMjmin,BGMjmax
    DO kBGM = BGMkmin,BGMkmax
      ALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%FIBGM(iBGM,jBGM,kBGM)%nElem))
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!--- map elements to background cells
DO iElem=1,nTotalElems
  xmin = HUGE(1.0)
  xmax =-HUGE(1.0)
  ymin = HUGE(1.0)
  ymax =-HUGE(1.0)
  zmin = HUGE(1.0)
  zmax =-HUGE(1.0)

  ! use XCL_NGeo of each element :)
  IF(DoRefMapping)THEN
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  ELSE
    ! get min,max of BezierControlPoints of Element
    DO iLocSide = 1,6
      SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, iElem)
      xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,SideID)))
      xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,SideID)))
      ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,SideID)))
      ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,SideID)))
      zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,SideID)))
      zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,SideID)))
    END DO ! ilocSide
  END IF ! DoRefMapping

  ! same as above
  BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmax = MIN(BGMCellXmax,BGMimax)
  BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  BGMCellXmin = MAX(BGMCellXmin,BGMimin)
  BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmax = MIN(BGMCellYmax,BGMjmax)
  BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  BGMCellYmin = MAX(BGMCellYmin,BGMjmin)
  BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmax = MIN(BGMCellZmax,BGMkmax)
  BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  BGMCellZmin = MAX(BGMCellZmin,BGMkmin)     
  ! add current Element to BGM-Elem
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1    
        GEO%FIBGM(iBGM,jBGM,kBGM)%Element(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

DO iElem=1,nTotalElems
  IF(.NOT.ElementFound(iElem))THEN
    IPWRITE(*,*) ' FIBGM , iElem'
    IF(DoRefMapping)THEN
     ! IF(PartMPI%MyRank.EQ.1)THEN
       IPWRITE(*,*) 'xmin',GEO%xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem))
       IPWRITE(*,*) 'xmax',GEO%xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem))
       IPWRITE(*,*) 'ymin',GEO%ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem))
       IPWRITE(*,*) 'ymax',GEO%ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem))
       IPWRITE(*,*) 'zmin',GEO%zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem))
       IPWRITE(*,*) 'zmax',GEO%zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem))
        xmin=MINVAL(XCL_NGeo(1,:,:,:,iElem))
        xmax=MAXVAL(XCL_NGeo(1,:,:,:,iElem))
        ymin=MINVAL(XCL_NGeo(2,:,:,:,iElem))
        ymax=MAXVAL(XCL_NGeo(2,:,:,:,iElem))
        zmin=MINVAL(XCL_NGeo(3,:,:,:,iElem))
        zmax=MAXVAL(XCL_NGeo(3,:,:,:,iElem))
       IPWRITE(*,*) ' BGM , iBGM'
       IPWRITE(*,*) 'xmin', BGMimin,CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
       IPWRITE(*,*) 'xmax', BGMimax,CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
       IPWRITE(*,*) 'ymin', BGMjmin,CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
       IPWRITE(*,*) 'ymax', BGMjmax,CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
       IPWRITE(*,*) 'zmin', BGMkmin,CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
       IPWRITE(*,*) 'zmax', BGMkmax,CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
     ! END IF
    END IF
    CALL abort(__STAMP__&
    ,' Element not located in FIBGM! iElem, myRank',iElem,REAL(PartMPI%MyRank))
  END IF
END DO ! iElem

END SUBROUTINE AddHALOCellsToFIBGM
#endif /*MPI*/


SUBROUTINE InitElemVolumes()
!===================================================================================================================================
! Calculate Element volumes for later use in particle routines
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_Mesh_Vars,          ONLY:nElems,NGeo,sJ
USE MOD_Particle_Mesh_Vars, ONLY:GEO
USE MOD_Interpolation_Vars, ONLY:wGP
USE MOD_Particle_Vars,      ONLY:usevMPF
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


SUBROUTINE ReShapeBezierSides()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalBCSides,PartBCSideList
USE MOD_Mesh_Vars,              ONLY:nSides,nBCSides,NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars, ONLY:SlabNormals,SlabIntervalls,BoundingBoxIsEmpty
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ALLOCSTAT
INTEGER           :: iSide,nPeriodicSides,nOldBCSides,newBCSideID,BCInc
REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummySlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummySlabIntervalls               ! intervalls beta1, beta2, beta3
LOGICAL,ALLOCATABLE,DIMENSION(:)   :: DummyBoundingBoxIsEmpty
REAL,ALLOCATABLE                   :: DummyBezierControlPoints3D(:,:,:,:)                                
!===================================================================================================================================


nPeriodicSides=0
DO iSide=nBCSides+1,nSides
  IF(SidePeriodicType(iSide).NE.0) nPeriodicSides=nPeriodicSides+1
END DO ! iSide

! now, shrink partbcsidelist
nOldBCSides  =nTotalBCSides
nTotalBCSides=nTotalBCSides+nPeriodicSides

! allocate & fill dummy
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
! SlabNormals
ALLOCATE(DummySlabNormals(1:3,1:3,1:nOldBCSides))
IF (.NOT.ALLOCATED(DummySlabNormals)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummySlabNormals=SlabNormals
DEALLOCATE(SlabNormals)
ALLOCATE(SlabNormals(1:3,1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
! SlabIntervalls
ALLOCATE(DummySlabIntervalls(1:6,1:nOldBCSides))
IF (.NOT.ALLOCATED(DummySlabIntervalls)) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
DummySlabIntervalls=SlabIntervalls
DEALLOCATE(SlabIntervalls)
ALLOCATE(SlabIntervalls(1:6,1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,& !wunderschoen!!!
  'Could not allocate ElemIndex')
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

BCInc=0
DO iSide=1,nSides
  IF((iSide.LE.nBCSides).OR.(SidePeriodicType(iSide).NE.0))THEN
    IF(iSide.LE.nBCSides)THEN
      newBCSideID=iSide
    ELSE IF(SidePeriodicType(iSide).NE.0)THEN
      BCInc=BCInc+1
      newBCSideID=nBCSides+BCInc
      PartBCSideList(iSide)=newBCSideID
    END IF
    BezierControlPoints3d(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
    SlabNormals          (1:3,1:3,          newBCSideID) =DummySlabNormals         (1:3,1:3,           iSide)
    SlabIntervalls       (1:6,              newBCSideID) =DummySlabIntervalls      (1:6,               iSide)
    BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
  END IF
END DO ! iSide

! deallocate dummy buffer
DEALLOCATE(DummyBezierControlPoints3D)
DEALLOCATE(DummySlabNormals)
DEALLOCATE(DummySlabIntervalls)
DEALLOCATE(DummyBoundingBoxIsEmpty)


END SUBROUTINE ReShapeBezierSides

END MODULE MOD_Particle_Mesh
