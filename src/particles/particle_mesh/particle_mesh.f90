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

MODULE MOD_Particle_Mesh
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
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

!INTERFACE InitFIBGM
!  MODULE PROCEDURE InitFIBGM
!END INTERFACE

!INTERFACE InitElemBoundingBox
!  MODULE PROCEDURE InitElemBoundingBox
!END INTERFACE

PUBLIC::DefineParametersParticleMesh
PUBLIC::InitParticleMesh
PUBLIC::FinalizeParticleMesh
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle tracking
!==================================================================================================================================
SUBROUTINE DefineParametersParticleMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('Tracking')

CALL prms%CreateIntFromStringOption('TrackingMethod', "Define Method that is used for tracking of particles:\n"//&
                                                      "refmapping (1): reference mapping of particle position"//&
                                                      " with (bi-)linear and bezier (curved) description of sides.\n"//&
                                                      "tracing (2): tracing of particle path "//&
                                                      "with (bi-)linear and bezier (curved) description of sides.\n"//&
                                                      "triatracking (3): tracing of particle path "//&
                                                      "with triangle-aproximation of (bi-)linear sides.\n", &
                                                      "triatracking")
CALL addStrListEntry('TrackingMethod' , 'refmapping'      , REFMAPPING)
CALL addStrListEntry('TrackingMethod' , 'tracing'         , TRACING)
CALL addStrListEntry('TrackingMethod' , 'triatracking'    , TRIATRACKING)
CALL addStrListEntry('TrackingMethod' , 'default'         , TRIATRACKING)

CALL prms%CreateLogicalOption( 'TriaSurfaceFlux'&
  , 'Using Triangle-aproximation [T] or (bi-)linear and bezier (curved) description [F] of sides for surfaceflux.'//&
  ' Default is set to TriaTracking')
CALL prms%CreateLogicalOption( 'DisplayLostParticles'&
  , 'Display position, velocity, species and host element of particles lost during particle tracking (TrackingMethod = '//&
    'triatracking, tracing)','.FALSE.')
CALL prms%CreateLogicalOption( 'CountNbrOfLostParts'&
    , 'Count the number of lost particles during tracking that cannot be found with fallbacks. Additionally, the lost particle '//&
    'information is stored in a PartStateLost*.h5 file. When particles are not found during restart in their host cell '//&
    '(sanity check), they are marked missing and are also written to PartStateLost*.h5 file even if they are re-located '//&
    'on a different processor.','.TRUE.')
CALL prms%CreateIntOption(     'PartOut'&
  , 'If compiled with CODE_ANALYZE flag: For This particle number every tracking information is written as STDOUT.','0')
CALL prms%CreateIntOption(     'MPIRankOut'&
  , 'If compiled with CODE_ANALYZE flag: This MPI-Proc writes the tracking information for the defined PartOut.','0')
CALL prms%CreateLogicalOption( 'MeasureTrackTime'&
  , 'If .TRUE. then the time how long the tracking routines are called are sampled and written for each MPI-Proc.','.FALSE.')
CALL prms%CreateLogicalOption( 'CartesianPeriodic'&
    , ' Simplified treatment for periodic box with Refmapping. Not computation of intersection points at periodic BCs.','.FALSE.')
CALL prms%CreateLogicalOption( 'FastPeriodic'&
  , ' Further simplification by directly moving particle into grid. Instead of moving the particle several times the periodic'//&
    ' displacements, the particle is mapped directly back into the domain. ','.FALSE.')
CALL prms%CreateIntOption(     'RefMappingGuess'&
  , ' Initial guess of the Newton for mapping the particle into reference coordinates.\n'//&
    '1 -linear pseudo-Cartesian coordinates\n'//&
    '2 - Xi of closest Gauss point\n'//&
    '3 - Xi of closest XCL_ngeo point\n'//&
    '4 -trival guess (0,0,0)^t')
CALL prms%CreateRealOption(    'RefMappingEps'&
  , ' Tolerance for mapping particle into reference element measured as L2-norm of deltaXi' , '1e-4')
CALL prms%CreateIntOption(     'BezierElevation'&
  , ' Use BezierElevation>0 to tighten the bounding box. Typical values>10','0')
CALL prms%CreateIntOption(     'BezierSampleN'&
  , 'TODO-DEFINE-PARAMETER\n'//&
    'Default value: NGeo equidistant sampling of bezier surface for emission','0')

CALL prms%CreateLogicalOption( 'CalcHaloInfo',         'Output halo element information to ElemData for each processor'//&
                                                       ' "MyRank_ElemHaloInfo"\n'//&
                                                       ' ElemHaloInfo = -1            : element not in list\n'//&
                                                       '              = 0             : halo elements\n'//&
                                                       '              = 1 to PP_nElems: local elements','.FALSE.')
CALL prms%CreateRealOption(    'BezierNewtonAngle'      , ' BoundingBox intersection angle for switching between '//&
'Bezierclipping and BezierNewton.' , '1.570796326')
CALL prms%CreateRealOption(    'BezierClipTolerance'    , ' Tolerance for BezierClipping' , '1e-8')
CALL prms%CreateRealOption(    'BezierNewtonTolerance'  , ' Tolerance for BezierNewton' , '1e-4')
CALL prms%CreateIntOption(     'BezierNewtonGuess'      , ' Initial guess for BezierNewton\n'// &
    '1 - linear projected face\n'//&
    '2 - closest projected BeziercontrolPoint\n'//&
    '4 - (0,0)^t' , '1')
CALL prms%CreateIntOption(     'BezierNewtonMaxIter'    , ' TODO-DEFINE-PARAMETER' , '100')
CALL prms%CreateRealOption(    'BezierSplitLimit'       , ' Limit for splitting in BezierClipping.'// &
   ' Value allows to detect multiple intersections and speed up computation. Parameter is multiplied by 2' , '0.6')
CALL prms%CreateIntOption(     'BezierClipMaxIter'      , ' Max iteration of BezierClipping' , '100')
CALL prms%CreateIntOption(     'BezierClipLineVectorMethod' , ' TODO-DEFINE-PARAMETER' , '2')
CALL prms%CreateRealOption(    'epsilontol'             , 'TODO-DEFINE-PARAMETER' , '0.')
CALL prms%CreateRealOption(    'BezierClipHit'          , ' Tolerance in [-1,1] of BezierFace' , '0.')
CALL prms%CreateRealOption(    'BezierNewtonHit'        , ' Tolerance in [-1,1] of BezierNewton' , '0.')
CALL prms%CreateIntOption(     'BezierClipMaxIntersec'  , ' Max. number of multiple intersections. Default: 2*(NGeo+1)')

END SUBROUTINE DefineParametersParticleMesh


SUBROUTINE InitParticleMesh()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Tools             ,ONLY: InitGetGlobalElemID,InitGetCNElemID,GetCNElemID
USE MOD_Mesh_Tools             ,ONLY: InitGetGlobalSideID,InitGetCNSideID,GetGlobalSideID
USE MOD_Mesh_Vars              ,ONLY: deleteMeshPointer,NodeCoords
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Mesh_Vars              ,ONLY: useCurveds
USE MOD_Analyze_Vars           ,ONLY: CalcHaloInfo
USE MOD_Particle_BGM           ,ONLY: BuildBGMAndIdentifyHaloRegion
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitPEM_LocalElemID,InitPEM_CNElemID,GetMeshMinMax,IdentifyElemAndSideType
USE MOD_Particle_Mesh_Tools    ,ONLY: CalcParticleMeshMetrics,InitElemNodeIDs,InitParticleGeometry,CalcBezierControlPoints
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated,SideSlabNormals,SideSlabIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars ,ONLY: MeasureTrackTime,FastPeriodic,CountNbrOfLostParts,CartesianPeriodic
USE MOD_Particle_Tracking_Vars ,ONLY: NbrOfLostParticles,NbrOfLostParticlesTotal,NbrOfLostParticlesTotal_old
USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLostVecLength,PartStateLost,PartLostDataSize
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod, DisplayLostParticles
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition,DepositionType
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GetRealArray, GETINTFROMSTR
USE MOD_Particle_Vars          ,ONLY: Symmetry
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: SideBoundingBoxVolume
USE MOD_Mesh_Vars              ,ONLY: nSides
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
! TODO
! USE MOD_MPI_Vars               ,ONLY: offsetMPISides_YOUR
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Particle_BGM           ,ONLY: WriteHaloInfo
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /* USE_MPI */
USE MOD_Particle_Mesh_Build    ,ONLY: BuildElementRadiusTria,BuildElemTypeAndBasisTria,BuildEpsOneCell,BuildBCElemDistance
USE MOD_Particle_Mesh_Build    ,ONLY: BuildNodeNeighbourhood,BuildElementOriginShared,BuildElementBasisAndRadius
USE MOD_Particle_Mesh_Build    ,ONLY: BuildSideOriginAndRadius,BuildLinearSideBaseVectors
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!USE MOD_DSMC_Vars              ,ONLY: DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: RefMappingGuessProposal
INTEGER          :: iSample
INTEGER          :: firstSide,lastSide,iSide,SideID
CHARACTER(LEN=2) :: tmpStr
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER          :: ALLOCSTAT
#endif
#ifdef CODE_ANALYZE
! TODO
! REAL             :: dx,dy,dz
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH ...'
IF(ParticleMeshInitIsDone) CALL abort(&
__STAMP__&
, ' Particle-Mesh is already initialized.')

! Potentially curved elements. FIBGM needs to be built on BezierControlPoints rather than NodeCoords to avoid missing elements
IF (TrackingMethod.EQ.TRACING .OR. TrackingMethod.EQ.REFMAPPING) THEN
  ! Bezier elevation now more important than ever, also determines size of FIBGM extent
  BezierElevation = GETINT('BezierElevation')
  NGeoElevated    = NGeo + BezierElevation

  CALL CalcParticleMeshMetrics()

  CALL CalcBezierControlPoints()
END IF

! Mesh min/max must be built on BezierControlPoint for possibly curved elements
CALL GetMeshMinMax()

! Build BGM to Element mapping and identify which of the elements, sides and nodes are in the compute-node local and halo region
CALL BuildBGMAndIdentifyHaloRegion()

#if USE_MPI
CalcHaloInfo = GETLOGICAL('CalcHaloInfo')
IF (CalcHaloInfo) CALL WriteHaloInfo
#endif /*USE_MPI*/

! Initialize mapping function: GetGlobalElemID(1:nComputeNodeTotalElems)
CALL InitGetGlobalElemID()

! Initialize mapping function: GetCNElemID(1:GlobalElemID)
CALL InitGetCNElemID()

! Initialize mapping function: GetGlobalSideID(1:nComputeNodeTotalSides)
CALL InitGetGlobalSideID()

! Initialize mapping function: GetCNSideID(1:GlobalSideID)
CALL InitGetCNSideID()

! Initialize mapping function: PEM%LocalElemID(1:PDM%ParticleVecLength)
CALL InitPEM_LocalElemID()

! Initialize mapping function: PEM%CNElemID(1:PDM%ParticleVecLength)
CALL InitPEM_CNElemID()

CountNbrOfLostParts  = GETLOGICAL('CountNbrOfLostParts')
IF(CountNbrOfLostParts)THEN
  ! Nullify and reset lost parts container after write out
  PartStateLostVecLength = 0

  ! Allocate PartStateLost for a small number of particles and double the array size each time the
  ! maximum is reached
  ALLOCATE(PartStateLost(1:PartLostDataSize,1:10))
  PartStateLost=0.
END IF ! CountNbrOfLostParts
NbrOfLostParticles      = 0
#if USE_LOADBALANCE
! Nullify only once at the beginning of a simulation (not during load balance steps!)
IF(.NOT.PerformLoadBalance)THEN
#endif
  NbrOfLostParticlesTotal = 0
  NbrOfLostParticlesTotal_old = 0
#if USE_LOADBALANCE
END IF
#endif
DisplayLostParticles    = GETLOGICAL('DisplayLostParticles')

#ifdef CODE_ANALYZE
PARTOUT            = GETINT('PartOut','0')
MPIRankOut         = GETINT('MPIRankOut','0')
#endif /*CODE_ANALYZE*/

MeasureTrackTime  = GETLOGICAL('MeasureTrackTime','.FALSE.')
CartesianPeriodic = GETLOGICAL('CartesianPeriodic','.FALSE.')
IF(CartesianPeriodic) FastPeriodic = GETLOGICAL('FastPeriodic','.FALSE.')

! method from xPhysic to parameter space
IF(UseCurveds)THEN ! don't use RefMappingGuess=1, because RefMappingGuess is only best for linear cubical elements
  ! curved elements can be stronger deformed, hence, a better guess can be used
  ! RefMappingGuess 2,3 searches the closest Gauss/CL points of the considered element. This point is used as the initial value for
  ! the mapping. Note, that the position of the CL points can still be advantageous for the initial guess.
  RefMappingGuessProposal=2
  IF(PP_N.GT.NGeo)THEN ! there are more Gauss points within an element then CL-points, Gauss points sample the element finer
    RefMappingGuessProposal=2
  ELSE ! more CL-points than Gauss points, hence, better sampling of the element
    RefMappingGuessProposal=3
  END IF
ELSE
  RefMappingGuessProposal=1 ! default for linear meshes. Guess is exact for cubical, non-twisted elements
END IF
WRITE(tmpStr,'(I2.2)') RefMappingGuessProposal
RefMappingGuess = GETINT('RefMappingGuess',tmpStr)
IF((RefMappingGuess.LT.1).AND.(UseCurveds)) THEN ! Linear intial guess on curved meshes might cause problems.
  SWRITE(UNIT_stdOut,'(A)')' WARNING: read-in [RefMappingGuess=1] when using [UseCurveds=T] may create problems!'
END IF
RefMappingEps   = GETREAL('RefMappingEps','1e-4')

epsInCell       = SQRT(3.0*RefMappingEps)

IF((RefMappingGuess.LT.1).OR.(RefMappingGuess.GT.4))THEN
   CALL abort(&
__STAMP__ &
,'Wrong guessing method for mapping from physical space in reference space.',RefMappingGuess,999.)
END IF

WRITE(tmpStr,'(L1)') (TrackingMethod.EQ.TRIATRACKING)
TriaSurfaceFlux = GETLOGICAL('TriaSurfaceFlux',TRIM(tmpStr))
IF((TrackingMethod.EQ.TRIATRACKING).AND.(.NOT.TriaSurfaceFlux)) THEN
  CALL ABORT(__STAMP__,'TriaSurfaceFlux must be for TriaTracking!')
END IF
IF (Symmetry%Order.LE.2) THEN
  SWRITE(UNIT_stdOut,'(A)') "Surface Flux set to triangle approximation due to Symmetry2D."
  TriaSurfaceFlux = .TRUE.
END IF

! Set logical for building node neighbourhood
SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean','shape_function_adaptive')
    FindNeighbourElems = .TRUE.
  CASE DEFAULT
    FindNeighbourElems = .FALSE.
END SELECT

SELECT CASE(TrackingMethod)
  CASE(TRIATRACKING)
    CALL InitParticleGeometry()
    CALL InitElemNodeIDs()
    ! Compute convex element radius^2
    CALL BuildElementRadiusTria()

    ! Interpolation needs coordinates in reference system
    !IF (DoInterpolation.OR.DSMC%UseOctree) THEN ! use this in future if possible
    IF (DoInterpolation.OR.DoDeposition) THEN
      CALL CalcParticleMeshMetrics()

      CALL BuildElemTypeAndBasisTria()
    END IF ! DoInterpolation.OR.DSMC%UseOctree

    IF (DoDeposition) CALL BuildEpsOneCell()

CASE(TRACING,REFMAPPING)
    IF(TriaSurfaceFlux) CALL InitParticleGeometry()
    IF(FindNeighbourElems) CALL InitElemNodeIDs()

!    CALL CalcParticleMeshMetrics()

!    CALL CalcBezierControlPoints()

#if USE_MPI
    MPISharedSize = INT((3**2*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/3,3,nComputeNodeTotalSides/),SideSlabNormals_Shared_Win,SideSlabNormals_Shared)
    CALL MPI_WIN_LOCK_ALL(0,SideSlabNormals_Shared_Win,IERROR)
    MPISharedSize = INT((6*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalSides/),SideSlabIntervals_Shared_Win,SideSlabIntervals_Shared)
    CALL MPI_WIN_LOCK_ALL(0,SideSlabIntervals_Shared_Win,IERROR)
    MPISharedSize = INT((nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),BoundingBoxIsEmpty_Shared_Win,BoundingBoxIsEmpty_Shared)
    CALL MPI_WIN_LOCK_ALL(0,BoundingBoxIsEmpty_Shared_Win,IERROR)
    firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
    lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
    SideSlabNormals    => SideSlabNormals_Shared
    SideSlabIntervals  => SideSlabIntervals_Shared
    BoundingBoxIsEmpty => BoundingBoxIsEmpty_Shared
    CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
    ALLOCATE(SideSlabNormals(1:3,1:3,1:nNonUniqueGlobalSides) &
            ,SideSlabIntervals(  1:6,1:nNonUniqueGlobalSides) &
            ,BoundingBoxIsEmpty(     1:nNonUniqueGlobalSides) &
            ,STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate SideMetrics arrays!')
    firstSide = 1
    lastSide  = nNonUniqueGlobalSides
#endif /* USE_MPI */
#ifdef CODE_ANALYZE
    ALLOCATE(SideBoundingBoxVolume(nSides))
#endif /*CODE_ANALYZE*/

    IF (BezierElevation.GT.0) THEN
      DO iSide = firstSide,LastSide
        ! ignore sides that are not on the compute node
        ! IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

        SideID = GetGlobalSideID(iSide)

        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,SideID).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,SideID).GT.6) CYCLE

        ! BezierControlPoints are always on nonUniqueGlobalSide
        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,SideID) &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                                       &
                                           ,SideSlabInterVals( 1:6    ,iSide)                                       &
                                           ,BoundingBoxIsEmpty(iSide))
      END DO
    ELSE
      DO iSide=firstSide,LastSide
        ! ignore sides that are not on the compute node
        ! IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

        SideID = GetGlobalSideID(iSide)

        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,SideID).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,SideID).GT.6) CYCLE

        ! BezierControlPoints are always on nonUniqueGlobalSide
        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)                         &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                                       &
                                           ,SideSlabInterVals( 1:6    ,iSide)                                       &
                                           ,BoundingBoxIsEmpty(iSide))
      END DO
  END IF
#if USE_MPI
    CALL MPI_WIN_SYNC(SideSlabNormals_Shared_Win,IERROR)
    CALL MPI_WIN_SYNC(SideSlabIntervals_Shared_Win,IERROR)
    CALL MPI_WIN_SYNC(BoundingBoxIsEmpty_Shared_Win,IERROR)
    CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */
#ifdef CODE_ANALYZE
    ! TODO: bounding box volumes must be calculated for all unique sides.
    !               offsetSideID = ElemInfo_Shared(SideIf
    !               DO iSide=offsetMPISides_YOUR,LastSide
    !                 dx=ABS(SideSlabIntervals(2)-SideSlabIntervals(1))
    !                 dy=ABS(SideSlabIntervals(4)-SideSlabIntervals(3))
    !                 dz=ABS(SideSlabIntervals(6)-SideSlabIntervals(5))
    !                 SideID = SideInfo
    !                 SideBoundingBoxVolume(SideID)=dx*dy*dz
    !               END DO
#endif /*CODE_ANALYZE*/

    ! Compute element bary and element radius for node elements (with halo region)
    CALL BuildElementOriginShared()

    ! Check the side type (planar, bilinear, curved)
    CALL IdentifyElemAndSideType()

    ! Compute the element XiEtaZetaBasis and the radius of the convex hull
    CALL BuildElementBasisAndRadius()

    ! Get basevectors for (bi-)linear sides
    CALL BuildLinearSideBaseVectors()

    IF (TrackingMethod.EQ.REFMAPPING) THEN
      ! Identify BCSides and build side origin and radius
      CALL BuildSideOriginAndRadius()

      ! Identify BCElems
      CALL BuildBCElemDistance()
    END IF

    CALL BuildEpsOneCell()

  CASE DEFAULT
    CALL ABORT(__STAMP__,'Invalid tracking method in particle_mesh.f90!')

END SELECT

! Build mappings UniqueNodeID->CN Element IDs and CN Element ID -> CN Element IDs
IF(FindNeighbourElems) CALL BuildNodeNeighbourhood()

! BezierAreaSample stuff:
IF (TriaSurfaceFlux) THEN
  BezierSampleN = 1
  SurfFluxSideSize=(/1,2/)
ELSE
  WRITE(tmpStr,'(I2.2)') NGeo
  BezierSampleN = GETINT('BezierSampleN',tmpStr)
  SurfFluxSideSize=BezierSampleN
  ALLOCATE(BezierSampleXi(0:BezierSampleN))!,STAT=ALLOCSTAT)
  DO iSample=0,BezierSampleN
    BezierSampleXi(iSample)=-1.+2.0/BezierSampleN*iSample
  END DO
END IF

ParticleMeshInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)') " NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()
DEALLOCATE(NodeCoords)

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE FinalizeParticleMesh()
!===================================================================================================================================
! Deallocates variables for the particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_BGM           ,ONLY: FinalizeBGM
USE MOD_Particle_Mesh_Readin   ,ONLY: FinalizeMeshReadin
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod,Distance,ListDistance,PartStateLost
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
#if USE_MPI
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
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

! Particle mesh readin happens during mesh readin, finalize with gathered routine here
CALL FinalizeMeshReadin()

CALL FinalizeBGM()

SELECT CASE (TrackingMethod)

  ! RefMapping, Tracing
  CASE(REFMAPPING,TRACING)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! BuildSideOriginAndRadius()
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL UNLOCK_AND_FREE(BCSide2SideID_Shared_Win)
      CALL UNLOCK_AND_FREE(SideID2BCSide_Shared_Win)
      CALL UNLOCK_AND_FREE(BCSideMetrics_Shared_Win)
    END IF

#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! CalcParticleMeshMetrics()
      CALL UNLOCK_AND_FREE(XCL_NGeo_Shared_Win)
      CALL UNLOCK_AND_FREE(Elem_xGP_Shared_Win)
      CALL UNLOCK_AND_FREE(dXCL_NGeo_Shared_Win)

      ! CalcBezierControlPoints()
      CALL UNLOCK_AND_FREE(BezierControlPoints3D_Shared_Win)
      IF (BezierElevation.GT.0) CALL UNLOCK_AND_FREE(BezierControlPoints3DElevated_Shared_Win)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! GetSideSlabNormalsAndIntervals() (allocated in particle_mesh.f90)
    CALL UNLOCK_AND_FREE(SideSlabNormals_Shared_Win)
    CALL UNLOCK_AND_FREE(SideSlabIntervals_Shared_Win)
    CALL UNLOCK_AND_FREE(BoundingBoxIsEmpty_Shared_Win)

    ! BuildElementOriginShared()
    CALL UNLOCK_AND_FREE(ElemBaryNGeo_Shared_Win)

    ! IdentifyElemAndSideType()
    CALL UNLOCK_AND_FREE(ElemCurved_Shared_Win)
    CALL UNLOCK_AND_FREE(SideType_Shared_Win)
    CALL UNLOCK_AND_FREE(SideDistance_Shared_Win)
    CALL UNLOCK_AND_FREE(SideNormVec_Shared_Win)

    ! BuildElementBasisAndRadius()
    CALL UNLOCK_AND_FREE(ElemRadiusNGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemRadius2NGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(XiEtaZetaBasis_Shared_Win)
    CALL UNLOCK_AND_FREE(slenXiEtaZetaBasis_Shared_Win)

    ! BuildLinearSideBaseVectors()
    CALL UNLOCK_AND_FREE(BaseVectors0_Shared_Win)
    CALL UNLOCK_AND_FREE(BaseVectors1_Shared_Win)
    CALL UNLOCK_AND_FREE(BaseVectors2_Shared_Win)
    CALL UNLOCK_AND_FREE(BaseVectors3_Shared_Win)
    CALL UNLOCK_AND_FREE(BaseVectorsScale_Shared_Win)

    ! BuildBCElemDistance
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL UNLOCK_AND_FREE(ElemToBCSides_Shared_Win)
      CALL UNLOCK_AND_FREE(SideBCMetrics_Shared_Win)
    END IF

    ! BuildEpsOneCell()
    CALL UNLOCK_AND_FREE(ElemEpsOneCell_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemsJ_Shared_Win)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
    ! BuildSideOriginAndRadius()
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      ADEALLOCATE(BCSide2SideID_Shared)
      ADEALLOCATE(SideID2BCSide_Shared)
      ADEALLOCATE(BCSideMetrics_Shared)
      ADEALLOCATE(BCSide2SideID)
      ADEALLOCATE(SideID2BCSide)
      ADEALLOCATE(BCSideMetrics)
    END IF

#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! CalcParticleMeshMetrics()
      ADEALLOCATE(XCL_NGeo_Array)
      ADEALLOCATE(Elem_xGP_Array)
      ADEALLOCATE(dXCL_NGeo_Array)
      IF(ASSOCIATED(XCL_NGeo_Shared))  NULLIFY(XCL_NGeo_Shared)
      IF(ASSOCIATED(Elem_xGP_Shared))  NULLIFY(Elem_xGP_Shared)
      IF(ASSOCIATED(dXCL_NGeo_Shared)) NULLIFY(dXCL_NGeo_Shared)

      ! CalcBezierControlPoints()
      ADEALLOCATE(BezierControlPoints3D_Shared)
      ADEALLOCATE(BezierControlPoints3DElevated_Shared)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! GetSideSlabNormalsAndIntervals() (allocated in particle_mesh.f90)
    ADEALLOCATE(SideSlabNormals_Shared)
    ADEALLOCATE(SideSlabIntervals_Shared)
    ADEALLOCATE(BoundingBoxIsEmpty_Shared)

    ! BuildElementOriginShared()
    ADEALLOCATE(ElemBaryNGeo_Shared)

    ! IdentifyElemAndSideType()
    ADEALLOCATE(ElemCurved)
    ADEALLOCATE(ElemCurved_Shared)
    ADEALLOCATE(SideType_Shared)
    ADEALLOCATE(SideDistance_Shared)
    ADEALLOCATE(SideNormVec_Shared)

    ! BuildElementBasisAndRadius()
    ADEALLOCATE(ElemRadiusNGeo_Shared)
    ADEALLOCATE(ElemRadius2NGeo_Shared)
    ADEALLOCATE(XiEtaZetaBasis_Shared)
    ADEALLOCATE(slenXiEtaZetaBasis_Shared)

    ! BuildLinearSideBaseVectors()
    ADEALLOCATE(BaseVectors0_Shared)
    ADEALLOCATE(BaseVectors1_Shared)
    ADEALLOCATE(BaseVectors2_Shared)
    ADEALLOCATE(BaseVectors3_Shared)
    ADEALLOCATE(BaseVectorsScale_Shared)

    ! BuildBCElemDistance()
    IF (TrackingMethod.EQ.1) THEN
      ADEALLOCATE(ElemToBCSides_Shared)
      ADEALLOCATE(SideBCMetrics_Shared)
    END IF

    ! BuildEpsOneCell()
    ADEALLOCATE(ElemsJ_Shared)
    ADEALLOCATE(ElemEpsOneCell_Shared)

!  ! Tracing
!  CASE(TRACING)

  ! TriaTracking
  CASE(TRIATRACKING)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! InitParticleGeometry()
    CALL UNLOCK_AND_FREE(ConcaveElemSide_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemSideNodeID_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemMidPoint_Shared_Win)

    ! BuildElementRadiusTria()
    CALL UNLOCK_AND_FREE(ElemBaryNGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemRadius2NGeo_Shared_Win)

    IF(StringBeginsWith(DepositionType,'shape_function'))THEN
      CALL UNLOCK_AND_FREE(ElemRadiusNGeo_Shared_Win)
    END IF

    !IF (DoInterpolation.OR.DSMC%UseOctree) THEN ! use this in future if possible
    IF (DoInterpolation.OR.DoDeposition) THEN
#if USE_LOADBALANCE
      IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
        ! CalcParticleMeshMetrics()
        CALL UNLOCK_AND_FREE(XCL_NGeo_Shared_Win)
        CALL UNLOCK_AND_FREE(Elem_xGP_Shared_Win)
        CALL UNLOCK_AND_FREE(dXCL_NGeo_Shared_Win)
#if USE_LOADBALANCE
      END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

      ! BuildElemTypeAndBasisTria()
      CALL UNLOCK_AND_FREE(ElemCurved_Shared_Win)
      CALL UNLOCK_AND_FREE(XiEtaZetaBasis_Shared_Win)
      CALL UNLOCK_AND_FREE(slenXiEtaZetaBasis_Shared_Win)
    END IF ! DoInterpolation

    ! BuildEpsOneCell()
    IF (DoDeposition) CALL UNLOCK_AND_FREE(ElemsJ_Shared_Win)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! Then, free the pointers or arrays
    ! CalcParticleMeshMetrics
#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ADEALLOCATE(XCL_NGeo_Array)
      ADEALLOCATE(Elem_xGP_Array)
      ADEALLOCATE(dXCL_NGeo_Array)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! CalcParticleMeshMetrics
      IF(ASSOCIATED(XCL_NGeo_Shared))  NULLIFY(XCL_NGeo_Shared)
      IF(ASSOCIATED(Elem_xGP_Shared))  NULLIFY(Elem_xGP_Shared)
      IF(ASSOCIATED(dXCL_NGeo_Shared)) NULLIFY(dXCL_NGeo_Shared)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! InitParticleGeometry
    ADEALLOCATE(ConcaveElemSide_Shared)
    ADEALLOCATE(ElemSideNodeID_Shared)
    ADEALLOCATE(ElemMidPoint_Shared)

    ! BuildElementRadiusTria
    ADEALLOCATE(ElemBaryNGeo_Shared)
    ADEALLOCATE(ElemRadius2NGEO_Shared)
    ADEALLOCATE(ElemRadiusNGeo_Shared)!only shape function

    ! BuildElemTypeAndBasisTria()
    ADEALLOCATE(XiEtaZetaBasis_Shared)
    ADEALLOCATE(slenXiEtaZetaBasis_Shared)
    ADEALLOCATE(ElemCurved)
    ADEALLOCATE(ElemCurved_Shared)

    ! BuildEpsOneCell
    SNULLIFY(ElemsJ)
    ADEALLOCATE(ElemsJ_Shared)
END SELECT

IF(FindNeighbourElems.OR.TrackingMethod.EQ.TRIATRACKING)THEN
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  ! From InitElemNodeIDs
  CALL UNLOCK_AND_FREE(ElemNodeID_Shared_Win)

  !FindNeighbourElems = .FALSE. ! THIS IS SET TO FALSE CURRENTLY in InitParticleMesh()
  ! TODO: fix when FindNeighbourElems is not always set false
  IF(FindNeighbourElems)THEN
    ! From BuildNodeNeighbourhood
    CALL UNLOCK_AND_FREE(NodeToElemMapping_Shared_Win)
    CALL UNLOCK_AND_FREE(NodeToElemInfo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemToElemMapping_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemToElemInfo_Shared_Win)
  END IF ! FindNeighbourElems

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

  ADEALLOCATE(ElemNodeID_Shared)
  IF(FindNeighbourElems)THEN
    ADEALLOCATE(NodeToElemMapping_Shared)
    ADEALLOCATE(NodeToElemInfo_Shared)
    ADEALLOCATE(ElemToElemMapping_Shared)
    ADEALLOCATE(ElemToElemInfo_Shared)
  END IF
END IF
SDEALLOCATE(GEO%DirPeriodicVectors)
SDEALLOCATE(GEO%PeriodicVectors)
SDEALLOCATE(GEO%FIBGM)
SDEALLOCATE(GEO%TFIBGM)

ADEALLOCATE(XiEtaZetaBasis)
ADEALLOCATE(slenXiEtaZetaBasis)
ADEALLOCATE(ElemBaryNGeo)
ADEALLOCATE(ElemRadiusNGeo)
ADEALLOCATE(ElemRadius2NGeo)
ADEALLOCATE(ElemEpsOneCell)
SDEALLOCATE(Distance)
SDEALLOCATE(ListDistance)
SDEALLOCATE(PartStateLost)
SDEALLOCATE(ElemTolerance)
SDEALLOCATE(ElemToGlobalElemID)

ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh


END MODULE MOD_Particle_Mesh
