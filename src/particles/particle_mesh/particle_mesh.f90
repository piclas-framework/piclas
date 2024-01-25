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
CALL prms%CreateLogicalOption( 'DisplayLostParticles' , 'Display position, velocity, species and host element of particles lost during particle tracking (TrackingMethod = triatracking, tracing)','.FALSE.')
CALL prms%CreateLogicalOption( 'CountNbrOfLostParts'&
    , 'Count the number of lost particles during tracking that cannot be found with fallbacks. Additionally, the lost particle '//&
    'information is stored in a PartStateLost*.h5 file. When particles are not found during restart in their host cell '//&
    '(sanity check), they are marked missing and are also written to PartStateLost*.h5 file even if they are re-located '//&
    'on a different processor.','.TRUE.')
CALL prms%CreateIntOption( 'PhotonModeBPO' , 'Output mode to store position, direction, host element etc. of rays/photons in PartStateBoundary.h5 (only radiation transport or ray tracing solver):\n'&
                                             '0: Output nothing to PartStateBoundary.h5\n'&
                                             '1: Output the initial position of the rays and their direction vector\n'&
                                             '2: Output initial position and all calculated intersection points calculated in radtrans tracking\n'&
                                             ,'0')
                                         CALL prms%CreateLogicalOption( 'UsePhotonTriaTracking', 'Activates usage of TriaTracking methods for photon tracking or Bilinear methods (default is True). Can only be selected when ray tracing is actually performed.','.TRUE.')
CALL prms%CreateLogicalOption( 'DoBoundaryParticleOutputRay', 'Activates output of emission particles by ray tracing SEE and ray tracing volume ionization to PartStateBoundary.h5 (with negative species IDs to indicate creation)','.FALSE.')
CALL prms%CreateIntOption(     'PartOut'&
  , 'If compiled with CODE_ANALYZE flag: For This particle number every tracking information is written as STDOUT.','0')
CALL prms%CreateIntOption(     'MPIRankOut'&
  , 'If compiled with CODE_ANALYZE flag: This MPI-Proc writes the tracking information for the defined PartOut.','-1')
CALL prms%CreateLogicalOption( 'MeasureTrackTime'&
  , 'If .TRUE. then the time how long the tracking routines are called are sampled and written for each MPI-Proc.','.FALSE.')
CALL prms%CreateLogicalOption( 'CartesianPeriodic'&
    , ' Simplified treatment for periodic box with Refmapping. Not computation of intersection points at periodic BCs.','.FALSE.')
CALL prms%CreateLogicalOption( 'FastPeriodic'&
  , ' Further simplification by directly moving particle into grid. Instead of moving the particle several times the periodic'//&
    ' displacements, the particle is mapped directly back into the domain. ','.FALSE.')
CALL prms%CreateLogicalOption( 'meshCheckWeirdElements'&
  , 'Abort when weird elements are found: it means that part of the element is turned inside-out. ','.TRUE.')

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
USE MOD_Mesh_Tools             ,ONLY: InitGetGlobalSideID,InitGetCNSideID,GetGlobalSideID,InitElemNodeIDs
USE MOD_Mesh_Vars              ,ONLY: deleteMeshPointer!,NodeCoords
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Mesh_Vars              ,ONLY: useCurveds
#if USE_MPI
USE MOD_Analyze_Vars           ,ONLY: CalcHaloInfo
#endif /*USE_MPI*/
USE MOD_Particle_BGM           ,ONLY: BuildBGMAndIdentifyHaloRegion
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitPEM_LocalElemID,InitPEM_CNElemID,GetMeshMinMax,IdentifyElemAndSideType
USE MOD_Particle_Mesh_Tools    ,ONLY: CalcParticleMeshMetrics,InitParticleGeometry,CalcBezierControlPoints
USE MOD_Particle_Mesh_Tools    ,ONLY: CalcXCL_NGeo
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Tracking_Vars ,ONLY: MeasureTrackTime,FastPeriodic,CountNbrOfLostParts,CartesianPeriodic
USE MOD_Particle_Tracking_Vars ,ONLY: NbrOfLostParticles,NbrOfLostParticlesTotal,NbrOfLostParticlesTotal_old
USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLostVecLength,PartStateLost,PartLostDataSize
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod, DisplayLostParticles
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition,DepositionType
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GetRealArray, GETINTFROMSTR
USE MOD_Particle_Vars          ,ONLY: Symmetry, DoVirtualCellMerge
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
#ifdef CODE_ANALYZE
!USE MOD_Particle_Surfaces_Vars ,ONLY: SideBoundingBoxVolume
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
! TODO
! USE MOD_MPI_Vars               ,ONLY: offsetMPISides_YOUR
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_BGM           ,ONLY: WriteHaloInfo
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_Particle_MPI_Vars      ,ONLY: DoParticleLatencyHiding
#endif /* USE_MPI */
USE MOD_Particle_Mesh_Build    ,ONLY: BuildElementRadiusTria,BuildElemTypeAndBasisTria,BuildEpsOneCell,BuildBCElemDistance
USE MOD_Particle_Mesh_Build    ,ONLY: BuildNodeNeighbourhood,BuildElementOriginShared,BuildElementBasisAndRadius
USE MOD_Particle_Mesh_Build    ,ONLY: BuildSideOriginAndRadius,BuildLinearSideBaseVectors,BuildSideSlabAndBoundingBox
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_PICDepo_Shapefunction_Tools, ONLY:InitShapeFunctionDimensionalty
USE MOD_Particle_Boundary_Init ,ONLY: InitPartStateBoundary
USE MOD_Particle_Boundary_Vars ,ONLY: DoBoundaryParticleOutputHDF5,nSurfSample,DoBoundaryParticleOutputRay
USE MOD_Photon_TrackingVars    ,ONLY: PhotonModeBPO,UsePhotonTriaTracking
USE MOD_RayTracing_Vars        ,ONLY: UseRayTracing
USE MOD_RayTracing_Vars        ,ONLY: PerformRayTracing
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
CHARACTER(LEN=2) :: tmpStr
#if !USE_MPI
INTEGER          :: ALLOCSTAT
#endif
#ifdef CODE_ANALYZE
! TODO
! REAL             :: dx,dy,dz
#endif /*CODE_ANALYZE*/
CHARACTER(3)      :: hilf
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH ...'
IF(ParticleMeshInitIsDone) CALL abort(__STAMP__, ' Particle-Mesh is already initialized.')

WRITE(UNIT=hilf,FMT='(I0)') NGeo
nSurfSample = GETINT('DSMC-nSurfSample',TRIM(hilf))

#if USE_MPI
IF(DoParticleLatencyHiding)THEN
  ! Exchange elements may receive particles during MPI communication and cannot be used for latency hiding
  ALLOCATE(IsExchangeElem(nElems))
  IsExchangeElem = .FALSE.
  CALL AddToElemData(ElementOut,'IsExchangeElem',LogArray=IsExchangeElem)
END IF ! DoParticleLatencyHiding
#endif /*USE_MPI*/

! Check if Bezier control points are required for high-order surface sampling
nSurfSampleAndTriaTracking = .FALSE. ! default
IF((TrackingMethod.EQ.TRIATRACKING).AND.(Symmetry%Order.EQ.3).AND.(nSurfSample.GT.1)) nSurfSampleAndTriaTracking = .TRUE.

! Potentially curved elements. FIBGM needs to be built on BezierControlPoints rather than NodeCoords to avoid missing elements
IF (TrackingMethod.EQ.TRACING .OR. TrackingMethod.EQ.REFMAPPING .OR. nSurfSampleAndTriaTracking .OR. UseRayTracing) THEN
  UseBezierControlPoints = .TRUE.
  ! Bezier elevation now more important than ever, also determines size of FIBGM extent
  BezierElevation = GETINT('BezierElevation')
  NGeoElevated    = NGeo + BezierElevation

  CALL CalcParticleMeshMetrics() ! Required for Elem_xGP_Shared and dXCL_NGeo_Shared
  CALL CalcXCL_NGeo()            ! Required for XCL_NGeo_Shared
  CALL CalcBezierControlPoints() ! Required for BezierControlPoints3D and BezierControlPoints3DElevated (requires XCL_NGeo_Shared)
ELSE
  UseBezierControlPoints = .FALSE.
END IF

! Mesh min/max must be built on BezierControlPoint for possibly curved elements
CALL GetMeshMinMax()

! Set shape function dimension (1D, 2D or 3D)
! This function requires GetMeshMinMax() and values calculated in it are used in BuildBGMAndIdentifyHaloRegion()
IF(StringBeginsWith(DepositionType,'shape_function')) CALL InitShapeFunctionDimensionalty()

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
DisplayLostParticles = GETLOGICAL('DisplayLostParticles')

! Ray tracing information to .h5 for debugging when using the radiation transport model or pure ray tracing
PhotonModeBPO               = GETINT('PhotonModeBPO')
! Use TriaTracking methods for photon tracking or Bilinear methods (default is UsePhotonTriaTracking=T)
IF(PerformRayTracing)THEN
  UsePhotonTriaTracking = GETLOGICAL('UsePhotonTriaTracking')
ELSE
  UsePhotonTriaTracking = .TRUE.
END IF ! PerformRayTracing
! Activate output of emission particles by ray tracing SEE and ray tracing volume ionization to PartStateBoundary.h5 (with negative species IDs to indicate creation
DoBoundaryParticleOutputRay = GETLOGICAL('DoBoundaryParticleOutputRay')
! Check if DoBoundaryParticleOutputHDF5 is already activated and PartStateBoundary therefore already allocated
IF((PhotonModeBPO.GE.1)         .AND.(.NOT.DoBoundaryParticleOutputHDF5)) DoBoundaryParticleOutputHDF5 = .TRUE.
IF((DoBoundaryParticleOutputRay).AND.(.NOT.DoBoundaryParticleOutputHDF5)) DoBoundaryParticleOutputHDF5 = .TRUE.

IF(DoBoundaryParticleOutputHDF5) CALL InitPartStateBoundary()

#ifdef CODE_ANALYZE
PARTOUT            = GETINT('PartOut')
MPIRankOut         = GETINT('MPIRankOut')
#endif /*CODE_ANALYZE*/

MeasureTrackTime  = GETLOGICAL('MeasureTrackTime')
CartesianPeriodic = GETLOGICAL('CartesianPeriodic')
IF(CartesianPeriodic) FastPeriodic = GETLOGICAL('FastPeriodic')

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
  LBWRITE(UNIT_stdOut,'(A)')' WARNING: read-in [RefMappingGuess=1] when using [UseCurveds=T] may create problems!'
END IF
RefMappingEps   = GETREAL('RefMappingEps','1e-4')

epsInCell       = SQRT(3.0*RefMappingEps)

IF((RefMappingGuess.LT.1).OR.(RefMappingGuess.GT.4))THEN
   CALL abort(__STAMP__,'Wrong guessing method for mapping from physical space in reference space.',RefMappingGuess,999.)
END IF

WRITE(tmpStr,'(L1)') (TrackingMethod.EQ.TRIATRACKING)
TriaSurfaceFlux = GETLOGICAL('TriaSurfaceFlux',TRIM(tmpStr))
IF((TrackingMethod.EQ.TRIATRACKING).AND.(.NOT.TriaSurfaceFlux)) THEN
  CALL ABORT(__STAMP__,'TriaSurfaceFlux must be for TriaTracking!')
END IF
IF (Symmetry%Order.LE.2) THEN
  LBWRITE(UNIT_stdOut,'(A)') "Surface Flux set to triangle approximation due to Symmetry2D."
  TriaSurfaceFlux = .TRUE.
END IF

! Set logical for building node neighbourhood
FindNeighbourElems = .FALSE.
! PIC deposition types require the neighbourhood
SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean','shape_function_adaptive')
    FindNeighbourElems = .TRUE.
END SELECT
! Rotational periodic BC requires the neighbourhood to add elements of the BC nodes
IF(PartBound%UseRotPeriodicBC) FindNeighbourElems = .TRUE.

IF(DoVirtualCellMerge) FindNeighbourElems = .TRUE.

! Build ConcaveElemSide_Shared, ElemSideNodeID_Shared, ElemMidPoint_Shared
CALL InitParticleGeometry()

SELECT CASE(TrackingMethod)
  CASE(TRIATRACKING)
    CALL InitElemNodeIDs()
    ! Compute convex element radius^2
    CALL BuildElementRadiusTria() ! Required for ElemBaryNGeo_Shared, ElemRadius2NGEO_Shared, ElemRadiusNGEO_Shared (only for shape function)

    ! Interpolation needs coordinates in reference system
    !IF (DoInterpolation.OR.DSMC%UseOctree) THEN ! use this in future if possible
    IF (DoInterpolation.OR.DoDeposition.OR.UseRayTracing) THEN
      ! Do not call these functions twice. This is already done above
      IF(.NOT.UseBezierControlPoints)THEN
        CALL CalcParticleMeshMetrics()   ! Required for Elem_xGP_Shared and dXCL_NGeo_Shared
        CALL CalcXCL_NGeo()              ! Required for XCL_NGeo_Shared
      END IF ! .NOT.UseBezierControlPoints
      CALL BuildElemTypeAndBasisTria() ! Required for ElemCurved_Shared, XiEtaZetaBasis_Shared, slenXiEtaZetaBasis_Shared. Needs XCL_NGeo_Shared
    END IF ! DoInterpolation.OR.DSMC%UseOctree

    IF (DoDeposition) CALL BuildEpsOneCell()

    IF(.NOT.UsePhotonTriaTracking)THEN
      ! Build stuff required for bilinear tracing algorithms
      CALL BuildSideSlabAndBoundingBox() ! Required for SideSlabNormals_Shared, SideSlabIntervals_Shared, BoundingBoxIsEmpty_Shared

      ! Check the side type (planar, bilinear, curved)
      CALL IdentifyElemAndSideType() ! Builds ElemCurved_Shared, SideType_Shared, SideDistance_Shared, SideNormVec_Shared

      ! Get basevectors for (bi-)linear sides
      CALL BuildLinearSideBaseVectors() ! Required for BaseVectors0_Shared, BaseVectors1_Shared, BaseVectors2_Shared, BaseVectors3_Shared, BaseVectorsScale_Shared
    END IF ! UsePhotonTriaTracking

  CASE(TRACING,REFMAPPING)
    ! Build stuff required for tracing algorithms
    CALL BuildSideSlabAndBoundingBox() ! Required for SideSlabNormals_Shared, SideSlabIntervals_Shared, BoundingBoxIsEmpty_Shared

    ! ElemNodeID_Shared required
    IF(FindNeighbourElems) CALL InitElemNodeIDs()

    ! Compute element bary and element radius for node elements (with halo region)
    CALL BuildElementOriginShared()

    ! Check the side type (planar, bilinear, curved)
    CALL IdentifyElemAndSideType() ! Required for  ElemCurved_Shared, SideType_Shared, SideDistance_Shared, SideNormVec_Shared

    ! Compute the element XiEtaZetaBasis and the radius of the convex hull
    CALL BuildElementBasisAndRadius() ! Required for ElemRadiusNGeo_Shared, ElemRadius2NGeo_Shared, XiEtaZetaBasis_Shared, slenXiEtaZetaBasis_Shared

    ! Get basevectors for (bi-)linear sides
    CALL BuildLinearSideBaseVectors() ! Required for BaseVectors0_Shared, BaseVectors1_Shared, BaseVectors2_Shared, BaseVectors3_Shared, BaseVectorsScale_Shared

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

LBWRITE(UNIT_stdOut,'(A)') " InitParticleMesh: NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()
!DEALLOCATE(NodeCoords)

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE FinalizeParticleMesh()
!===================================================================================================================================
! Deallocates variables for the particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_RayTracing_Vars        ,ONLY: UseRayTracing
#if USE_MPI
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
#endif /*USE_MPI*/
USE MOD_Particle_BGM           ,ONLY: FinalizeBGM
USE MOD_Mesh_ReadIn            ,ONLY: FinalizeMeshReadin
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod,Distance,ListDistance,PartStateLost
#if USE_MPI
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
#else
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
USE MOD_Photon_TrackingVars    ,ONLY: UsePhotonTriaTracking
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
CALL FinalizeMeshReadin(2)

CALL FinalizeBGM()

SELECT CASE (TrackingMethod)

  ! =============================================================================
  ! RefMapping, Tracing
  CASE(REFMAPPING,TRACING)
  ! =============================================================================
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! InitParticleGeometry()
    CALL UNLOCK_AND_FREE(ConcaveElemSide_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemSideNodeID_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemMidPoint_Shared_Win)

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

    ! BuildSideSlabAndBoundingBox() builds SideSlabNormals_Shared, SideSlabIntervals_Shared, BoundingBoxIsEmpty_Shared
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

    ! BuildSideSlabAndBoundingBox() builds SideSlabNormals_Shared, SideSlabIntervals_Shared, BoundingBoxIsEmpty_Shared
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

  ! =============================================================================
  ! TriaTracking
  CASE(TRIATRACKING)
  ! =============================================================================
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! InitParticleGeometry()
    CALL UNLOCK_AND_FREE(ConcaveElemSide_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemSideNodeID_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemMidPoint_Shared_Win)

    IF(.NOT.UsePhotonTriaTracking)THEN
      ! GetSideSlabNormalsAndIntervals()
      CALL UNLOCK_AND_FREE(SideSlabNormals_Shared_Win)
      CALL UNLOCK_AND_FREE(SideSlabIntervals_Shared_Win)
      CALL UNLOCK_AND_FREE(BoundingBoxIsEmpty_Shared_Win)

      ! IdentifyElemAndSideType()
      CALL UNLOCK_AND_FREE(SideType_Shared_Win)
      CALL UNLOCK_AND_FREE(SideDistance_Shared_Win)
      CALL UNLOCK_AND_FREE(SideNormVec_Shared_Win)

      ! BuildLinearSideBaseVectors()
      CALL UNLOCK_AND_FREE(BaseVectors0_Shared_Win)
      CALL UNLOCK_AND_FREE(BaseVectors1_Shared_Win)
      CALL UNLOCK_AND_FREE(BaseVectors2_Shared_Win)
      CALL UNLOCK_AND_FREE(BaseVectors3_Shared_Win)
      CALL UNLOCK_AND_FREE(BaseVectorsScale_Shared_Win)
    END IF ! .NOT.UsePhotonTriaTracking

    ! BuildElementRadiusTria()
    CALL UNLOCK_AND_FREE(ElemBaryNGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemRadius2NGeo_Shared_Win)

    IF(StringBeginsWith(DepositionType,'shape_function'))THEN
      CALL UNLOCK_AND_FREE(ElemRadiusNGeo_Shared_Win)
    END IF

    !IF (DoInterpolation.OR.DSMC%UseOctree) THEN ! use this in future if possible
    IF (DoInterpolation.OR.DoDeposition.OR.UseRayTracing.OR.nSurfSampleAndTriaTracking) THEN
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

      IF (DoInterpolation.OR.DoDeposition.OR.UseRayTracing) THEN
        ! BuildElemTypeAndBasisTria()
        CALL UNLOCK_AND_FREE(ElemCurved_Shared_Win)
        CALL UNLOCK_AND_FREE(XiEtaZetaBasis_Shared_Win)
        CALL UNLOCK_AND_FREE(slenXiEtaZetaBasis_Shared_Win)
      END IF ! DoInterpolation.OR.DoDeposition.OR.UseRayTracing
    END IF ! DoInterpolation.OR.DoDeposition.OR.UseRayTracing.OR.nSurfSampleAndTriaTracking

    ! BuildEpsOneCell()
    IF (DoDeposition) CALL UNLOCK_AND_FREE(ElemsJ_Shared_Win)

#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      IF(UseBezierControlPoints)THEN
        ! CalcBezierControlPoints()
        CALL UNLOCK_AND_FREE(BezierControlPoints3D_Shared_Win)
        IF (BezierElevation.GT.0) CALL UNLOCK_AND_FREE(BezierControlPoints3DElevated_Shared_Win)
      END IF ! UseBezierControlPoints
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! Then, free the pointers or arrays
    ! CalcParticleMeshMetrics
#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ADEALLOCATE(XCL_NGeo_Array)
      ADEALLOCATE(Elem_xGP_Array)
      ADEALLOCATE(dXCL_NGeo_Array)

      ! CalcBezierControlPoints()
      ADEALLOCATE(BezierControlPoints3D_Shared)
      ADEALLOCATE(BezierControlPoints3DElevated_Shared)
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

    IF(.NOT.UsePhotonTriaTracking)THEN
      ! GetSideSlabNormalsAndIntervals()
      ADEALLOCATE(SideSlabNormals_Shared)
      ADEALLOCATE(SideSlabIntervals_Shared)
      ADEALLOCATE(BoundingBoxIsEmpty_Shared)

      ! IdentifyElemAndSideType()
      ADEALLOCATE(SideType_Shared)
      ADEALLOCATE(SideDistance_Shared)
      ADEALLOCATE(SideNormVec_Shared)

      ! BuildLinearSideBaseVectors()
      ADEALLOCATE(BaseVectors0_Shared)
      ADEALLOCATE(BaseVectors1_Shared)
      ADEALLOCATE(BaseVectors2_Shared)
      ADEALLOCATE(BaseVectors3_Shared)
      ADEALLOCATE(BaseVectorsScale_Shared)
    END IF ! .NOT.UsePhotonTriaTracking

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
#if USE_LOADBALANCE
  !IF(.NOT.(PerformLoadBalance.AND.DoDeposition.AND.DoDielectricSurfaceCharge))THEN
  ! Note that no inquiry for DoDeposition is made here because the surface charging container is to be preserved
  IF(.NOT.(PerformLoadBalance.AND.DoDielectricSurfaceCharge))THEN
#endif /*USE_LOADBALANCE*/
    ! From InitElemNodeIDs
    CALL UNLOCK_AND_FREE(ElemNodeID_Shared_Win)
#if USE_LOADBALANCE
  END IF ! .NOT.(PerformLoadBalance.AND.DoDeposition.AND.DoDielectricSurfaceCharge)
#endif /*USE_LOADBALANCE*/

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

#if USE_LOADBALANCE
  !IF(.NOT.(PerformLoadBalance.AND.DoDeposition.AND.DoDielectricSurfaceCharge))THEN
  ! Note that no inquiry for DoDeposition is made here because the surface charging container is to be preserved
  IF(.NOT.(PerformLoadBalance.AND.DoDielectricSurfaceCharge))THEN
#endif /*USE_LOADBALANCE*/
      ADEALLOCATE(ElemNodeID_Shared)
#if USE_LOADBALANCE
  END IF ! .NOT.(PerformLoadBalance.AND.DoDeposition.AND.DoDielectricSurfaceCharge)
#endif /*USE_LOADBALANCE*/
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
SDEALLOCATE(GEO%XMinMax)
#if USE_MPI
SDEALLOCATE(IsExchangeElem)
#endif /*USE_MPI*/

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

! InitParticleGeometry
ADEALLOCATE(ConcaveElemSide_Shared)
ADEALLOCATE(ElemSideNodeID_Shared)
ADEALLOCATE(ElemMidPoint_Shared)

! Load Balance
!#if !USE_LOADBALANCE
!SDEALLOCATE(ElemTime)
!#endif /* !USE_LOADBALANCE */

ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh


END MODULE MOD_Particle_Mesh
