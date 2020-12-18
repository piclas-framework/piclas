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

INTERFACE GetMeshMinMax
  MODULE PROCEDURE GetMeshMinMax
END INTERFACE

INTERFACE MapRegionToElem
  MODULE PROCEDURE MapRegionToElem
END INTERFACE

INTERFACE MarkAuxBCElems
  MODULE PROCEDURE MarkAuxBCElems
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
PUBLIC::MapRegionToElem
PUBLIC::MarkAuxBCElems
PUBLIC::GetMeshMinMax
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

CALL prms%CreateLogicalOption( 'Write-Tria-DebugMesh'&
  , 'Writes per proc triangulated Surfacemesh used for Triatracking. Requires TriaTracking=T.'&
  ,'.FALSE.')
CALL prms%CreateLogicalOption( 'TriaSurfaceFlux'&
  , 'Using Triangle-aproximation [T] or (bi-)linear and bezier (curved) description [F] of sides for surfaceflux.'//&
  ' Default is set to TriaTracking')
CALL prms%CreateLogicalOption( 'Write-TriaSurfaceFlux-DebugMesh'&
  , 'Writes per proc triangulated Surfacemesh used for TriaSurfaceFlux. Requires TriaSurfaceFlux=T.'&
  ,'.FALSE.')

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


! Background mesh init variables
!CALL prms%CreateLogicalOption( 'printMPINeighborWarnings'&
    !,  ' Print warning if the MPI-Halo-region between to procs are not overlapping. Only one proc find the other in halo ' &
    !,'.FALSE.')
CALL prms%CreateLogicalOption( 'CalcHaloInfo',         'Output halo element information to ElemData for each processor'//&
                                                       ' "MyRank_ElemHaloInfo"\n'//&
                                                       ' ElemHaloInfo = -1            : element not in list\n'//&
                                                       '              = 0             : halo elements\n'//&
                                                       '              = 1 to PP_nElems: local elements','.FALSE.')
CALL prms%CreateLogicalOption( 'printBezierControlPointsWarnings'&
    ,  ' Print warning if MINVAL(BezierControlPoints3D(iDir,:,:,newSideID)) and global boundaries are too close ' &
    ,'.FALSE.')

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
USE MOD_Mesh_Tools             ,ONLY: InitGetGlobalSideID,InitGetCNSideID
USE MOD_Mesh_Vars              ,ONLY: deleteMeshPointer,NodeCoords
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Mesh_Vars              ,ONLY: useCurveds
USE MOD_Analyze_Vars           ,ONLY: CalcHaloInfo
USE MOD_Particle_BGM           ,ONLY: BuildBGMAndIdentifyHaloRegion
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitPEM_LocalElemID,InitPEM_CNElemID
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated,SideSlabNormals,SideSlabIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars ,ONLY: MeasureTrackTime,FastPeriodic,CountNbrOfLostParts,CartesianPeriodic
USE MOD_Particle_Tracking_Vars ,ONLY: NbrOfLostParticles,NbrOfLostParticlesTotal
USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLostVecLength,PartStateLost
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod, DisplayLostParticles
!USE MOD_Particle_Tracking_Vars ,ONLY: WriteTriaDebugMesh
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
INTEGER          :: firstSide,lastSide,iSide
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
IF (CalcHaloInfo) THEN
  CALL WriteHaloInfo
END IF
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
  ALLOCATE(PartStateLost(1:14,1:10))
  PartStateLost=0.
END IF ! CountNbrOfLostParts
NbrOfLostParticles      = 0
NbrOfLostParticlesTotal = 0
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

SELECT CASE(TrackingMethod)
  CASE(TRIATRACKING)
    CALL InitParticleGeometry()
    CALL InitElemNodeIDs()
    ! Compute convex element radius^2
    CALL BuildElementRadiusTria()

    ! Interpolation needs coordinates in reference system
    !IF (DoInterpolation.OR.DSMC%UseOctree) THEN ! use this in future if possible
    IF (DoInterpolation) THEN
      CALL CalcParticleMeshMetrics()

      CALL BuildElemTypeAndBasisTria()
    END IF ! DoInterpolation.OR.DSMC%UseOctree

    IF (DoDeposition) CALL BuildEpsOneCell()

CASE(TRACING,REFMAPPING)
    IF(TriaSurfaceFlux) CALL InitParticleGeometry()
    IF(TRIM(DepositionType).EQ.'shape_function_adaptive') CALL InitElemNodeIDs()

!    CALL CalcParticleMeshMetrics()

!    CALL CalcBezierControlPoints()

#if USE_MPI
    MPISharedSize = INT((3**2*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/3,3,nNonUniqueGlobalSides/),SideSlabNormals_Shared_Win,SideSlabNormals_Shared)
    CALL MPI_WIN_LOCK_ALL(0,SideSlabNormals_Shared_Win,IERROR)
    MPISharedSize = INT((6*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/6,nNonUniqueGlobalSides/),SideSlabIntervals_Shared_Win,SideSlabIntervals_Shared)
    CALL MPI_WIN_LOCK_ALL(0,SideSlabIntervals_Shared_Win,IERROR)
    MPISharedSize = INT((nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
    CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),BoundingBoxIsEmpty_Shared_Win,BoundingBoxIsEmpty_Shared)
    CALL MPI_WIN_LOCK_ALL(0,BoundingBoxIsEmpty_Shared_Win,IERROR)
    firstSide = INT(REAL (myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
    lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
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
      DO iSide=firstSide,LastSide
        ! ignore sides that are not on the compute node
        IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE

        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide)     &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                          &
                                           ,SideSlabInterVals( 1:6    ,iSide)                          &
                                           ,BoundingBoxIsEmpty(iSide))
      END DO
    ELSE
      DO iSide=firstSide,LastSide
        ! ignore sides that are not on the compute node
        IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE

        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)             &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                          &
                                           ,SideSlabInterVals( 1:6    ,iSide)                          &
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
    CALL GetLinearSideBaseVectors()

    IF (TrackingMethod.EQ.REFMAPPING) THEN
      ! Identify BCSides and build side origin and radius
      CALL GetBCSidesAndOrgin()

      ! Identify BCElems
      CALL BuildBCElemDistance()
    END IF

    CALL BuildEpsOneCell()

  CASE DEFAULT
    CALL ABORT(__STAMP__,'Invalid tracking method in particle_mesh.f90!')

END SELECT

! Build mappings UniqueNodeID->CN Element IDs and CN Element ID -> CN Element IDs
FindNeighbourElems = .FALSE.
IF(FindNeighbourElems.OR.TRIM(DepositionType).EQ.'shape_function_adaptive') CALL BuildNodeNeighbourhood()

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


SUBROUTINE CalcParticleMeshMetrics()
!===================================================================================================================================
!> calculates XCL_Ngeo and dXCL_Ngeo for global mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis                  ,ONLY: BuildBezierVdm,BuildBezierDMat
USE MOD_Basis                  ,ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
!USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
!USE MOD_Interpolation          ,ONLY: GetDerivativeMatrix
!USE MOD_Interpolation          ,ONLY: GetVandermonde
!USE MOD_Interpolation_Vars     ,ONLY: NodeType,NodeTypeCL,NodeTypeVISU
!USE MOD_Mesh_Vars              ,ONLY: useCurveds
USE MOD_Mesh_Vars              ,ONLY: Elem_xGP
USE MOD_Mesh_Vars              ,ONLY: NGeo,XCL_NGeo,wBaryCL_NGeo,XiCL_NGeo,dXCL_NGeo,Xi_NGeo
USE MOD_Mesh_Vars              ,ONLY: wBaryCL_NGeo1,Vdm_CLNGeo1_CLNGeo,XiCL_NGeo1,nElems
!USE MOD_Mesh_Tools             ,ONLY: GetCNElemID,GetGlobalElemID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: Vdm_Bezier,sVdm_Bezier,D_Bezier
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems,offsetElem
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                        :: iElem
REAL                           :: Vdm_NGeo_CLNGeo(0:NGeo,0:NGeo)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================

! small wBaryCL_NGEO
ALLOCATE( wBaryCL_NGeo1(            0:1)    ,&
          XiCL_NGeo1(               0:1)    ,&
          Vdm_CLNGeo1_CLNGeo(0:NGeo,0:1)    ,&
! new for curved particle sides
          Vdm_Bezier(        0:NGeo,0:NGeo) ,&
          sVdm_Bezier(       0:NGeo,0:NGeo) ,&
          D_Bezier(          0:NGeo,0:NGeo))

CALL ChebyGaussLobNodesAndWeights(1,XiCL_NGeo1)
CALL BarycentricWeights(1,XiCL_NGeo1,wBaryCL_NGeo1)
CALL InitializeVandermonde(1, NGeo,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo ,Vdm_CLNGeo1_CLNGeo)
! initialize Vandermonde for Bezier basis surface representation (particle tracking with curved elements)
CALL BuildBezierVdm(NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier) !CHANGETAG
CALL BuildBezierDMat(NGeo,Xi_NGeo,D_Bezier)
CALL InitializeVandermonde(NGeo,NGeo,wBaryCL_NGeo,Xi_NGeo,XiCL_NGeo,Vdm_NGeo_CLNGeo)

#if USE_LOADBALANCE
! XCL and dXCL are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
MPISharedSize = INT((3*(NGeo+1)**3*nGlobalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/), XCL_NGeo_Shared_Win,XCL_NGeo_Array)
MPISharedSize = INT((3*(PP_N+1)**3*nGlobalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (PP_N+1)*(PP_N+1)*(PP_N+1)*nGlobalElems/), Elem_xGP_Shared_Win,Elem_xGP_Array)
MPISharedSize = INT((3*3*(NGeo+1)**3*nGlobalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/),dXCL_NGeo_Shared_Win,dXCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,dXCL_NGeo_Shared_Win,IERROR)
XCL_NGeo_Shared (1:3    ,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => XCL_NGeo_Array
Elem_xGP_Shared (1:3    ,0:PP_N,0:PP_N,0:PP_N,1:nGlobalElems) => Elem_xGP_Array
dXCL_NGeo_Shared(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => dXCL_NGeo_Array

! Copy local XCL and dXCL into shared memory
!IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  DO iElem = 1,nElems
    XCL_NGeo_Shared (:  ,:,:,:,offsetElem+iElem) = XCL_NGeo (:  ,:,:,:,iElem)
    Elem_xGP_Shared (:  ,:,:,:,offsetElem+iElem) = Elem_xGP (:  ,:,:,:,iElem)
    dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
  END DO ! iElem = 1, nElems
!ELSE
!  DO iElem = 1, nElems
!    XCL_NGeo_Shared (:,  :,:,:,offsetElem+iElem) = XCL_NGeo (:,  :,:,:,iElem)
!    Elem_xGP_Shared (:,  :,:,:,offsetElem+iElem) = Elem_xGP (:,  :,:,:,iElem)
!    dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
!  END DO ! iElem = 1, nElems
!END IF

! Communicate XCL and dXCL between compute node roots instead of calculating globally
CALL MPI_WIN_SYNC(XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

IF (nComputeNodeProcessors.NE.nProcessors_Global .AND. myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , XCL_NGeo_Shared               &
                     , 3*(NGeo+1)**3*recvcountElem   &
                     , 3*(NGeo+1)**3*displsElem      &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)

  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , Elem_xGP_Shared               &
                     , 3*(PP_N+1)**3*recvcountElem   &
                     , 3*(PP_N+1)**3*displsElem      &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)

  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , dXCL_NGeo_Shared              &
                     , 3*3*(NGeo+1)**3*recvcountElem &
                     , 3*3*(NGeo+1)**3*displsElem    &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)
END IF

!nComputeNodeHaloElems = nComputeNodeTotalElems - nComputeNodeElems
!IF (nComputeNodeHaloElems.GT.nComputeNodeProcessors) THEN
!  firstHaloElem = INT(REAL( myComputeNodeRank   *nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))+1
!  lastHaloElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))
!ELSE
!  firstHaloElem = myComputeNodeRank + 1
!  IF (myComputeNodeRank.LT.nComputeNodeHaloElems) THEN
!    lastHaloElem = myComputeNodeRank + 1
!  ELSE
!    lastHaloElem = 0
!  END IF
!END IF
!
!! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!!       Important for curved meshes if NGeo<N, no effect for N>=NGeo
!CALL GetVandermonde(    NGeo, NodeTypeVISU, PP_N, NodeTypeCL, Vdm_EQNGeo_CLN,  modal=.FALSE.)
!CALL GetVandermonde(    PP_N, NodeTypeCL  , PP_N, NodeType  , Vdm_CLNloc_N,    modal=.FALSE.)
!
!!Transform from EQUI_NGeo to solution points on Nloc
!Vdm_EQNGeo_CLN = MATMUL(Vdm_CLNloc_N,Vdm_EQNGeo_CLN)
!
!! Build XCL and dXCL for compute node halo region (each proc of compute-node build only its fair share)
!  CALL GetDerivativeMatrix(NGeo  , NodeTypeCL  , DCL_Ngeo)
!
!  DO iElem = firstHaloElem, lastHaloElem
!    ElemID = GetGlobalElemID(nComputeNodeElems+iElem)
!!    firstNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+1
!    firstNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)
!!    nodeID = 0
!    nodeID = 1
!    IF (useCurveds) THEN
!      DO k = 0, NGeo; DO j = 0, NGeo; DO i = 0, NGeo
!        NodeCoordstmp(:,i,j,k) = NodeCoords_Shared(:,firstNodeID+NodeID)
!        nodeID = nodeID + 1
!      END DO; END DO; END DO ! i = 0, NGeo
!    ELSE
!      NodeCoordstmp(:,0,0,0) = NodeCoords_Shared(:,firstNodeID+1)
!      NodeCoordstmp(:,1,0,0) = NodeCoords_Shared(:,firstNodeID+2)
!      NodeCoordstmp(:,0,1,0) = NodeCoords_Shared(:,firstNodeID+3)
!      NodeCoordstmp(:,1,1,0) = NodeCoords_Shared(:,firstNodeID+4)
!      NodeCoordstmp(:,0,0,1) = NodeCoords_Shared(:,firstNodeID+5)
!      NodeCoordstmp(:,1,0,1) = NodeCoords_Shared(:,firstNodeID+6)
!      NodeCoordstmp(:,0,1,1) = NodeCoords_Shared(:,firstNodeID+7)
!      NodeCoordstmp(:,1,1,1) = NodeCoords_Shared(:,firstNodeID+8)
!    END IF
!    CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_NGeo_CLNGeo,NodeCoordstmp,XCL_NGeo_Shared(:,:,:,:,nComputeNodeElems+iElem))
!    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_EQNGeo_CLN ,NodeCoordstmp,Elem_xGP_Shared(:,:,:,:,nComputeNodeElems+iElem))
!
!    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
!      ! Matrix-vector multiplication
!      DO ll=0,Ngeo
!        dXCL_NGeo_Shared(1,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(1,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(i,ll) &
!                                                            *  XCL_NGeo_Shared(: ,ll,j,k,nComputeNodeElems+iElem)
!        dXCL_NGeo_Shared(2,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(2,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(j,ll) &
!                                                            *  XCL_NGeo_Shared(: ,i,ll,k,nComputeNodeElems+iElem)
!        dXCL_NGeo_Shared(3,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(3,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(k,ll) &
!                                                            *  XCL_NGeo_Shared(: ,i,j,ll,nComputeNodeElems+iElem)
!      END DO
!    END DO; END DO; END DO !i,j,k=0,Ngeo
!    END DO ! iElem = firstHaloElem, lastHaloElem

CALL MPI_WIN_SYNC(XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
XCL_NGeo_Shared  => XCL_NGeo
Elem_xGP_Shared  => Elem_xGP
dXCL_NGeo_Shared => dXCL_NGeo
#endif /*USE_MPI*/

END SUBROUTINE CalcParticleMeshMetrics


SUBROUTINE CalcBezierControlPoints()
!===================================================================================================================================
!> calculate the Bezier control point (+elevated) for shared compute-node mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Mappings               ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared
!USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: GetBezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3DElevated,BezierElevation
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
USE MOD_Particle_Mesh_Vars     ,ONLY: BezierControlPoints3D_Shared,BezierControlPoints3D_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: BezierControlPoints3DElevated_Shared,BezierControlPoints3DElevated_Shared_Win
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
!USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iSide,ilocSide
INTEGER                        :: SideID
INTEGER                        :: firstElem,lastElem,firstSide,lastSide
INTEGER                        :: p,q,pq(2)
REAL                           :: tmp(3,0:NGeo,0:NGeo)
REAL                           :: tmp2(3,0:NGeo,0:NGeo)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER                        :: ALLOCSTAT
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_LOADBALANCE
! BezierControlPoints are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

SWRITE(UNIT_stdOut,'(A)') ' CALCULATING BezierControlPoints ...'

! Build BezierControlPoints3D (compute-node local+halo)
#if USE_MPI
MPISharedSize = INT((3*(NGeo+1)**2*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared(MPISharedSize,(/3*(NGeo+1)*(NGeo+1)*nNonUniqueGlobalSides/),BezierControlPoints3D_Shared_Win,BezierControlPoints3D_Shared)
CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3D_Shared_Win,IERROR)
BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nNonUniqueGlobalSides) => BezierControlPoints3D_Shared
IF (myComputeNodeRank.EQ.0) THEN
  BezierControlPoints3D         = 0.
END IF
IF (BezierElevation.GT.0) THEN
  MPISharedSize = INT((3*(NGeoElevated+1)**2*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3*(NGeoElevated+1)*(NGeoElevated+1)*nNonUniqueGlobalSides/),BezierControlPoints3DElevated_Shared_Win,BezierControlPoints3DElevated_Shared)
  CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3DElevated_Shared_Win,IERROR)
  BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nNonUniqueGlobalSides) => BezierControlPoints3DElevated_Shared
  IF (myComputeNodeRank.EQ.0) THEN
    BezierControlPoints3DElevated = 0.
  END IF
END IF
#else
ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nNonUniqueGlobalSides) &
        ,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3D!')
BezierControlPoints3D         = 0.

IF (BezierElevation.GT.0) THEN
  ALLOCATE(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nNonUniqueGlobalSides) &
          ,STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3DElevated!')
  BezierControlPoints3DElevated = 0.
END IF
#endif /*USE_MPI*/

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
IF (BezierElevation.GT.0) THEN
  CALL MPI_WIN_SYNC(BezierControlPoints3DElevated_Shared_Win,IERROR)
END IF
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank*   nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nGlobalElems)/REAL(nComputeNodeProcessors))
!firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
!lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
firstSide = INT(REAL (myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! iElem is CN elem
DO iElem = firstElem,lastElem
  DO ilocSide=1,6
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , 0    , :    , :   ,iElem )
    CASE(XI_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , NGeo , :    , :   ,iElem )
    CASE(ETA_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , 0    , :   ,iElem )
    CASE(ETA_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , NGeo , :   ,iElem )
    CASE(ZETA_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , :    , 0   ,iElem )
    CASE(ZETA_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , :    , NGeo,iElem )
    END SELECT
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,tmp,tmp2)

    ! get global SideID of local side
    SideID = GetGlobalNonUniqueSideID(iElem,iLocSide)

    DO q=0,NGeo; DO p=0,NGeo
      ! turn into right hand system of side
      pq = CGNS_SideToVol2(NGeo,p,q,iLocSide)
      BezierControlPoints3D(1:3,p,q,SideID) = tmp2(1:3,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! ilocSide=1,6
END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! Calculate elevated BezierControlPoints
IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,LastSide
    ! Ignore small mortar sides attached to big mortar sides
    IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE
    ! Indices in shared arrays are shifted by 1
    CALL GetBezierControlPoints3DElevated( NGeo,NGeoElevated                                                       &
                                         , BezierControlPoints3D        (1:3,0:NGeo        ,0:NGeo        ,iSide)  &
                                         , BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide))
  END DO

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3DElevated_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/
END IF

END SUBROUTINE CalcBezierControlPoints


SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
! Subroutine for particle geometry initialization (GEO container)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_ReadInTools
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,FirstElem,LastElem,GlobalElemID
INTEGER            :: GlobalSideID,nlocSides,localSideID,iLocSide
INTEGER            :: iNode
INTEGER            :: nStart, NodeNum
INTEGER            :: NodeMap(4,6)
REAL               :: A(3,3),detcon
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
INTEGER            :: CornerNodeIDswitch(8)
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE GEOMETRY INFORMATION...'

! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
CornerNodeIDswitch(1)=1
CornerNodeIDswitch(2)=(Ngeo+1)
CornerNodeIDswitch(3)=(Ngeo+1)**2
CornerNodeIDswitch(4)=(Ngeo+1)*Ngeo+1
CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2
CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1

! New crazy corner node switch (philipesque)
ASSOCIATE(CNS => CornerNodeIDswitch )
  ! CGNS Mapping
  NodeMap(:,1)=(/CNS(1),CNS(4),CNS(3),CNS(2)/)
  NodeMap(:,2)=(/CNS(1),CNS(2),CNS(6),CNS(5)/)
  NodeMap(:,3)=(/CNS(2),CNS(3),CNS(7),CNS(6)/)
  NodeMap(:,4)=(/CNS(3),CNS(4),CNS(8),CNS(7)/)
  NodeMap(:,5)=(/CNS(1),CNS(5),CNS(8),CNS(4)/)
  NodeMap(:,6)=(/CNS(5),CNS(6),CNS(7),CNS(8)/)

#if USE_MPI
  MPISharedSize = INT(6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),ConcaveElemSide_Shared_Win,ConcaveElemSide_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ConcaveElemSide_Shared_Win,IERROR)
  firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
  lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

  MPISharedSize = INT(4*6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/4,6,nComputeNodeTotalElems/),ElemSideNodeID_Shared_Win,ElemSideNodeID_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemSideNodeID_Shared_Win,IERROR)

  MPISharedSize = INT(3*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemMidPoint_Shared_Win,ElemMidPoint_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemMidPoint_Shared_Win,IERROR)
#else
  ALLOCATE(ConcaveElemSide_Shared(   1:6,1:nElems))
  ALLOCATE(ElemSideNodeID_Shared(1:4,1:6,1:nElems))
  ALLOCATE(ElemMidPoint_Shared(      1:3,1:nElems))
  firstElem = 1
  lastElem  = nElems
#endif  /*USE_MPI*/

#if USE_MPI
  IF (myComputeNodeRank.EQ.0) THEN
#endif
    ElemSideNodeID_Shared  = 0
    ConcaveElemSide_Shared = .FALSE.
#if USE_MPI
  END IF
  CALL MPI_WIN_SYNC(ElemSideNodeID_Shared_Win,IERROR)
  CALL MPI_WIN_SYNC(ConcaveElemSide_Shared_Win,IERROR)
  CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

  ! iElem is CNElemID
  DO iElem = firstElem,lastElem
    GlobalElemID = GetGlobalElemID(iElem)

    nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
    DO iLocSide = 1,nlocSides
      ! Get global SideID
      GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide

      IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
      localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
      ! Find start of CGNS mapping from flip
      IF (SideInfo_Shared(SIDE_ID,GlobalSideID).GT.0) THEN
        nStart = 0
      ELSE
        nStart = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)-1)
      END IF
      ! Shared memory array starts at 1, but NodeID at 0
      ElemSideNodeID_Shared(1:4,localSideID,iElem) = (/ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart  ,4)+1,localSideID)-1, &
          ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+1,4)+1,localSideID)-1, &
          ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+2,4)+1,localSideID)-1, &
          ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+3,4)+1,localSideID)-1/)

    END DO
  END DO

END ASSOCIATE
#if USE_MPI
CALL MPI_WIN_SYNC(ElemSideNodeID_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

!--- Save whether Side is concave or convex
DO iElem = firstElem,lastElem
  ! iElem is CNElemID
  GlobalElemID = GetGlobalElemID(iElem)

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    !--- Check whether the bilinear side is concave
    !--- Node Number 4 and triangle 1-2-3
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    DO NodeNum = 1,3               ! for all 3 nodes of triangle
       A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,iElem)+1) &
                    - NodeCoords_Shared(:,ElemSideNodeID_Shared(4      ,localSideID,iElem)+1)
    END DO
    !--- concave if detcon < 0:
    detcon = ((A(2,1) * A(3,2) - A(3,1) * A(2,2)) * A(1,3) +     &
              (A(3,1) * A(1,2) - A(1,1) * A(3,2)) * A(2,3) +     &
              (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * A(3,3))
    IF (detcon.LT.0) THEN
      ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
    ELSE IF (detcon.EQ.0.0) THEN
      IF (GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
    END IF
  END DO
END DO

!-- write debug-mesh
! IF (WriteTriaDebugMesh) THEN
!   nSides=6
!   WRITE(UNIT=hilf,FMT='(I4.4)') myRank
!   FileString='TRIA-DebugMesh_PROC'//TRIM(hilf)//'.vtu'
!   ALLOCATE(Coords(1:3,1:4,1:nSides,1:nElems))
!   DO iElem = 1,nElems ; DO iLocSide = 1,nSides ; DO iNode = 1,4
!     Coords(:,iNode,iLocSide,iElem)=GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,iLocSide,iElem))
!   END DO ; END DO ; END DO
!   CALL WriteTriaDataToVTK(nSides,nElems,Coords(1:3,1:4,1:6,1:nElems),FileString)
!   SDEALLOCATE(Coords)
! END IF !WriteTriaDebugMesh

DO iElem = firstElem,lastElem
  ! iElem is CNElemID
  GlobalElemID = GetGlobalElemID(iElem)

  ElemMidPoint_Shared(:,iElem) = 0.
  DO iNode = 1,8
    ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+iNode)
  END DO
  ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) / 8.
END DO

#if USE_MPI
CALL MPI_WIN_SYNC(ConcaveElemSide_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemMidPoint_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

!--- check for elements with intersecting sides (e.g. very flat elements)
CALL WeirdElementCheck()

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GEOMETRY INFORMATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleGeometry


SUBROUTINE InitElemNodeIDs()
!===================================================================================================================================
! Subroutine for particle geometry initialization (GEO container)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_ReadInTools
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,firstElem,lastElem,GlobalElemID
INTEGER            :: iNode
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
INTEGER            :: CornerNodeIDswitch(8)
!===================================================================================================================================

! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
CornerNodeIDswitch(1)=1
CornerNodeIDswitch(2)=(Ngeo+1)
CornerNodeIDswitch(3)=(Ngeo+1)**2
CornerNodeIDswitch(4)=(Ngeo+1)*Ngeo+1
CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2
CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1

! New crazy corner node switch (philipesque)
ASSOCIATE(CNS => CornerNodeIDswitch )
#if USE_MPI
  firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
  lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
  MPISharedSize = INT(8*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/8,nComputeNodeTotalElems/),ElemNodeID_Shared_Win,ElemNodeID_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemNodeID_Shared_Win,IERROR)
#else
  ALLOCATE(ElemNodeID_Shared(1:8,1:nElems))
  firstElem = 1
  lastElem  = nElems
#endif  /*USE_MPI*/

#if USE_MPI
  IF (myComputeNodeRank.EQ.0) THEN
#endif
    ElemNodeID_Shared = 0
#if USE_MPI
  END IF
  CALL MPI_WIN_SYNC(ElemNodeID_Shared_Win,IERROR)
  CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

  ! iElem is CNElemID
  DO iElem = firstElem,lastElem
    GlobalElemID = GetGlobalElemID(iElem)
    DO iNode = 1,8
      ElemNodeID_Shared(iNode,iElem) = ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID) + CNS(iNode)
    END DO
  END DO
END ASSOCIATE
#if USE_MPI
CALL MPI_WIN_SYNC(ElemNodeID_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif
END SUBROUTINE InitElemNodeIDs


!SUBROUTINE WriteTriaDataToVTK(nSides,nElems,Coord,FileString)
!!===================================================================================================================================
!!> Routine writing data to VTK Triangles (cell type = 5)
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT / OUTPUT VARIABLES
!INTEGER,INTENT(IN)          :: nSides               !< Number of sides per element
!INTEGER,INTENT(IN)          :: nElems               !< Number of elements
!REAL   ,INTENT(IN)          :: Coord(1:3,1:4,1:nSides,1:nElems)
!CHARACTER(LEN=*),INTENT(IN) :: FileString           ! < Output file name
!!----------------------------------------------------------------------------------------------------------------------------------!
!! LOCAL VARIABLES
!INTEGER            :: iElem,nVTKElems,nVTKCells,ivtk=44,iSide!,iVal,iVar,str_len
!INTEGER(KIND=8)    :: Offset, nBytes
!INTEGER            :: IntegerType
!INTEGER            :: Vertex(3,nSides*nElems*2)
!INTEGER            :: NodeID,CellID,CellType
!CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
!CHARACTER(LEN=200) :: Buffer
!CHARACTER(LEN=1)   :: lf!,components_string
!!CHARACTER(LEN=255) :: VarNameString
!REAL(KIND=4)       :: FloatType
!!===================================================================================================================================
!SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE TRIA DATA TO VTX XML BINARY (VTU) FILE..."
!IF(nSides.LT.1)THEN
!  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
!  RETURN
!END IF
!
!! Line feed character
!lf = char(10)
!
!! Write file
!OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
!! Write header
!Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
!Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
!nVTKElems=nSides*nElems*4 ! number of Nodes
!nVTKCells=nSides*2*nElems ! number of Triangles
!
!Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
!WRITE(TempStr1,'(I16)')nVTKElems
!WRITE(TempStr2,'(I16)')nVTKCells
!Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//&
!'" NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
!! Specify point data
!Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
!Offset=0
!WRITE(StrOffset,'(I16)')Offset
!!IF (nVal .GT.0)THEN
!!  DO iVar=1,nVal
!!    IF (VarNamePartCombine(iVar).EQ.0) THEN
!!      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNamePartVisu(iVar))//&
!!      '" NumberOfComponents="1" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
!!      Offset=Offset+SIZEOF(IntegerType)+nVTKElems*SIZEOF(FloatType)
!!      WRITE(StrOffset,'(I16)')Offset
!!    ELSE IF (VarNamePartCombine(iVar).EQ.1) THEN
!!      str_len = LEN_TRIM(VarNamePartVisu(iVar))
!!      write(components_string,'(I1)') VarNamePartCombineLen(iVar)
!!      !IF(FileType.EQ.'DSMCHOState')THEN
!!      !  VarNameString = VarNamePartVisu(iVar)(1:str_len-4)//VarNamePartVisu(iVar)(str_len-2:str_len)
!!      !ELSE
!!        VarNameString = VarNamePartVisu(iVar)(1:str_len-1)
!!      !END IF
!!      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameString)//&
!!      '" NumberOfComponents="'//components_string//'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf
!!      WRITE(ivtk) TRIM(Buffer)
!!      Offset=Offset+SIZEOF(IntegerType)+nVTKElems*SIZEOF(FloatType)*VarNamePartCombineLen(iVar)
!!      WRITE(StrOffset,'(I16)')Offset
!!    END IF
!!  END DO
!!END IF
!Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
!! Specify cell data
!Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
!! Specify coordinate data
!Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
!Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"'// &
!       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
!Offset=Offset+SIZEOF(IntegerType)+3*nVTKElems*SIZEOF(FloatType)
!WRITE(StrOffset,'(I16)')Offset
!Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
!! Specify necessary cell data
!Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
!! Connectivity
!Buffer='        <DataArray type="Int32" Name="connectivity" format="appended"'// &
!       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
!Offset=Offset+SIZEOF(IntegerType)+nVTKCells*3*SIZEOF(IntegerType)
!WRITE(StrOffset,'(I16)')Offset
!! Offsets
!Buffer='        <DataArray type="Int32" Name="offsets" format="appended"'// &
!       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
!Offset=Offset+SIZEOF(IntegerType)+nVTKCells*SIZEOF(IntegerType)
!WRITE(StrOffset,'(I16)')Offset
!! Elem types
!Buffer='        <DataArray type="Int32" Name="types" format="appended"'// &
!       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
!Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
!Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
!Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
!! Prepare append section
!Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
!! Write leading data underscore
!Buffer='_';WRITE(ivtk) TRIM(Buffer)
!
!! Write binary raw data into append section
!! Point data
!nBytes = nVTKElems*SIZEOF(FloatType)
!!DO iVal=1,nVal
!!  IF (VarNamePartCombine(iVal).EQ.0) THEN
!!    WRITE(ivtk) nBytes,REAL(Value(1:nParts,iVal),4)
!!  ELSEIF(VarNamePartCombine(iVal).EQ.1) THEN
!!    WRITE(ivtk) nBytes*VarNamePartCombineLen(iVal),REAL(Value(1:nParts,iVal:iVal+VarNamePartCombineLen(iVal)-1),4)
!!  ENDIF
!!END DO
!! Points
!nBytes = nBytes * 3
!WRITE(ivtk) nBytes
!WRITE(ivtk) REAL(Coord,4)
!! Connectivity
!NodeID = -1
!CellID = 1
!DO iElem=1,nElems
!  DO iSide=1,6
!    ! nodes 1,2,3 and nodes 1,3,4 forming one triangle
!    ! nodes indexes start with 0 in vtk
!    Vertex(:,CellID) = (/NodeID+1,NodeID+2,NodeID+3/)
!    Vertex(:,CellID+1) = (/NodeID+1,NodeID+3,NodeID+4/)
!    CellID=CellID+2
!    NodeID=NodeID+4
!  END DO
!END DO
!nBytes = 3*nVTKCells*SIZEOF(IntegerType)
!WRITE(ivtk) nBytes
!WRITE(ivtk) Vertex(:,:)
!! Offset
!nBytes = nVTKCells*SIZEOF(IntegerType)
!WRITE(ivtk) nBytes
!WRITE(ivtk) (Offset,Offset=3,3*nVTKCells,3)
!! Cell type
!CellType = 5  ! VTK_TRIANGLE
!!CellType = 6  ! VTK_TRIANGLE_STRIP
!WRITE(ivtk) nBytes
!WRITE(ivtk) (CellType,iElem=1,nVTKCells)
!! Write footer
!Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
!Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
!CLOSE(ivtk)
!SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
!
!END SUBROUTINE WriteTriaDataToVTK


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
USE MOD_Globals
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared,ConcaveElemSide_Shared,ElemSideNodeID_Shared
USE MOD_Particle_Mesh_Vars        ,ONLY: WeirdElems
USE MOD_Mesh_Tools                ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared_Vars           ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank
#else
USE MOD_Mesh_Vars                 ,ONLY: nElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, iLocSide, kLocSide, iNode
INTEGER,ALLOCATABLE :: WeirdElemNbrs(:)
REAL              :: vec(1:3), Node(1:3,1:4),det(1:3)
LOGICAL           :: WEIRD, TRICHECK, TRIABSCHECK
INTEGER           :: firstElem,lastElem
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' CHECKING FOR WEIRD ELEMENTS...'

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

ALLOCATE(WeirdElemNbrs(1:lastElem-firstElem+1))

WeirdElems = 0

! go through all CN elements
DO iElem = firstElem,lastElem
  WEIRD = .FALSE.
  DO iLocSide = 1,5  ! go through local sides
    IF (.not.WEIRD) THEN  ! if one is found there is no need to continue
      IF (ConcaveElemSide_Shared(iLocSide,iElem)) THEN  ! only concave elements need to be checked
        ! build vector from node 1 to node 3
        vec(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
               - NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,iElem)+1)
        ! check all other sides
        DO kLocSide = iLocSide + 1, 6
          IF (ConcaveElemSide_Shared(kLocSide,iElem)) THEN  ! only concave elements need to be checked
            ! build 4 vectors from point 1 of edge to 4 nodes of kLocSide
            DO iNode = 1,4
              Node(:,iNode) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1    ,iLocSide,iElem)+1) &
                            - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
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
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1    ,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3    ,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
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
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(4,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(4,kLocSide,iElem)+1)
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
    WeirdElemNbrs(WeirdElems) = GetGlobalElemID(iElem)
  END IF
END DO

IF(WeirdElems.GT.0) THEN
  IPWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
  IPWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
  DO iElem = 1,WeirdElems
    IPWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
  END DO
END IF

SWRITE(UNIT_stdOut,'(A)')' CHECKING FOR WEIRD ELEMENTS DONE!'

DEALLOCATE(WeirdElemNbrs)

SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WeirdElementCheck


SUBROUTINE MapRegionToElem()
!----------------------------------------------------------------------------------------------------------------------------------!
! map a particle region to element
! check only element barycenter, nothing else
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars ,ONLY: NbrOfRegions, RegionBounds,GEO
USE MOD_Mesh_Vars          ,ONLY: ElemBaryNGeo
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
 INTEGER                :: iElem, iRegions
!===================================================================================================================================
SDEALLOCATE(GEO%ElemToRegion)
ALLOCATE(GEO%ElemToRegion(1:PP_nElems))
GEO%ElemToRegion=0

DO iElem=1,PP_nElems
  DO iRegions=1,NbrOfRegions
    IF ((ElemBaryNGeo(1,iElem).LT.RegionBounds(1,iRegions)).OR.(ElemBaryNGEO(1,iElem).GE.RegionBounds(2,iRegions))) CYCLE
    IF ((ElemBaryNGeo(2,iElem).LT.RegionBounds(3,iRegions)).OR.(ElemBaryNGEO(2,iElem).GE.RegionBounds(4,iRegions))) CYCLE
    IF ((ElemBaryNGeo(3,iElem).LT.RegionBounds(5,iRegions)).OR.(ElemBaryNGEO(3,iElem).GE.RegionBounds(6,iRegions))) CYCLE
    IF (GEO%ElemToRegion(iElem).EQ.0) THEN
      GEO%ElemToRegion(iElem)=iRegions
    ELSE
      CALL ABORT(__STAMP__,'Defined regions are overlapping')
    END IF
  END DO ! iRegions=1,NbrOfRegions
END DO ! iElem=1,PP_nElems


END SUBROUTINE MapRegionToElem


SUBROUTINE BuildElementBasisAndRadius()
!================================================================================================================================
! Build the element local basis system, where the origin is located at xi=(0,0,0)^T and each local coordinate system is pointing
! to an element side
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis                  ,ONLY: DeCasteljauInterpolation
USE MOD_Basis                  ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars              ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
#if USE_MPI
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadiusNGEO_Shared,ElemRadiusNGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo_Shared,ElemRadius2NGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: nELems
USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,SideID
INTEGER                        :: i,j,k,ilocSide
INTEGER                        :: iDir
REAL                           :: Xi(3,6),xPos(3),Radius
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!================================================================================================================================
#if USE_MPI
  MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadiusNGeo_Shared_Win,ElemRadiusNGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadiusNGeo_Shared_Win,IERROR)
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadius2NGeo_Shared_Win,ElemRadius2NGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadius2NGeo_Shared_Win,IERROR)
  MPISharedSize = INT((3*6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,6,nComputeNodeTotalElems/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
  CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
  MPISharedSize = INT((6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
  CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
  ElemRadiusNGeo     => ElemRadiusNGeo_Shared
  ElemRadius2NGeo    => ElemRadius2NGeo_Shared
  XiEtaZetaBasis     => XiEtaZetaBasis_Shared
  slenXiEtaZetaBasis => slenXiEtaZetaBasis_Shared

ASSOCIATE(XCL_NGeo     => XCL_NGeo_Shared)

#else
  ALLOCATE(ElemRadiusNGeo(          nElems) &
          ,ElemRadius2NGeo(         nElems) &
          ,XiEtaZetaBasis(1:3,1:6,1:nElems) &
          ,slenXiEtaZetaBasis(1:6,1:nElems))
#endif /*USE_MPI*/

#if USE_MPI
IF(myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ElemRadiusNGeo =0.
  ElemRadius2NGeo=0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

#if USE_MPI
firstElem=INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif /*USE_MPI*/

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

! iElem is CN elem
DO iElem=firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos=0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos=xPos+XCL_NGeo(:,i,j,k,ElemID)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    XiEtaZetaBasis(1:3,iDir,iElem)=xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(:,iDir,iElem)=XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
  END DO ! iDir = 1, 6

  Radius=0.
  DO ilocSide=1,6
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)
    IF(SideID.EQ.-1) CYCLE
    DO j=0,NGeo
      DO i=0,NGeo
        xPos=BezierControlPoints3D(:,i,j,SideID)-ElemBaryNGeo(:,iElem)
        Radius=MAX(Radius,SQRT(DOT_PRODUCT(xPos,xPos)))
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO ! ilocSide
  ElemRadiusNGeo (iElem)=Radius
  ElemRadius2NGeo(iElem)=Radius*Radius
END DO ! iElem

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(XiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(slenXiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementBasisAndRadius


SUBROUTINE BuildElementRadiusTria()
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo,ElemRadius2NGeo, ElemRadiusNGeo
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#if USE_MPI
USE MOD_MPI_Shared!            ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo_Shared,ElemRadius2NGeo_Shared,ElemRadiusNGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared_Win,ElemRadius2NGeo_Shared_Win, ElemRadiusNGeo_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: nELems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,iNode
REAL                           :: xPos(3),Radius
INTEGER                        :: firstElem, lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!================================================================================================================================

#if USE_MPI
  MPISharedSize = INT((3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)
  MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadius2NGeo_Shared_Win,ElemRadius2NGEO_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadius2NGeo_Shared_Win,IERROR)
  ElemRadius2NGeo    => ElemRadius2NGeo_Shared
  ElemBaryNGeo       => ElemBaryNGeo_Shared
  IF(StringBeginsWith(DepositionType,'shape_function'))THEN
    CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadiusNGeo_Shared_Win,ElemRadiusNGEO_Shared)
    CALL MPI_WIN_LOCK_ALL(0,ElemRadiusNGeo_Shared_Win,IERROR)
    ElemRadiusNGeo    => ElemRadiusNGeo_Shared
  END IF

#else
ALLOCATE(ElemBaryNGeo(1:3,nElems) &
        ,ElemRadius2NGeo( nElems))
IF(StringBeginsWith(DepositionType,'shape_function')) ALLOCATE(ElemRadiusNGeo(nElems))
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemRadius2NGeo = 0.
  IF(StringBeginsWith(DepositionType,'shape_function')) ElemRadiusNGeo = 0.
#if USE_MPI
END IF
IF(StringBeginsWith(DepositionType,'shape_function')) CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

#if USE_MPI
firstElem=INT(REAL(myComputeNodeRank*    nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif /*USE_MPI*/

DO iElem=firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  Radius=0.
  xPos  =0.

  DO iNode=1,8
    xPos = xPos + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode)
  END DO
    ElemBaryNGeo(:,iElem) = xPos/8.
  DO iNode=1,8
    xPos   = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode) - ElemBaryNGeo(:,iElem)
    Radius = MAX(Radius,VECNORM(xPos))
  END DO
  ElemRadius2NGeo(iElem) = Radius*Radius
  IF(StringBeginsWith(DepositionType,'shape_function')) ElemRadiusNGeo(iElem) = Radius
END DO ! iElem

#if USE_MPI
IF(StringBeginsWith(DepositionType,'shape_function')) CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementRadiusTria


SUBROUTINE BuildElemTypeAndBasisTria()
!===================================================================================================================================
!> Dummy routine to fill the ElemCurved array with TriaTracking
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis                   ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Mesh_Vars               ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
USE MOD_MPI_Shared!             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars               ,ONLY: nElems
USE MOD_Mesh_Vars               ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,iDir
INTEGER                        :: i,j,k
REAL                           :: xPos(3)
REAL                           :: Xi(3,6),Lag(1:3,0:NGeo)
INTEGER                        :: firstElem, lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /* USE_MPI */
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved ...'

! elements
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
MPISharedSize = INT((3*6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,6,nComputeNodeTotalElems/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
MPISharedSize = INT((6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
ElemCurved         => ElemCurved_Shared
XiEtaZetaBasis     => XiEtaZetaBasis_Shared
slenXiEtaZetaBasis => slenXiEtaZetaBasis_Shared

ASSOCIATE(XCL_NGeo     => XCL_NGeo_Shared)

#else
ALLOCATE(ElemCurved(            1:nElems) &
        ,XiEtaZetaBasis(1:3,1:6,1:nElems) &
        ,slenXiEtaZetaBasis(1:6,1:nElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ElemCurved   = .FALSE.
#if USE_MPI
END IF
#endif /*USE_MPI*/

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos = 0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos = xPos+XCL_NGeo(:,i,j,k,ElemID)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    XiEtaZetaBasis(1:3,iDir,iElem) = xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(:,iDir,iElem)   = XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem) = 1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
  END DO ! iDir = 1, 6
END DO

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemCurved_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(XiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(slenXiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElemTypeAndBasisTria


SUBROUTINE PointsEqual(N,Points1,Points2,IsNotEqual)
!===================================================================================================================================
! compute the distance between two data sets
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: N
REAL,INTENT(IN)           :: Points1(1:3,1:N)
REAL,INTENT(IN)           :: Points2(1:3,1:N)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL                   :: IsNotEqual
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i
!===================================================================================================================================

IsNotEqual=.FALSE.

DO i=1,N
  IF( ABS(Points1(1,i)-Points2(1,i)).GT.1e-14 .OR. &
      ABS(Points1(2,i)-Points2(2,i)).GT.1e-14 .OR. &
      ABS(Points1(3,i)-Points2(3,i)).GT.1e-14 ) THEN
    IsNotEqual=.TRUE.
    RETURN
  END IF
END DO ! i=0,N

END SUBROUTINE PointsEqual


SUBROUTINE BuildElementOriginShared()
!================================================================================================================================
! compute the element origin at xi=(0,0,0)^T and set it as ElemBaryNGeo
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars          ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars ,ONLY: ElemBaryNGeo
#if USE_MPI
USE MOD_MPI_Shared!        ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars, ONLY: ElemBaryNGeo_Shared,ElemBaryNGeo_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Mesh_Vars ,ONLY: nComputeNodeElems
USE MOD_Mesh_Vars          ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,i,j,k
REAL                           :: XPos(3),buf
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!================================================================================================================================
#if USE_MPI
MPISharedSize = INT((3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)
ElemBaryNGeo => ElemBaryNGeo_Shared

ASSOCIATE(XCL_NGeo => XCL_NGeo_Shared)

! Set ranges
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(ElemBaryNGeo(1:3,nComputeNodeElems))
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemBaryNGeo = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

! evaluate the polynomial at origin: Xi=(/0.0,0.0,0.0/)
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  xPos=0.

  DO k=0,NGeo
    DO j=0,NGeo
      buf=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        xPos=xPos+XCL_NGeo(1:3,i,j,k,ElemID)*Lag(1,i)*buf
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo

  ElemBaryNGeo(:,iElem)=xPos
END DO ! iElem

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementOriginShared


SUBROUTINE IdentifyElemAndSideType()
!===================================================================================================================================
!> get the element and side type of each element
!> 1) Get Elem Type (curved_elem)
!> 2) Get Side Type (planar_rect, planar_nonrect, bilineard, curved, planar_curved)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_ChangeBasis            ,ONLY: changeBasis3D
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides
USE MOD_Mesh_Vars              ,ONLY: Vdm_CLNGeo1_CLNGeo,NGeo,Vdm_CLNGeo1_CLNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared,ElemBaryNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared,ElemCurved
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,SideType,SideNormVec,SideDistance
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideDistance_Shared,SideDistance_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideType_Shared,SideType_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideNormVec_Shared,SideNormVec_Shared_Win
#else
USE MOD_Particle_Mesh_Vars     ,ONLY: nComputeNodeElems
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: ilocSide,SideID,flip
INTEGER                                  :: iElem,firstElem,lastElem,ElemID
REAL,DIMENSION(1:3)                      :: v1,v2,v3
LOGICAL,ALLOCATABLE                      :: SideIsDone(:)
REAL                                     :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                                     :: XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                                     :: BezierControlPoints_loc(1:3,0:NGeo,0:NGeo)
INTEGER                                  :: NGeo3,NGeo2
REAL                                     :: XCL_NGeoSideNew(1:3,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoSideOld(1:3,0:NGeo,0:NGeo)
LOGICAL                                  :: isCurvedSide,isRectangular
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND)           :: MPISharedSize
#endif /* USE_MPI */
! output and sanity check
INTEGER                                  :: nPlanarRectangular,   nPlanarNonRectangular,   nPlanarCurved,   nBilinear,   nCurved
INTEGER                                  :: nPlanarRectangularTot,nPlanarNonRectangularTot,nPlanarCurvedTot,nBilinearTot,nCurvedTot
INTEGER                                  :: nLinearElems,   nCurvedElems
INTEGER                                  :: nLinearElemsTot,nCurvedElemsTot
!#if USE_MPI
!INTEGER                                  :: nDummy
!#endif /* USE_MPI */
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved ...'

! elements
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
ElemCurved => ElemCurved_Shared
#else
ALLOCATE(ElemCurved(1:nComputeNodeElems))
#endif /*USE_MPI*/

! sides
#if USE_MPI
MPISharedSize = INT((nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),SideType_Shared_Win,SideType_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideType_Shared_Win,IERROR)
SideType => SideType_Shared
CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),SideDistance_Shared_Win,SideDistance_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideDistance_Shared_Win,IERROR)
SideDistance => SideDistance_Shared
MPISharedSize = INT((3*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),SideNormVec_Shared_Win,SideNormVec_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideNormVec_Shared_Win,IERROR)
SideNormVec => SideNormVec_Shared
#else
ALLOCATE(SideType(       nNonUniqueGlobalSides))
ALLOCATE(SideDistance(   nNonUniqueGlobalSides))
ALLOCATE(SideNormVec(1:3,nNonUniqueGlobalSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemCurved   = .FALSE.
  SideType     = -1
  SideDistance = -0.
  SideNormVec  = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemCurved_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideType_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideDistance_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideNormVec_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

ALLOCATE(SideIsDone(nNonUniqueGlobalSides))
SideIsDone = .FALSE.

NGeo2 = (NGeo+1)*(NGeo+1)
NGeo3 = NGeo2   *(NGeo+1)

! decide if element is (bi-)linear or curved
! decide if sides are planar-rect, planar-nonrect, planar-curved, bilinear or curved
#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

DO iElem=firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)

  XCL_NGeoLoc = XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID)
  ! 1) check if elem is curved
  !   a) get the coordinates of the eight nodes of the hexahedral
  XCL_NGeo1(1:3,0,0,0) = XCL_NGeoLoc(1:3, 0  , 0  , 0  )
  XCL_NGeo1(1:3,1,0,0) = XCL_NGeoLoc(1:3,NGeo, 0  , 0  )
  XCL_NGeo1(1:3,0,1,0) = XCL_NGeoLoc(1:3, 0  ,NGeo, 0  )
  XCL_NGeo1(1:3,1,1,0) = XCL_NGeoLoc(1:3,NGeo,NGeo, 0  )
  XCL_NGeo1(1:3,0,0,1) = XCL_NGeoLoc(1:3, 0  , 0  ,NGeo)
  XCL_NGeo1(1:3,1,0,1) = XCL_NGeoLoc(1:3,NGeo, 0  ,NGeo)
  XCL_NGeo1(1:3,0,1,1) = XCL_NGeoLoc(1:3, 0  ,NGeo,NGeo)
  XCL_NGeo1(1:3,1,1,1) = XCL_NGeoLoc(1:3,NGeo,NGeo,NGeo)

  !  b) interpolate from the nodes to NGeo
  !     Compare the bi-linear mapping with the used mapping
  !     For NGeo=1, this should always be true, because the mappings are identical
  CALL ChangeBasis3D(3,1,NGeo,Vdm_CLNGeo1_CLNGeo,XCL_NGeo1,XCL_NGeoNew)
  ! check the coordinates of all Chebychev-Lobatto geometry points between the bi-linear and used mapping
  CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo),ElemCurved(iElem))

  ! 2) check sides
  ! loop over all 6 sides of element
  ! a) check if the sides are straight
  ! b) use curved information to decide side type
  DO ilocSide=1,6
    SideID = GetGlobalNonUniqueSideID(ElemID,iLocSide)
    flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

    IF(.NOT.ElemCurved(iElem))THEN
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      ! linear element
      IF(BoundingBoxIsEmpty(SideID))THEN
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,SideID) = CROSSNORM(v1,v2)
        v1=0.25*(BezierControlPoints_loc(:,0,0      )     &
                +BezierControlPoints_loc(:,NGeo,0   )  &
                +BezierControlPoints_loc(:,0,NGeo   )  &
                +BezierControlPoints_loc(:,NGeo,NGeo))
        ! check if normal vector points outwards
        v2=v1-ElemBaryNGeo(:,iElem)
        IF(flip.EQ.0)THEN
          IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).LT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
        ELSE
          IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).GT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
        END IF
        SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
        ! check if it is rectangular
        isRectangular=.TRUE.
        v1=UNITVECTOR(BezierControlPoints_loc(:,0   ,NGeo)-BezierControlPoints_loc(:,0   ,0   ))
        v2=UNITVECTOR(BezierControlPoints_loc(:,NGeo,0   )-BezierControlPoints_loc(:,0   ,0   ))
        v3=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,0   ,NGeo))
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        IF(isRectangular)THEN
          v1=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,NGeo,0))
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        END IF
        IF(isRectangular)THEN
          SideType(SideID)=PLANAR_RECT
        ELSE
          SideType(SideID)=PLANAR_NONRECT
        END IF
      ELSE
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,SideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
        SideType(SideID)=BILINEAR
      END IF
    ELSE
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      ! possible curved face
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0,0:NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0,0:NGeo,0:NGeo)
      CASE(XI_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,NGeo,0:NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,NGeo,0:NGeo,0:NGeo)
      CASE(ETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0,0:NGeo)
      CASE(ETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,NGeo,0:NGeo)
      CASE(ZETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0)
      CASE(ZETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,NGeo)
        XCL_NGeoSideNew=XCL_NGeoNEw(1:3,0:NGeo,0:NGeo,NGeo)
      END SELECT
      CALL PointsEqual(NGeo2,XCL_NGeoSideNew,XCL_NGeoSideOld,isCurvedSide)
      IF(isCurvedSide)THEn
        IF(BoundingBoxIsEmpty(SideID))THEN
          SideType(SideID)=PLANAR_CURVED
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,SideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0   )  &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).LT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).GT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          END IF
          SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
        ELSE
          SideType(SideID)=CURVED
        END IF
      ELSE
        IF (BoundingBoxIsEmpty(SideID)) THEN
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,SideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0)     &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).LT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).GT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          END IF
          SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
          ! check if it is rectangular
          isRectangular=.TRUE.
          v1=UNITVECTOR(BezierControlPoints_loc(:,0   ,NGeo)-BezierControlPoints_loc(:,0   ,0   ))
          v2=UNITVECTOR(BezierControlPoints_loc(:,NGeo,0   )-BezierControlPoints_loc(:,0   ,0   ))
          v3=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,0   ,NGeo))
!          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
!          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
          IF(DOT_PRODUCT(v1,v2).GT.1E-14) isRectangular=.FALSE.
          IF(DOT_PRODUCT(v1,v3).GT.1E-14) isRectangular=.FALSE.
          IF(isRectangular)THEN
            v1=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,NGeo,0))
!            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
!            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
            IF(DOT_PRODUCT(v1,v2).GT.1E-14) isRectangular=.FALSE.
            IF(DOT_PRODUCT(v1,v3).GT.1E-14) isRectangular=.FALSE.
          END IF
          IF(isRectangular)THEN
            SideType(SideID)=PLANAR_RECT
          ELSE
            SideType(SideID)=PLANAR_NONRECT
          END IF
        ELSE
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          SideNormVec(:,SideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
          SideType(SideID)=BILINEAR
        END IF
      END IF
    END IF
    SideIsDone(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iElem=1,nTotalElems

#if USE_MPI
CALL MPI_WIN_SYNC(ElemCurved_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */

! zero counter for side types
nPlanarRectangular         = 0
nPlanarNonRectangular      = 0
nPlanarCurved              = 0
nBilinear                  = 0
nCurved                    = 0
! zero counter for elem types
nCurvedElems               = 0
nLinearElems               = 0

#if USE_MPI
firstElem = offsetElem + 1
lastElem  = offsetElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif


DO iElem = firstElem,lastElem
  ElemID = GetCNElemID(iElem)
  IF (ElemCurved(ElemID)) THEN
    nCurvedElems = nCurvedElems+1
  ELSE
    nLinearElems = nLinearElems+1
  END IF

  DO ilocSide = 1,6
    ! ignore small mortar sides attached to big mortar sides
    SideID = GetGlobalNonUniqueSideID(iElem,ilocSide)
    SELECT CASE(SideType(SideID))
      CASE (PLANAR_RECT)
        nPlanarRectangular    = nPlanarRectangular   +1
      CASE (PLANAR_NONRECT)
        nPlanarNonRectangular = nPlanarNonRectangular+1
      CASE (BILINEAR)
        nBilinear             = nBilinear            +1
      CASE (PLANAR_CURVED)
        nPlanarCurved         = nPlanarCurved        +1
      CASE (CURVED)
        nCurved               = nCurved              +1
    END SELECT
  END DO
END DO

#if USE_MPI
CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
CALL MPI_REDUCE(nLinearElems         ,nLinearElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
CALL MPI_REDUCE(nCurvedElems         ,nCurvedElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
#else
nPlanarRectangularTot    = nPlanarRectangular
nPlanarNonRectangularTot = nPlanarNonRectangular
nBilinearTot             = nBilinear
nPlanarCurvedTot         = nPlanarCurved
nCurvedTot               = nCurved
nLinearElemsTot          = nLinearElems
nCurvedElemsTot          = nCurvedElems
#endif /* USE_MPI */

! sanity check
! This only works if the full mesh is built on one node!
!#if USE_MPI
!IF (myComputeNodeRank.EQ.0) THEN
!#endif /* USE_MPI */
!  IF ((nComputeNodeTotalElems-nCurvedElemsTot).NE.nLinearElemsTot) &
!    CALL ABORT(__STAMP__, 'Error in particle mesh: lost elements while trying to dermine if elements are curved')
!#if USE_MPI
!END IF
!#endif /* USE_MPI */

SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-rectangular     faces: ', nPlanarRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-non-rectangular faces: ', nPlanarNonRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of bi-linear              faces: ', nBilineartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-curved          faces: ', nPlanarCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 faces: ', nCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of (bi-)linear            elems: ', nLinearElemsTot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 elems: ', nCurvedElemsTot

END SUBROUTINE IdentifyElemAndSideType


SUBROUTINE GetBCSidesAndOrgin()
!===================================================================================================================================
! Globally identifies all BC sides and build side origin and radius
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                  ,ONLY: DeCasteljauInterpolation
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueBCSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalSides
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Mesh_Tools             ,ONLY: GetGlobalSideID
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID_Shared,SideID2BCSide_Shared,BCSideMetrics_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID_Shared_Win,SideID2BCSide_Shared_Win,BCSideMetrics_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides
#else
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: p,q
INTEGER                        :: iSide,SideID,firstSide,lastSide,BCSideID
INTEGER                        :: nUniqueBCSidesProc,offsetUniqueBCSidesProc
REAL                           :: origin(1:3),xi(1:2),radius,radiusMax,vec(1:3)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! Count number of BC sides in range
nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
  SideID = GetGlobalSideID(iSide)

  ! ignore elements not on the compute node
  IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
END DO

! Find global number of BC sides and write side <=> BCSide mapping into shared array
#if USE_MPI
sendbuf = nUniqueBCSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetUniqueBCSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetUniqueBCSidesProc + nUniqueBCSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nUniqueBCSides = sendbuf

MPISharedSize = INT((nUniqueBCSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nUniqueBCSides/),BCSide2SideID_Shared_Win,BCSide2SideID_Shared)
CALL MPI_WIN_LOCK_ALL(0,BCSide2SideID_Shared_Win,IERROR)
MPISharedSize = INT((nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),SideID2BCSide_Shared_Win,SideID2BCSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideID2BCSide_Shared_Win,IERROR)
BCSide2SideID => BCSide2SideID_Shared
SideID2BCSide => SideID2BCSide_Shared

! Also allocate array to hold BC Side metrics
MPISharedSize = INT((4*nUniqueBCSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/4,nUniqueBCSides/),BCSideMetrics_Shared_Win,BCSideMetrics_Shared)
CALL MPI_WIN_LOCK_ALL(0,BCSideMetrics_Shared_Win,IERROR)
BCSideMetrics => BCSideMetrics_Shared
#else
offsetUniqueBCSidesProc = 0
nUniqueBCSides = nUniqueBCSidesProc

ALLOCATE(BCSide2SideID(    1:nUniqueBCSides)        &
        ,SideID2BCSide(    1:nNonUniqueGlobalSides))
! Also allocate array to hold BC Side metrics
ALLOCATE(BCSideMetrics(1:4,1:nUniqueBCSides))
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  BCSide2SideID = -1
  SideID2BCSide = -1
  BCSideMetrics = -1
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(BCSide2SideID_Shared_Win,iError)
CALL MPI_WIN_SYNC(SideID2BCSide_Shared_Win,iError)
CALL MPI_WIN_SYNC(BCSideMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
  ! ignore elements not on the compute node
  IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
  BCSideID           = offsetUniqueBCSidesProc + nUniqueBCSidesProc
  BCSide2SideID(BCSideID) = iSide
  SideID2BCSide(iSide)    = BCSideID

  ! calculate origin, radius for all BC sides
  !> build side origin
  xi     = 0.
  ! TODO: BezierControlPoints are allocated with global side ID, so this SHOULD work. Breaks if we reduce the halo region
  CALL DeCasteljauInterpolation(NGeo,xi,iSide,origin)
  BCSideMetrics(1:3,BCSideID) = origin(1:3)

  !> build side radius
  radiusMax = 0.
  DO q = 0,NGeo
    DO p = 0,NGeo
      vec(1:3) = BezierControlPoints3D(:,p,q,iSide) - origin
      radius   = DOT_PRODUCT(Vec,Vec)
      radiusMax= MAX(radiusMax,radius)
    END DO
  END DO
  BCSideMetrics(4,BCSideID) = SQRT(RadiusMax)
END DO

#if USE_MPI
CALL MPI_WIN_SYNC(BCSide2SideID_Shared_Win,iError)
CALL MPI_WIN_SYNC(SideID2BCSide_Shared_Win,iError)
CALL MPI_WIN_SYNC(BCSideMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE GetBCSidesAndOrgin


SUBROUTINE BuildBCElemDistance()
!===================================================================================================================================
! get the distance of each BC face
!> 1) identify BC elements to be handled by Tracing
!> 2) build mapping global elem ID to BC elem ID
!> 3) calculate distance from element origin to each side
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: c
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides,SideBCMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo,ElemRadiusNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,nUniqueBCSides
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Vars          ,ONLY: ManualTimeStep
USE MOD_Utils                  ,ONLY: InsertionSort
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides_Shared,SideBCMetrics_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides_Shared_Win,SideBCMetrics_Shared_Win
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps,halo_eps_velo
#if ! (USE_HDG)
USE MOD_CalcTimeStep           ,ONLY: CalcTimeStep
#endif
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars          ,ONLY: nRKStages,RK_c
#endif
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: nComputeNodeElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: p,q,r,s
INTEGER                        :: firstBezierPoint,lastBezierPoint
INTEGER                        :: ElemID,SideID
INTEGER                        :: iElem,firstElem,lastElem
INTEGER                        :: iSide,firstSide,lastSide,iLocSide
INTEGER                        :: nComputeNodeBCSides
INTEGER                        :: nBCSidesElem,nBCSidesProc,offsetBCSidesProc,offsetBCSides
INTEGER                        :: iBCSide,BCElemID,BCSideID
INTEGER                        :: CNElemID,BCCNElemID
INTEGER,ALLOCATABLE            :: ElemToBCSidesProc(:,:)
REAL                           :: dX,dY,dZ
REAL                           :: origin(1:3),vec(1:3)
REAL                           :: BC_halo_eps
LOGICAL                        :: fullMesh
REAL,ALLOCATABLE               :: tmpSideBCMetrics(:,:)
REAL,ALLOCATABLE               :: tmpSideBCDistance(:)
INTEGER,ALLOCATABLE            :: intSideBCMetrics(:)
#if USE_MPI
REAL                           :: BC_halo_eps_velo,BC_halo_diag,deltaT
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
INTEGER                        :: iStage
#endif
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying BC sides and calculating side metrics ...'

! elements
#if USE_MPI
MPISharedSize = INT((2*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,nComputeNodeTotalElems/),ElemToBCSides_Shared_Win,ElemToBCSides_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBCSides_Shared_Win,IERROR)
ElemToBCSides => ElemToBCSides_Shared
#else
ALLOCATE(ElemToBCSides(1:2,1:nComputeNodeElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemToBCSides = -1
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemToBCSides_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
firstSide = 1
lastSide  = nNonUniqueGlobalSides
ALLOCATE(ElemToBCSidesProc(1:2,1:nComputeNodeTotalElems))

! if running on one node, halo_eps is meaningless. Get a representative BC_halo_eps for BC side identification
fullMesh = .FALSE.
IF (halo_eps.EQ.0) THEN
  ! reconstruct halo_eps_velo
  IF (halo_eps_velo.EQ.0) THEN
    BC_halo_eps_velo = c
  ELSE
    BC_halo_eps_velo = halo_eps_velo
  END IF

  ! reconstruct deltaT
  deltaT = 0.
  IF (ManualTimeStep.GT.0.) THEN
    deltaT    = ManualTimeStep
#if ! (USE_HDG)
  ELSE
    deltaT    = CalcTimeStep()
#endif
  END IF

  ! calculate halo_eps
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  BC_halo_eps = RK_c(2)
  DO iStage=2,nRKStages-1
    BC_halo_eps = MAX(BC_halo_eps,RK_c(iStage+1)-RK_c(iStage))
  END DO
  BC_halo_eps = MAX(BC_halo_eps,1.-RK_c(nRKStages))
  BC_halo_eps = BC_halo_eps*BC_halo_eps_velo*deltaT
#else
  BC_halo_eps = BC_halo_eps_velo*deltaT
#endif

  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  BC_halo_diag = VECNORM(vec)

  ! compare halo_eps against global diagonal and reduce if necessary
  IF (.NOT.ALMOSTZERO(BC_halo_eps).AND.(BC_halo_diag.GE.BC_halo_eps)) THEN
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to ',BC_halo_eps
  ELSEIF (.NOT.ALMOSTZERO(BC_halo_eps).AND.(BC_halo_diag.LT.BC_halo_eps)) THEN
    fullMesh = .TRUE.
    BC_halo_eps = BC_halo_diag
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to global diag with ',BC_halo_eps
  ! halo_eps still at zero. Set it to global diagonal
  ELSE
    fullMesh = .TRUE.
    BC_halo_eps = BC_halo_diag
    SWRITE(UNIT_stdOUt,'(A,F11.3)') ' | No halo_eps given and could not be reconstructed. Using global diag with ',BC_halo_eps
  END IF
ELSE
  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  BC_halo_diag = VECNORM(vec)

  IF (BC_halo_diag.LE.halo_eps) fullMesh = .TRUE.
  BC_halo_eps = halo_eps
END IF

#else
! get distance of diagonal of mesh
vec(1)   = GEO%xmaxglob-GEO%xminglob
vec(2)   = GEO%ymaxglob-GEO%yminglob
vec(3)   = GEO%zmaxglob-GEO%zminglob
BC_halo_eps = DOT_PRODUCT(vec,vec)
fullMesh = .TRUE.

firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nNonUniqueGlobalSides
ALLOCATE(ElemToBCSidesProc(1:2,1:nElems))
#endif /*USE_MPI*/

ElemToBCSidesProc = -1
nBCSidesProc      = 0
offsetBCSides     = 0

! for fullMesh, each element requires ALL BC faces
IF (fullMesh) THEN
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

    DO iBCSide = 1,nUniqueBCSides
      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO ! iBCSide

    ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
    IF (nBCSidesElem.GT.0) THEN
      ElemToBCSidesProc(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
      ElemToBCSidesProc(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides
    END IF
  END DO ! iElem

  offsetBCSides = nBCSidesProc

! .NOT. fullMesh
ELSE
  ! sum up all BC sides in range of BC_halo_eps
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)
      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

    ! loop over all sides. Check distance from every local side to total sides.
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)
      BCCNElemID = GetCNElemID(BCElemID)

      ! Ignore elements not on the compute node
      IF (BCCNElemID.EQ.-1) CYCLE

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,iElem) - ElemBaryNGeo(:,BCCNElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(iElem) + ElemRadiusNGeo(BCCNElemID))) CYCLE

      ! loop over all local sides of the element. Use a named loop so the entire element can be cycled
Check1: DO ilocSide = 1,6
      SideID = GetGlobalNonUniqueSideID(ElemID,ilocSide)

        ! compare all nodes between local and BC side to check if within BC_halo_eps. Once one node pair is in range, flag the entire
        ! side and stop checking. First, get BezierControlPoints for local side. BezierControlPoints3D available for ALL sides in shared memory
        SELECT CASE(ilocSide)
          CASE(XI_MINUS,XI_PLUS)
            firstBezierPoint = 0
            lastBezierPoint  = NGeo
          CASE DEFAULT
            firstBezierPoint = 1
            lastBezierPoint  = NGeo-1
        END SELECT

        ! finally compare the node coords
        DO q = firstBezierPoint,lastBezierPoint
          DO p = firstBezierPoint,lastBezierPoint
  !           ! get all nodes for BC side
  !           NodeBCSide(:) = BezierControlPoints3D(:,p,q,BCSideID)
            ! finally compare the node coords
            DO s = firstBezierPoint,lastBezierPoint
              DO r = firstBezierPoint,lastBezierPoint
                dX = ABS(BezierControlPoints3D(1,r,s,SideID)-BezierControlPoints3D(1,p,q,BCSideID))
                IF (dX.GT.BC_halo_eps) CYCLE
                dY = ABS(BezierControlPoints3D(2,r,s,SideID)-BezierControlPoints3D(2,p,q,BCSideID))
                IF (dY.GT.BC_halo_eps) CYCLE
                dZ = ABS(BezierControlPoints3D(3,r,s,SideID)-BezierControlPoints3D(3,p,q,BCSideID))
                IF (dZ.GT.BC_halo_eps) CYCLE

                IF (SQRT(dX*dX+dY*dY+dZ*dZ).LE.BC_halo_eps) THEN
                  nBCSidesElem = nBCSidesElem + 1
                  nBCSidesProc = nBCSidesProc + 1
                  EXIT Check1
                END IF
              END DO ! r
            END DO ! s
          END DO ! p
        END DO ! q
      END DO Check1 ! ilocSide
    END DO ! iBCSide

    ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
    IF (nBCSidesElem.GT.0) THEN
      ElemToBCSidesProc(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
      ElemToBCSidesProc(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides
    END IF

    offsetBCSides = nBCSidesProc
  END DO ! iElem
END IF ! fullMesh

! Find CN global number of BC sides and write Elem to BC Side mapping into shared array
#if USE_MPI
sendbuf = nBCSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetBCSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetBCSidesProc + nBCSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeBCSides = sendbuf

ElemToBCSides(ELEM_NBR_BCSIDES ,firstElem:lastElem) = ElemToBCSidesProc(ELEM_NBR_BCSIDES ,firstElem:lastElem)
ElemToBCSides(ELEM_FIRST_BCSIDE,firstElem:lastElem) = ElemToBCSidesProc(ELEM_FIRST_BCSIDE,firstElem:lastElem) + offsetBCSidesProc
#else
offsetBCSidesProc   = 0
nComputeNodeBCSides = nBCSidesProc

ElemToBCSides(:,firstElem:lastElem) = ElemToBCSidesProc(:,firstElem:lastElem)
#endif /*USE_MPI*/
DEALLOCATE(ElemToBCSidesProc)

! Allocate shared array for BC sides
#if USE_MPI
MPISharedSize = INT((7*nComputeNodeBCSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/7,nComputeNodeBCSides/),SideBCMetrics_Shared_Win,SideBCMetrics_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideBCMetrics_Shared_Win,IERROR)
SideBCMetrics => SideBCMetrics_Shared
#else
ALLOCATE(SideBCMetrics(1:7,1:nComputeNodeBCSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  SideBCMetrics = -1.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemToBCSides_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

nBCSidesProc      = 0

! We did not know the number of BC sides before. Therefore, we need to do the check again and build the final mapping
! for fullMesh, each element requires ALL BC faces
IF (fullMesh) THEN
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(iElem)
    END DO

    DO iBCSide = 1,nUniqueBCSides
      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(BCElemID)
    END DO ! iBCSide
  END DO ! iElem

! .NOT. fullMesh
ELSE
  ! sum up all BC sides in range of BC_halo_eps
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(ElemID)
    END DO

    ! loop over all sides. Check distance from every local side to total sides. Once a side has been flagged,
    ! it must not be counted again
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)
      BCCNElemID = GetCNElemID(BCElemID)

      ! Ignore elements not on the compute node
      IF (BCCNElemID.EQ.-1) CYCLE

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,iElem) - ElemBaryNGeo(:,BCCNElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(iElem) + ElemRadiusNGeo(BCCNElemID))) CYCLE

      ! loop over all local sides of the element. Use a named loop so the entire element can be cycled
Check2: DO ilocSide = 1,6
        SideID = GetGlobalNonUniqueSideID(ElemID,ilocSide)

        ! compare all nodes between local and BC side to check if within BC_halo_eps. Once one node pair is in range, flag the entire
        ! side and stop checking. First, get BezierControlPoints for local side. BezierControlPoints3D available for ALL sides in shared memory
        SELECT CASE(ilocSide)
          CASE(XI_MINUS,XI_PLUS)
            firstBezierPoint = 0
            lastBezierPoint  = NGeo
          CASE DEFAULT
            firstBezierPoint = 1
            lastBezierPoint  = NGeo-1
        END SELECT

        ! finally compare the node coords
        DO q = firstBezierPoint,lastBezierPoint
          DO p = firstBezierPoint,lastBezierPoint
!           ! get all nodes for BC side
!           NodeBCSide(:) = BezierControlPoints3D(:,p,q,BCSideID)
            ! finally compare the node coords
            DO s = firstBezierPoint,lastBezierPoint
              DO r = firstBezierPoint,lastBezierPoint
                dX = ABS(BezierControlPoints3D(1,r,s,SideID)-BezierControlPoints3D(1,p,q,BCSideID))
                IF (dX.GT.BC_halo_eps) CYCLE
                dY = ABS(BezierControlPoints3D(2,r,s,SideID)-BezierControlPoints3D(2,p,q,BCSideID))
                IF (dY.GT.BC_halo_eps) CYCLE
                dZ = ABS(BezierControlPoints3D(3,r,s,SideID)-BezierControlPoints3D(3,p,q,BCSideID))
                IF (dZ.GT.BC_halo_eps) CYCLE

                IF (SQRT(dX*dX+dY*dY+dZ*dZ).LE.BC_halo_eps) THEN
                  nBCSidesProc = nBCSidesProc + 1
                  SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
                  SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(BCElemID)
                  EXIT Check2
                END IF
              END DO ! r
            END DO ! s
          END DO ! p
        END DO ! q
      END DO Check2 ! ilocSide
    END DO ! iBCSide
  END DO ! iElem
END IF ! fullMesh

#if USE_MPI
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeBCSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeBCSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nComputeNodeBCSides
#endif /*USE_MPI*/

! calculate origin, radius and distance to sides
DO iSide = firstSide,lastSide
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,iSide))
  BCSideID = SideID2BCSide(SideID)
  CNElemID = GetCNElemID(INT(SideBCMetrics(BCSIDE_ELEMID,iSide)))

  !> get origin and radius from BC Side
  SideBCMetrics(5:7          ,iSide) = BCSideMetrics(1:3,BCSideID)
  SideBCMetrics(BCSIDE_RADIUS,iSide) = BCSideMetrics(4  ,BCSideID)

  !> build side distance
  origin(1:3) = ElemBaryNGeo(1:3,CNElemID)
  vec(1:3)    = origin(1:3) - SideBCMetrics(5:7,iSide)
  SideBCMetrics(BCSIDE_DISTANCE,iSide) = SQRT(DOT_PRODUCT(vec,vec))-ElemRadiusNGeo(CNElemID)-SideBCMetrics(BCSIDE_RADIUS,iSide)
END DO ! iSide

#if USE_MPI
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

! finally, sort by distance to help speed up BC tracking
!> allocate dummy array to hold variables
ALLOCATE(tmpSideBCMetrics(1:7,1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))
ALLOCATE(tmpSideBCDistance(   1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))
ALLOCATE(intSideBCMetrics(    1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))

DO iElem = firstElem,lastElem
  ! skip elements with no BC sides
  IF (ElemToBCSides(ELEM_NBR_BCSIDES,iElem).LE.0) CYCLE

  ! save values in temporary array
  firstSide    = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) + 1
  lastSide     = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) + ElemToBCSides(ELEM_NBR_BCSIDES,iElem)
  nBCSidesElem = ElemToBCSides(ELEM_NBR_BCSIDES,iElem)

  tmpSideBCMetrics(:,1:nBCSidesElem) = SideBCMetrics(:,firstSide:lastSide)
  tmpSideBCDistance( 1:nBCSidesElem) = SideBCMetrics(BCSIDE_DISTANCE,firstSide:lastSide)
  intSideBCMetrics(  1:nBCSidesElem) = (/((p),p=1,nBCSidesElem)/)

  ! sort SideID according to distance
  CALL InsertionSort(tmpSideBCDistance(1:nBCSidesElem),intSideBCMetrics(1:nBCSidesElem),nBCSidesElem)

  ! write back dummy array with variables
  DO iSide = 1,nBCSidesElem
    SideID = intSideBCMetrics(iSide)
    SideBCMetrics(:,firstSide+iSide-1) = tmpSideBCMetrics(:,SideID)
  END DO
END DO

DEALLOCATE(tmpSideBCMetrics)
DEALLOCATE(tmpSideBCDistance)
DEALLOCATE(intSideBCMetrics)

#if USE_MPI
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

END SUBROUTINE BuildBCElemDistance


SUBROUTINE BuildEpsOneCell()
!===================================================================================================================================
! Build epsOneCell for each element
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis              ,ONLY: ChangeBasis3D
USE MOD_Interpolation            ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars       ,ONLY: NodeTypeCL,NodeType
USE MOD_Mesh_Vars                ,ONLY: NGeo,NGeoRef
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Mesh_Vars                ,ONLY: nElems
USE MOD_Particle_Mesh_Vars       ,ONLY: nComputeNodeElems
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemsJ,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars       ,ONLY: RefMappingEps
USE MOD_Mesh_Tools               ,ONLY: GetGlobalElemID
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoRef
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: dXCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ_Shared,ElemEpsOneCell_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ_Shared_Win,ElemEpsOneCell_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,firstElem,lastElem
REAL                           :: scaleJ,maxScaleJ
#if USE_MPI
INTEGER                        :: ElemID
INTEGER                        :: i,j,k
! Vandermonde matrices
REAL                           :: Vdm_CLNGeo_NGeoRef(0:NgeoRef,0:NGeo)
REAL                           :: Vdm_NGeoRef_N(     0:PP_N   ,0:NGeoRef)
! Jacobian on CL N and NGeoRef
REAL                           :: detJac_Ref(1  ,0:NGeoRef,0:NGeoRef,0:NGeoRef)
REAL                           :: DetJac_N(  1  ,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on CL N
REAL                           :: dX_NGeoRef(3,3,0:NGeoRef,0:NGeoRef,0:NGeoRef)

INTEGER                        :: ElemLocID
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Building EpsOneCell for all elements ...'

! build sJ for all elements not on local proc
#if USE_MPI
MPISharedSize = INT(((PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/(PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems/),ElemsJ_Shared_Win,ElemsJ_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemsJ_Shared_Win,IERROR)
ElemsJ(0:PP_N,0:PP_N,0:PP_N,1:nComputeNodeTotalElems) => ElemsJ_Shared

IF (myComputeNodeRank.EQ.0) THEN
  ElemsJ_Shared = 0.
END IF

CALL MPI_WIN_SYNC(ElemsJ_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

#if USE_MPI
! Calculate sJ for elements not inside current proc, otherwise copy local values
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNGeo_NGeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NGeoRef_N     , modal=.TRUE.)

DO iElem = firstElem,lastElem
  ElemID    = GetGlobalElemID(iElem)
  ElemLocID = ElemID-offsetElem
  ! element on local proc, sJ already calculated in metrics.f90
  IF ((ElemLocID.GT.0) .AND. (ElemLocID.LE.nElems)) THEN
    ElemsJ(:,:,:,iElem) = sJ(:,:,:,ElemLocID)

  ! element not on local proc, calculate sJ frm dXCL_NGeo_Shared
  ELSE
    detJac_Ref = 0.
    ! Compute Jacobian on NGeo and then interpolate:
    ! required to guarantee conservativity when restarting with N<NGeo
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,1,:,:,:,iElem),dX_NGeoRef(:,1,:,:,:))
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,2,:,:,:,iElem),dX_NGeoRef(:,2,:,:,:))
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,3,:,:,:,iElem),dX_NGeoRef(:,3,:,:,:))
    DO k=0,NGeoRef; DO j=0,NGeoRef; DO i=0,NGeoRef
      detJac_Ref(1,i,j,k)=detJac_Ref(1,i,j,k) &
        + dX_NGeoRef(1,1,i,j,k)*(dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(3,3,i,j,k) - dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(2,3,i,j,k))  &
        + dX_NGeoRef(2,1,i,j,k)*(dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(1,3,i,j,k) - dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(3,3,i,j,k))  &
        + dX_NGeoRef(3,1,i,j,k)*(dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(2,3,i,j,k) - dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(1,3,i,j,k))
    END DO; END DO; END DO !i,j,k=0,NgeoRef

    ! interpolate detJac_ref to the solution points
    CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:),DetJac_N)

    ! assign to global Variable sJ
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ElemsJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
    END DO; END DO; END DO !i,j,k=0,PP_N
  END IF
END DO

CALL MPI_WIN_SYNC(ElemsJ_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ElemsJ => sJ
#endif /* USE_MPI*/
IF (TrackingMethod.EQ.TRIATRACKING) RETURN
! allocate epsOneCell
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemEpsOneCell_Shared_Win,ElemEpsOneCell_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemEpsOneCell_Shared_Win,IERROR)
ElemEpsOneCell => ElemEpsOneCell_Shared
#else
ALLOCATE(ElemEpsOneCell(1:nComputeNodeElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemEpsOneCell = -1.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemEpsOneCell_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

maxScaleJ = 0.
DO iElem = firstElem,lastElem
  scaleJ = MAXVAL(ElemsJ(:,:,:,iElem))/MINVAL(ElemsJ(:,:,:,iElem))
  ElemepsOneCell(iElem) = 1.0 + SQRT(3.0*scaleJ*RefMappingEps)
  maxScaleJ  =MAX(scaleJ,maxScaleJ)
END DO ! iElem = firstElem,lastElem

!IF(CalcMeshInfo)THEN
!  CALL AddToElemData(ElementOut,'epsOneCell',RealArray=epsOneCell(1:nElems))
!END IF

END SUBROUTINE BuildEpsOneCell


SUBROUTINE GetLinearSideBaseVectors()
!===================================================================================================================================
! computes the face base vector for linear (planar or bilinear) face intersection calculation
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,BaseVectorsScale
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectorsScale_Shared,BaseVectorsScale_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectors0_Shared,BaseVectors1_Shared,BaseVectors2_Shared,BaseVectors3_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectors0_Shared_Win,BaseVectors1_Shared_Win,BaseVectors2_Shared_Win,BaseVectors3_Shared_Win
#endif /* USE_MPI */
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,firstSide,lastSide
REAL                           :: crossVec(3)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' GET LINEAR SIDE BASEVECTORS...'
#if USE_MPI
MPISharedSize = INT((3*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors0_Shared_Win,BaseVectors0_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors0_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors1_Shared_Win,BaseVectors1_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors1_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors2_Shared_Win,BaseVectors2_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors2_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors3_Shared_Win,BaseVectors3_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors3_Shared_Win,IERROR)
MPISharedSize = INT((nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),BaseVectorsScale_Shared_Win,BaseVectorsScale_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectorsScale_Shared_Win,IERROR)
BaseVectors0 => BaseVectors0_Shared
BaseVectors1 => BaseVectors1_Shared
BaseVectors2 => BaseVectors2_Shared
BaseVectors3 => BaseVectors3_Shared
BaseVectorsScale => BaseVectorsScale_Shared

firstSide = INT(REAL (myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE( BaseVectors0(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors1(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors2(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors3(1:3,1:nNonUniqueGlobalSides),&
          BaseVectorsScale(1:nNonUniqueGlobalSides))

firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,lastSide
    BaseVectors0(:,iSide) = (+BezierControlPoints3DElevated(:,0,0           ,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             +BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3DElevated(:,0,0           ,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             -BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3DElevated(:,0,0           ,iSide)-BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             +BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3DElevated(:,0,0           ,iSide)-BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             -BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
ELSE
  DO iSide=firstSide,lastSide
    BaseVectors0(:,iSide) = (+BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
END IF

#if USE_MPI
CALL MPI_WIN_SYNC(BaseVectors0_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors1_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors2_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors3_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectorsScale_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */

SWRITE(UNIT_stdOut,'(A)')' GET LINEAR SIDE BASEVECTORS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE GetLinearSideBaseVectors


SUBROUTINE GetMeshMinMax()
!===================================================================================================================================
! computes the minimum and maximum value of the mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Mesh_Vars               ,ONLY: offsetElem,nElems
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeSide,nComputeNodeSides
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeNode,nComputeNodeNodes
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: offsetLocalSide,nLocalSides
INTEGER                        :: offsetLocalNode,nLocalNodes
!===================================================================================================================================

SELECT CASE(TrackingMethod)
  ! Build mesh min/max on BezierControlPoints for possibly curved elements
  CASE(REFMAPPING,TRACING)
    ! calculate all offsets
    offsetLocalSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1)
    nLocalSides     = ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1)

    ! proc local
    GEO%xmin     = MINVAL(BezierControlPoints3D(1,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%xmax     = MAXVAL(BezierControlPoints3D(1,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%ymin     = MINVAL(BezierControlPoints3D(2,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%ymax     = MAXVAL(BezierControlPoints3D(2,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%zmin     = MINVAL(BezierControlPoints3D(3,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%zmax     = MAXVAL(BezierControlPoints3D(3,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))

#if USE_MPI
    ! compute-node local
    GEO%CNxmin   = MINVAL(BezierControlPoints3D(1,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNxmax   = MAXVAL(BezierControlPoints3D(1,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNymin   = MINVAL(BezierControlPoints3D(2,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNymax   = MAXVAL(BezierControlPoints3D(2,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNzmin   = MINVAL(BezierControlPoints3D(3,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNzmax   = MAXVAL(BezierControlPoints3D(3,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))

    ! global
    GEO%xminglob = MINVAL(BezierControlPoints3D(1,:,:,:))
    GEO%xmaxglob = MAXVAL(BezierControlPoints3D(1,:,:,:))
    GEO%yminglob = MINVAL(BezierControlPoints3D(2,:,:,:))
    GEO%ymaxglob = MAXVAL(BezierControlPoints3D(2,:,:,:))
    GEO%zminglob = MINVAL(BezierControlPoints3D(3,:,:,:))
    GEO%zmaxglob = MAXVAL(BezierControlPoints3D(3,:,:,:))
#endif /*USE_MPI*/
  ! TriaTracking does not have curved elements, nodeCoords are sufficient
  CASE(TRIATRACKING)
    ! calculate all offsets
    offsetLocalNode = ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)
    nLocalNodes     = ElemInfo_Shared(ELEM_LASTNODEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)

    ! proc local
    GEO%xmin     = MINVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%xmax     = MAXVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%ymin     = MINVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%ymax     = MAXVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%zmin     = MINVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%zmax     = MAXVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))

#if USE_MPI
    ! compute-node local
    GEO%CNxmin   = MINVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNxmax   = MAXVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNymin   = MINVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNymax   = MAXVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNzmin   = MINVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNzmax   = MAXVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))

    ! global
    GEO%xminglob = MINVAL(NodeCoords_Shared(1,:))
    GEO%xmaxglob = MAXVAL(NodeCoords_Shared(1,:))
    GEO%yminglob = MINVAL(NodeCoords_Shared(2,:))
    GEO%ymaxglob = MAXVAL(NodeCoords_Shared(2,:))
    GEO%zminglob = MINVAL(NodeCoords_Shared(3,:))
    GEO%zmaxglob = MAXVAL(NodeCoords_Shared(3,:))
#endif /*USE_MPI*/
END SELECT

#if !USE_MPI
! compute-node local (dummy)
GEO%CNxmin   = GEO%xmin
GEO%CNxmax   = GEO%xmax
GEO%CNymin   = GEO%ymin
GEO%CNymax   = GEO%ymax
GEO%CNzmin   = GEO%zmin
GEO%CNzmax   = GEO%zmax

! global (dummy)
GEO%xminglob = GEO%xmin
GEO%xmaxglob = GEO%xmax
GEO%yminglob = GEO%ymin
GEO%ymaxglob = GEO%ymax
GEO%zminglob = GEO%zmin
GEO%zmaxglob = GEO%zmax
#endif /*USE_MPI*/

END SUBROUTINE GetMeshMinMax


SUBROUTINE BuildNodeNeighbourhood()
!----------------------------------------------------------------------------------------------------------------------------------!
! Build shared arrays for mapping
! 1. UniqueNodeID -> all connected CN element IDs:    NodeToElemMapping(1,:) = OffsetNodeToElemMapping(:)
!                                                     NodeToElemMapping(2,:) = NbrOfElemsOnUniqueNode(:)
!
!                                                     NodeToElemInfo(nNodeToElemMapping)
!
! 2. CN elment ID -> all connected CN element IDs:    ElemToElemMapping(1,:) = Offset
!                                                     ElemToElemMapping(2,:) = Number of elements
!
!                                                     ElemToElemInfo(nElemToElemMapping)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals            ,ONLY: abort!,myRank
USE MOD_Particle_Mesh_Vars ,ONLY: nUniqueGlobalNodes
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping,NodeToElemInfo,ElemToElemMapping,ElemToElemInfo
#if USE_MPI
USE MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping,NodeToElemInfo,ElemToElemMapping,ElemToElemInfo
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping_Shared,NodeToElemInfo_Shared,ElemToElemMapping_Shared,ElemToElemInfo_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: NodeToElemMapping_Shared_Win,NodeToElemInfo_Shared_Win
USE MOD_Particle_Mesh_Vars ,ONLY: ElemToElemMapping_Shared_Win,ElemToElemInfo_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE            :: NbrOfElemsOnUniqueNode(:),OffsetNodeToElemMapping(:)
INTEGER                        :: UniqueNodeID,NonUniqueNodeID,iElem,iNode,TestElemID,jElem,kElem
INTEGER                        :: OffsetCounter,OffsetElemToElemMapping,OffsetElemToElemCounter
INTEGER                        :: nNodeToElemMapping,iUniqueNode,firstElem,lastElem,iError,nElemToElemMapping,CountElems
INTEGER,ALLOCATABLE            :: CheckedElemIDs(:)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

! 1.1 Get number of CN elements attached to each UNIQUE node and store in NbrOfElemsOnUniqueNode(UniqueNodeID)
! 1.2 Store the total number of counted elements in nNodeToElemMapping = SUM(NbrOfElemsOnUniqueNode)
! 1.3 Store the element offset for each UNIQUE node in OffsetNodeToElemMapping(iUniqueNode) by summing up the number of elements
!     from the first the previous (iUniqueNode-1) node

! Only the node leader fills the arrays
#if USE_MPI
IF(myComputeNodeRank.EQ.0)THEN
#endif /*USE_MPI*/
  ALLOCATE(NbrOfElemsOnUniqueNode(nUniqueGlobalNodes))
  NbrOfElemsOnUniqueNode=0
  ! Loop all CN elements (iElem is CNElemID)
#if USE_MPI
  DO iElem = 1,nComputeNodeTotalElems
#else
  DO iElem = 1,nElems
#endif /*USE_MPI*/
    ! Loop all local nodes
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
      UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
      NbrOfElemsOnUniqueNode(UniqueNodeID) = NbrOfElemsOnUniqueNode(UniqueNodeID) + 1
    END DO ! iNode = 1, 8
  END DO ! iElem = 1, nComputeNodeTotalElems
  nNodeToElemMapping = SUM(NbrOfElemsOnUniqueNode)

  ALLOCATE(OffsetNodeToElemMapping(nUniqueGlobalNodes))
  OffsetNodeToElemMapping = 0
  DO iUniqueNode = 2, nUniqueGlobalNodes
    OffsetNodeToElemMapping(iUniqueNode) = SUM(NbrOfElemsOnUniqueNode(1:iUniqueNode-1))
  END DO ! iUniqueNode = 1, nUniqueGlobalNodes
#if USE_MPI
END IF ! myComputeNodeRank.EQ.0
CALL MPI_BCAST(nNodeToElemMapping,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/


! 2. Allocate shared arrays for mapping
!    UniqueNodeID -> all CN element IDs to which it is connected      : NodeToElemMapping = [OffsetNodeToElemMapping(:),
!                                                                                            NbrOfElemsOnUniqueNode(:)]
!    NodeToElemMapping (offset and number of elements) -> CN element IDs : NodeToElemInfo = [CN elem IDs]
#if USE_MPI
! NodeToElemMapping
MPISharedSize = INT((2*nUniqueGlobalNodes),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,nUniqueGlobalNodes/),NodeToElemMapping_Shared_Win,NodeToElemMapping_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeToElemMapping_Shared_Win,IERROR)
NodeToElemMapping => NodeToElemMapping_Shared
! NodeToElemInfo
MPISharedSize = INT((nNodeToElemMapping),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nNodeToElemMapping/),NodeToElemInfo_Shared_Win,NodeToElemInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeToElemInfo_Shared_Win,IERROR)
NodeToElemInfo => NodeToElemInfo_Shared
#else
ALLOCATE(NodeToElemMapping(2,nUniqueGlobalNodes))
ALLOCATE(NodeToElemInfo(nNodeToElemMapping))
#endif /*USE_MPI*/


! 3.1 Fill NodeToElemMapping = [OffsetNodeToElemMapping(:), NbrOfElemsOnUniqueNode(:)]
! 3.2 Store the CN element IDs in NodeToElemInfo(NodeToElemMapping(1,UniqueNodeID)+1 :
!                                                   NodeToElemMapping(1,UniqueNodeID)+NodeToElemMapping(2,UniqueNodeID))
! Now all CN elements attached to a UniqueNodeID can be accessed

! Only the node leader fills the arrays
#if USE_MPI
IF(myComputeNodeRank.EQ.0)THEN
#endif /*USE_MPI*/
  NodeToElemMapping = 0
  NodeToElemInfo = 0

  NodeToElemMapping(1,:) = OffsetNodeToElemMapping(:)
  DEALLOCATE(OffsetNodeToElemMapping)
  NodeToElemMapping(2,:) = NbrOfElemsOnUniqueNode(:)

  NbrOfElemsOnUniqueNode=0
  ! Loop all CN elements (iElem is CNElemID)
#if USE_MPI
  DO iElem = 1, nComputeNodeTotalElems
#else
  DO iElem = 1, nElems
#endif
    ! Loop all local nodes
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
      UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
      NbrOfElemsOnUniqueNode(UniqueNodeID) = NbrOfElemsOnUniqueNode(UniqueNodeID) + 1
      NodeToElemInfo(NodeToElemMapping(1,UniqueNodeID)+NbrOfElemsOnUniqueNode(UniqueNodeID)) = iElem
    END DO ! iNode = 1, 8
  END DO ! iElem = 1, nComputeNodeTotalElems
  DEALLOCATE(NbrOfElemsOnUniqueNode)
#if USE_MPI
END IF ! myComputeNodeRank.EQ.0

CALL MPI_WIN_SYNC(NodeToElemInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(NodeToElemMapping_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! 4. Allocate shared array for mapping
!    CN element ID -> all CN element IDs to which it is connected : ElemToElemMapping = [offset, Nbr of CN elements]
#if USE_MPI
! ElemToElemMapping
MPISharedSize = INT((2*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,nComputeNodeTotalElems/),ElemToElemMapping_Shared_Win,ElemToElemMapping_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToElemMapping_Shared_Win,IERROR)
ElemToElemMapping => ElemToElemMapping_Shared
#else
ALLOCATE(ElemToElemMapping(2,nElems))
#endif /*USE_MPI*/

! 5. Fill ElemToElemMapping = [offset, Nbr of CN elements]
!    Note that the number of elements stored in ElemToElemMapping(2,iElem) must be shifted after communication with other procs
#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

OffsetCounter = 0
ALLOCATE(CheckedElemIDs(500))

! Loop all CN elements (iElem is CNElemID)
DO iElem = firstElem,lastElem
  CountElems = 0
  CheckedElemIDs = 0
  ! Loop all local nodes
  DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)

    ! Loop 1D array [offset + 1 : offset + NbrOfElems]
    ! (all CN elements that are connected to the local nodes)
    Elemloop: DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
      TestElemID = NodeToElemInfo(jElem)

      IF(iElem.EQ.TestElemID) CYCLE Elemloop

      ! Check previously stored element IDs and cycle if an already stored ID is encountered
      DO kElem = 1, CountElems
        IF(CheckedElemIDs(kElem).EQ.TestElemID) CYCLE Elemloop
      END DO ! kElem = 1, CountElems

      CountElems = CountElems + 1

      IF(CountElems.GT.500) CALL abort(&
      __STAMP__&
      ,'CountElems > 500. Inrease the number and try again!')

      CheckedElemIDs(CountElems) = TestElemID
      ! Note that the number of elements stored in ElemToElemMapping(2,iElem) must be shifted after communication with other procs
      ElemToElemMapping(2,iElem) = CountElems
    END DO Elemloop
  END DO ! iNode = 1, 8
  ElemToElemMapping(1,iElem) = OffsetCounter
  OffsetCounter = OffsetCounter + CountElems
END DO ! iElem = firstElem, lastElem

! 6. Find CN global number of connected CN elements (=nElemToElemMapping) and write into shared array
#if USE_MPI
sendbuf = OffsetCounter
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
OffsetElemToElemMapping   = recvbuf
! last proc knows CN total number of connected CN elements
sendbuf = OffsetElemToElemMapping + OffsetCounter
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nElemToElemMapping = sendbuf

!ElemToElemMapping(1,firstElem:lastElem) = ElemToElemMapping(1,firstElem:lastElem)
ElemToElemMapping(1,firstElem:lastElem) = ElemToElemMapping(1,firstElem:lastElem) + OffsetElemToElemMapping
#else
OffsetElemToElemMapping = 0
nElemToElemMapping = OffsetCounter
!ElemToElemMapping(:,firstElem:lastElem) = ElemToElemMapping(:,firstElem:lastElem)
#endif /*USE_MPI*/


! 7. Allocate shared array for mapping
!    CN element ID -> all CN element IDs to which it is connected : ElemToElemInfo = [CN elem IDs]
#if USE_MPI
! ElemToElemInfo
MPISharedSize = INT((nElemToElemMapping),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nElemToElemMapping/),ElemToElemInfo_Shared_Win,ElemToElemInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToElemInfo_Shared_Win,IERROR)
ElemToElemInfo => ElemToElemInfo_Shared
#else
ALLOCATE(ElemToElemInfo(nElemToElemMapping))
#endif /*USE_MPI*/


! 8. Fill ElemToElemInfo = [CN elem IDs] by finding all nodes connected to a CN element
!    (and all subsequent CN elements that are connected to those nodes)
!    Store the CN element IDs in ElemToElemInfo(ElemToElemMapping(1,iElem)+1 :
!                                               ElemToElemMapping(1,iElem)+ElemToElemMapping(2,iElem))
OffsetElemToElemCounter = OffsetElemToElemMapping
! Loop all CN elements (iElem is CNElemID)
DO iElem = firstElem, lastElem
  CountElems = 0
  CheckedElemIDs = 0
  ! Loop all local nodes
  DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)

    ! Loop all CN elements that are connected to the local nodes
    Elemloop2: DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
      TestElemID = NodeToElemInfo(jElem)

      IF(iElem.EQ.TestElemID) CYCLE Elemloop2

      DO kElem = 1, CountElems
        IF(CheckedElemIDs(kElem).EQ.TestElemID) CYCLE Elemloop2
      END DO ! kElem = 1, CountElems

      CountElems = CountElems + 1
      OffsetElemToElemCounter = OffsetElemToElemCounter + 1

      IF(CountElems.GT.500) CALL abort(&
      __STAMP__&
      ,'CountElems > 500. Inrease the number and try again!')

      CheckedElemIDs(CountElems) = TestElemID
      ElemToElemInfo(OffsetElemToElemCounter) = TestElemID

    END DO ElemLoop2
  END DO ! iNode = 1, 8

END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL MPI_WIN_SYNC(ElemToElemInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemToElemMapping_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildNodeNeighbourhood


SUBROUTINE MarkAuxBCElems()
!===================================================================================================================================
! check if auxBCs are inside BoundingBox of Elems
! -- plane: use plane equation f=a1*x+a2*y+a3*z+a4=0 and insert corresponding intervals of box -> fmin and fmax
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemHasAuxBCs
USE MOD_Particle_Mesh_Vars     ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Boundary_Vars ,ONLY: nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone
#if USE_MPI
USE MOD_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,iAuxBC,icoord,dir(3),positiontype,positiontype_tmp
REAL                     :: r_vec(3),n_vec(3),fmin,fmax,radius,BoundsBC(1:2,1:3)
REAL                     :: lmin,lmax,deltamin,deltamax,origin(2),halfangle
LOGICAL                  :: cartesian, backwards
!===================================================================================================================================

ALLOCATE(ElemHasAuxBCs(1:PP_nElems , 1:nAuxBCs))
ElemHasAuxBCs=.FALSE.

DO iAuxBC=1,nAuxBCs
  SELECT CASE (TRIM(AuxBCType(iAuxBC)))
  CASE ('plane')
    r_vec=AuxBC_plane(AuxBCMap(iAuxBC))%r_vec
    n_vec=AuxBC_plane(AuxBCMap(iAuxBC))%n_vec
    radius=AuxBC_plane(AuxBCMap(iAuxBC))%radius
    ! loop over all  elements
    DO iElem=1,PP_nElems
      ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
        fmin=-DOT_PRODUCT(r_vec,n_vec)
        fmax=fmin
        DO icoord=1,3
          IF (n_vec(icoord).GE.0) THEN
            fmin = fmin + n_vec(icoord)*Bounds(1,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(2,icoord)
          ELSE
            fmin = fmin + n_vec(icoord)*Bounds(2,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(1,icoord)
          END IF
        END DO
        IF ((fmin.LE.0 .AND. fmax.GT.0).OR.(fmin.LT.0 .AND. fmax.GE.0)) THEN !plane intersects the box!
          !radius check needs to be implemented (compute intersection polygon and minimum radii): would sort out further elements!!!
          !quick, conservative solution: calculate bounding box of disc in space and compare with bb of element
          ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
          IF (radius .LT. 0.5*HUGE(radius)) THEN !huge was default
            BoundsBC(1,1:3) = r_vec - radius * SQRT(1.-(n_vec*n_vec))
            BoundsBC(2,1:3) = r_vec + radius * SQRT(1.-(n_vec*n_vec))
            DO icoord=1,3
              IF ( BoundsBC(2,icoord).LT.Bounds(1,icoord) .OR. BoundsBC(1,icoord).GT.Bounds(2,icoord) ) THEN
                ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
                EXIT
              END IF
            END DO
          END IF
        ELSE IF ((fmin.LT.0 .AND. fmax.LT.0).OR.(fmin.GT.0 .AND. fmax.GT.0)) THEN !plane does not intersect the box!
          ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
        ELSE !e.g. if elem has zero volume...
          CALL abort(&
            __STAMP__&
            ,'Error in MarkAuxBCElems for AuxBC:',iAuxBC)
        END IF
      END ASSOCIATE
    END DO
  CASE ('cylinder','cone')
    IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
      r_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec
      n_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%axis
      radius=AuxBC_cylinder(AuxBCMap(iAuxBC))%radius
      lmin=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin
      lmax=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax
    ELSE !cone
      r_vec=AuxBC_cone(AuxBCMap(iAuxBC))%r_vec
      n_vec=AuxBC_cone(AuxBCMap(iAuxBC))%axis
      halfangle=AuxBC_cone(AuxBCMap(iAuxBC))%halfangle
      lmin=AuxBC_cone(AuxBCMap(iAuxBC))%lmin
      lmax=AuxBC_cone(AuxBCMap(iAuxBC))%lmax
    END IF
    cartesian=.TRUE.
    backwards=.FALSE.
    IF (ABS(n_vec(1)).EQ.1.) THEN
      dir(1)=1
      dir(2)=2
      dir(3)=3
      IF (n_vec(1).LT.0.) backwards=.TRUE.
    ELSE IF (ABS(n_vec(2)).EQ.1.) THEN
      dir(1)=2
      dir(2)=3
      dir(3)=1
      IF (n_vec(2).LT.0.) backwards=.TRUE.
    ELSE IF (ABS(n_vec(3)).EQ.1.) THEN
      dir(1)=3
      dir(2)=1
      dir(3)=2
      IF (n_vec(3).LT.0.) backwards=.TRUE.
    ELSE
      cartesian=.FALSE.
      SWRITE(*,*) 'WARNING in MarkAuxBCElems: all Elems are set to ElemHasAuxBCs=.TRUE. for AuxBC:',iAuxBC
      ElemHasAuxBCs(:,iAuxBC)=.TRUE. !actual intersection with box check to-be implemented!!!
    END IF
    IF (cartesian) THEN
      IF (backwards) THEN
        deltamin = -lmax
        deltamax = -lmin
      ELSE
        deltamin = lmin
        deltamax = lmax
      END IF
      origin(1) = r_vec(dir(2))
      origin(2) = r_vec(dir(3))
      ! loop over all  elements
      DO iElem=1,PP_nElems
        ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
          ! check for lmin and lmax
          IF ( r_vec(dir(1))+deltamax.LT.Bounds(1,dir(1)) .OR. r_vec(dir(1))+deltamin.GT.Bounds(2,dir(1)) ) THEN
            ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
          ELSE !between lmin and lmax
            IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
            ELSE !cone
              !local minimum radius
              IF (backwards) THEN
                radius = MAX(-Bounds(2,dir(1))+r_vec(dir(1)),lmin)*TAN(halfangle)
              ELSE
                radius = MAX(Bounds(1,dir(1))-r_vec(dir(1)),lmin)*TAN(halfangle)
              END IF
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype_tmp)
              !local maximum radius
              IF (backwards) THEN
                radius = MIN(-Bounds(1,dir(1))+r_vec(dir(1)),lmax)*TAN(halfangle)
              ELSE
                radius = MIN(Bounds(2,dir(1))-r_vec(dir(1)),lmax)*TAN(halfangle)
              END IF
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
              !if both are type 0 or both are type 1 than the "total" type is not 2:
              IF ( .NOT.(positiontype_tmp.EQ.0 .AND. positiontype.EQ.0) &
                .AND. .NOT.(positiontype_tmp.EQ.1 .AND. positiontype.EQ.1) ) THEN
                positiontype=2
              END IF
            END IF
            IF (positiontype.EQ.2) THEN
              ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
            ELSE
              ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
            END IF
          END IF !check for lmin and lmax
        END ASSOCIATE
      END DO !iElem
    END IF !cartesian
  CASE('parabol')
    ElemHasAuxBCs(:,iAuxBC)=.TRUE. ! to be implemented!!!
  CASE DEFAULT
    SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
    CALL abort(&
      __STAMP__&
      ,'AuxBC does not exist')
  END SELECT
END DO

END SUBROUTINE MarkAuxBCElems


SUBROUTINE CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
!===================================================================================================================================
! checks how a cartesian bb is located with regard to a radius with cartesian axis (dir is cartesian axis and origin in orth. dirs)
!- positiontype=0 : complete bb is inside of radius
!- positiontype=1 : complete bb is outside of radius
!- positiontype=2 : bb is partly inside of radius
! (based on "check where the sides are located relative to rmax" in particle_emission for SimpleRadialVeloFit)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)           :: Bounds(1:2,1:3), origin(2), radius
INTEGER,INTENT(IN)        :: dir(3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: positiontype
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iDir1, iDir2, iDir3, iPoint
REAL                      :: BoundingBox(1:3,1:8), point(2), pointRadius
LOGICAL                   :: done, insideBound
!===================================================================================================================================
!-- convert minmax-values to bb-points
DO iDir1=0,1
  DO iDir2=0,1
      DO iDir3=0,1
        BoundingBox(1,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir1+1,1)
        BoundingBox(2,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir2+1,2)
        BoundingBox(3,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir3+1,3)
      END DO
  END DO
END DO

!-- check where the points are located relative to radius
done=.FALSE.
DO iDir1=0,1
  IF(done) EXIT
  DO iDir2=0,1
    IF(done) EXIT
    DO iDir3=0,1
      !-- coords orth. to axis of point:
      iPoint=iDir1*4 + iDir2*2 + iDir3+1
      point(1) = BoundingBox(dir(2),iPoint)-origin(1)
      point(2) = BoundingBox(dir(3),iPoint)-origin(2)
      pointRadius = SQRT( (point(1))**2+(point(2))**2 )
      IF (iPoint.EQ.1) THEN
        IF (pointRadius.LE.radius) THEN
          insideBound=.TRUE.
        ELSE !outside
          insideBound=.FALSE.
        END IF !in-/outside?
      ELSE !iPoint.GT.1: type must be 2 if state of point if different from last point
        IF (pointRadius.LE.radius) THEN
          IF (.NOT.insideBound) THEN !different from last point
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        ELSE !outside
          IF (insideBound) THEN !different from last point
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        END IF !in-/outside?
      END IF !iPoint.EQ.1
    END DO !iDir3
  END DO !iDir2
END DO !iDir1
IF (.NOT.done) THEN
  IF (insideBound) THEN
    positiontype=0
  ELSE
    ! all points are outside of radius, but when radius is smaller than box, it can intersect it:
    IF ( origin(1) + radius .GE. Bounds(1,dir(2)) .AND. &
         origin(1) - radius .LE. Bounds(2,dir(2)) .AND. &
         origin(2) + radius .GE. Bounds(1,dir(3)) .AND. &
         origin(2) - radius .LE. Bounds(2,dir(3)) ) THEN !circle completely or partly inside box
      positiontype=2
    ELSE !points are really outside
      positiontype=1
    END IF
  END IF
END IF

END SUBROUTINE CheckBoundsWithCartRadius


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
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#if USE_MPI
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
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

    ! GetBCSidesAndOrgin
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL UNLOCK_AND_FREE(BCSide2SideID_Shared_Win)
      CALL UNLOCK_AND_FREE(SideID2BCSide_Shared_Win)
      CALL UNLOCK_AND_FREE(BCSideMetrics_Shared_Win)
    END IF

#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! CalcParticleMeshMetrics
      CALL UNLOCK_AND_FREE(XCL_NGeo_Shared_Win)
      CALL UNLOCK_AND_FREE(Elem_xGP_Shared_Win)
      CALL UNLOCK_AND_FREE(dXCL_NGeo_Shared_Win)

      ! CalcBezierControlPoints
      CALL UNLOCK_AND_FREE(BezierControlPoints3D_Shared_Win)
      IF (BezierElevation.GT.0) CALL UNLOCK_AND_FREE(BezierControlPoints3DElevated_Shared_Win)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! GetSideSlabNormalsAndIntervals (allocated in particle_mesh.f90)
    CALL UNLOCK_AND_FREE(SideSlabNormals_Shared_Win)
    CALL UNLOCK_AND_FREE(SideSlabIntervals_Shared_Win)
    CALL UNLOCK_AND_FREE(BoundingBoxIsEmpty_Shared_Win)

    ! BuildElementOriginShared
    CALL UNLOCK_AND_FREE(ElemBaryNGeo_Shared_Win)

    ! IdentifyElemAndSideType
    CALL UNLOCK_AND_FREE(ElemCurved_Shared_Win)
    CALL UNLOCK_AND_FREE(SideType_Shared_Win)
    CALL UNLOCK_AND_FREE(SideDistance_Shared_Win)
    CALL UNLOCK_AND_FREE(SideNormVec_Shared_Win)

    ! BuildElementBasisAndRadius
    CALL UNLOCK_AND_FREE(ElemRadiusNGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemRadius2NGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(XiEtaZetaBasis_Shared_Win)
    CALL UNLOCK_AND_FREE(slenXiEtaZetaBasis_Shared_Win)

    ! GetLinearSideBaseVectors
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

    ! BuildEpsOneCell
    CALL UNLOCK_AND_FREE(ElemEpsOneCell_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemsJ_Shared_Win)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
    ! GetBCSidesAndOrgin
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
      ! CalcParticleMeshMetrics
      ADEALLOCATE(XCL_NGeo_Array)
      ADEALLOCATE(Elem_xGP_Array)
      ADEALLOCATE(dXCL_NGeo_Array)
      IF(ASSOCIATED(XCL_NGeo_Shared))  NULLIFY(XCL_NGeo_Shared)
      IF(ASSOCIATED(Elem_xGP_Shared))  NULLIFY(Elem_xGP_Shared)
      IF(ASSOCIATED(dXCL_NGeo_Shared)) NULLIFY(dXCL_NGeo_Shared)

      ! CalcBezierControlPoints
      ADEALLOCATE(BezierControlPoints3D_Shared)
      ADEALLOCATE(BezierControlPoints3DElevated_Shared)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! GetSideSlabNormalsAndIntervals (allocated in particle_mesh.f90)
    ADEALLOCATE(SideSlabNormals_Shared)
    ADEALLOCATE(SideSlabIntervals_Shared)
    ADEALLOCATE(BoundingBoxIsEmpty_Shared)

    ! BuildElementOriginShared
    ADEALLOCATE(ElemBaryNGeo_Shared)

    ! IdentifyElemAndSideType
    ADEALLOCATE(ElemCurved)
    ADEALLOCATE(ElemCurved_Shared)
    ADEALLOCATE(SideType_Shared)
    ADEALLOCATE(SideDistance_Shared)
    ADEALLOCATE(SideNormVec_Shared)

    ! BuildElementBasisAndRadius
    ADEALLOCATE(ElemRadiusNGeo_Shared)
    ADEALLOCATE(ElemRadius2NGeo_Shared)
    ADEALLOCATE(XiEtaZetaBasis_Shared)
    ADEALLOCATE(slenXiEtaZetaBasis_Shared)

    ! GetLinearSideBaseVectors
    ADEALLOCATE(BaseVectors0_Shared)
    ADEALLOCATE(BaseVectors1_Shared)
    ADEALLOCATE(BaseVectors2_Shared)
    ADEALLOCATE(BaseVectors3_Shared)
    ADEALLOCATE(BaseVectorsScale_Shared)

    ! BuildBCElemDistance
    IF (TrackingMethod.EQ.1) THEN
      ADEALLOCATE(ElemToBCSides_Shared)
      ADEALLOCATE(SideBCMetrics_Shared)
    END IF

    ! BuildEpsOneCell
    ADEALLOCATE(ElemsJ_Shared)
    ADEALLOCATE(ElemEpsOneCell_Shared)

!  ! Tracing
!  CASE(TRACING)

  ! TriaTracking
  CASE(TRIATRACKING)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
    !IF (DoInterpolation.OR.DSMC%UseOctree) THEN ! use this in future if possible
    IF (DoInterpolation) THEN
#if USE_LOADBALANCE
      IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
        ! CalcParticleMeshMetrics
        CALL UNLOCK_AND_FREE(XCL_NGeo_Shared_Win)
        CALL UNLOCK_AND_FREE(Elem_xGP_Shared_Win)
        CALL UNLOCK_AND_FREE(dXCL_NGeo_Shared_Win)
#if USE_LOADBALANCE
      END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

      ! BuildElemTypeAndBasisTria
      CALL UNLOCK_AND_FREE(ElemCurved_Shared_Win)
      CALL UNLOCK_AND_FREE(XiEtaZetaBasis_Shared_Win)
      CALL UNLOCK_AND_FREE(slenXiEtaZetaBasis_Shared_Win)
    END IF ! DoInterpolation

    ! InitParticleGeometry
    CALL UNLOCK_AND_FREE(ConcaveElemSide_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemSideNodeID_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemMidPoint_Shared_Win)

    ! BuildElementRadiusTria
    CALL UNLOCK_AND_FREE(ElemBaryNGeo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemRadius2NGeo_Shared_Win)

    ! BuildEpsOneCell
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

    ! BuildElemTypeAndBasisTria()
    ADEALLOCATE(XiEtaZetaBasis_Shared)
    ADEALLOCATE(slenXiEtaZetaBasis_Shared)
    ADEALLOCATE(ElemCurved)
    ADEALLOCATE(ElemCurved_Shared)

    ! BuildEpsOneCell
    SNULLIFY(ElemsJ)
    ADEALLOCATE(ElemsJ_Shared)
END SELECT

IF(TRIM(DepositionType).EQ.'shape_function_adaptive'.OR.TrackingMethod.EQ.TRIATRACKING)THEN
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  ! From InitElemNodeIDs
  CALL UNLOCK_AND_FREE(ElemNodeID_Shared_Win)

  !FindNeighbourElems = .FALSE. ! THIS IS SET TO FALSE CURRENTLY in InitParticleMesh()
  ! TODO: fix when FindNeighbourElems is not always set false
  IF(FindNeighbourElems.OR.TRIM(DepositionType).EQ.'shape_function_adaptive')THEN
    ! From BuildNodeNeighbourhood
    CALL UNLOCK_AND_FREE(NodeToElemMapping_Shared_Win)
    CALL UNLOCK_AND_FREE(NodeToElemInfo_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemToElemMapping_Shared_Win)
    CALL UNLOCK_AND_FREE(ElemToElemInfo_Shared_Win)
  END IF ! FindNeighbourElems.OR.TRIM(DepositionType).EQ.'shape_function_adaptive'

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

  ADEALLOCATE(ElemNodeID_Shared)
  IF(FindNeighbourElems.OR.TRIM(DepositionType).EQ.'shape_function_adaptive')THEN
    ADEALLOCATE(NodeToElemMapping_Shared)
    ADEALLOCATE(NodeToElemInfo_Shared)
    ADEALLOCATE(ElemToElemMapping_Shared)
    ADEALLOCATE(ElemToElemInfo_Shared)
  END IF
END IF

SDEALLOCATE(GEO%PeriodicVectors)
SDEALLOCATE(GEO%FIBGM)
SDEALLOCATE(GEO%TFIBGM)

ADEALLOCATE(XiEtaZetaBasis)
ADEALLOCATE(slenXiEtaZetaBasis)
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
