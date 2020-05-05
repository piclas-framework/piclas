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

INTERFACE InitParticleGeometry
  MODULE PROCEDURE InitParticleGeometry
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
PUBLIC::InitParticleGeometry
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

CALL prms%CreateLogicalOption( 'CountNbOfLostParts'&
  , 'Count number of lost particles during tracking that can not be found with fallbacks.','.FALSE.')
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
CALL prms%CreateLogicalOption( 'printMPINeighborWarnings'&
    ,  ' Print warning if the MPI-Halo-region between to procs are not overlapping. Only one proc find the other in halo ' &
    ,'.FALSE.')
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
USE MOD_Mesh_Vars              ,ONLY: deleteMeshPointer,NodeCoords
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Mesh_Vars              ,ONLY: useCurveds
USE MOD_Particle_BGM           ,ONLY: BuildBGMAndIdentifyHaloRegion
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitGetGlobalElemID,InitGetCNElemID,InitPEM_LocalElemID,InitPEM_CNElemID
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated,SideSlabNormals,SideSlabIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping,MeasureTrackTime,FastPeriodic,CountNbOfLostParts,nLostParts,CartesianPeriodic
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod, WriteTriaDebugMesh
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GetRealArray, GETINTFROMSTR
USE MOD_Particle_Vars          ,ONLY: Symmetry2D
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars
#endif /* USE_MPI */
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
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH ...'
IF(ParticleMeshInitIsDone) CALL abort(&
__STAMP__&
, ' Particle-Mesh is already initialized.')

! Build BGM to Element mapping and identify which of the elements, sides and nodes are in the compute-node local and halo region
CALL BuildBGMAndIdentifyHaloRegion()

! Initialize mapping function: GetGlobalElemID(1:nComputeNodeTotalElems)
CALL InitGetGlobalElemID()

! Initialize mapping function: GetCNElemID(1:GlobalElemID)
CALL InitGetCNElemID()

! Initialize mapping function: PEM%LocalElemID(1:PDM%ParticleVecLength)
CALL InitPEM_LocalElemID()

! Initialize mapping function: PEM%CNElemID(1:PDM%ParticleVecLength)
CALL InitPEM_CNElemID()

!IF ((TrackingMethod.EQ.REFMAPPING.OR.UseCurveds.OR.(NGeo.GT.1)).AND.(TrackingMethod.EQ.TRIATRACKING)) THEN
!  CALL CollectiveStop(__STAMP__, &
!         'TrackingMethod=REFMAPPING .OR. UseCurveds=T .OR. NGEO>1! Not possible with TrackingMethod=TRIATRACKING at the same time!')
!END IF

CountNbOfLostParts = GETLOGICAL('CountNbOfLostParts',".FALSE.")
nLostParts         = 0

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
  IF(PP_N.GT.NGeo)THEN ! there are more Gauss points within an element then CL-points
                       ! Gauss points sample the element finer
                       ! Note: the Gauss points does not exist for HALO elements, here, the closest CL point is used
    RefMappingGuessProposal=2
  ELSE ! more CL-points than Gauss points, hence, better sampling of the element
    RefMappingGuessProposal=3
  END IF
ELSE
  RefMappingGuessProposal=1 ! default for linear meshes. Guess is exact for cubical, non-twisted elements
END IF
WRITE(tmpStr,'(I2.2)') RefMappingGuessProposal
RefMappingGuess = GETINT('RefMappingGuess',tmpStr)
IF((RefMappingGuess.LT.1).AND.(UseCurveds)) THEN ! this might cause problems
  SWRITE(UNIT_stdOut,'(A)')' WARNING: read-in [RefMappingGuess=1] when using [UseCurveds=T] may create problems!'
END IF
RefMappingEps   = GETREAL('RefMappingEps','1e-4')

epsInCell       = SQRT(3.0*RefMappingEps)
!epsOneCell      = 1.0+epsInCell

IF((RefMappingGuess.LT.1).OR.(RefMappingGuess.GT.4))THEN
   CALL abort(&
__STAMP__ &
,'Wrong guessing method for mapping from physical space in reference space.',RefMappingGuess,999.)
END IF
!IF(DoRefMapping .AND. RefMappingGuess.EQ.2) THEN
!   CALL abort(&
!__STAMP__ &
!,' No-Elem_xGP allocated for Halo-Cells! Select other mapping guess',RefMappingGuess)
!END IF

SELECT CASE(TrackingMethod)

  CASE(TRIATRACKING)
    CALL InitParticleGeometry()
    ! Compute convex element radius^2
    CALL BuildElementRadiusTria()

    ! Interpolation needs coordinates in reference system
    IF (DoInterpolation) THEN
      CALL CalcParticleMeshMetrics()

      CALL BuildElemTypeAndBasisTria()
    END IF

CASE(TRACING,REFMAPPING)
    CALL CalcParticleMeshMetrics()

    BezierElevation = GETINT('BezierElevation')
    NGeoElevated    = NGeo + BezierElevation
    CALL CalcBezierControlPoints()

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
        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE
        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide)     &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                          &
                                           ,SideSlabInterVals( 1:6    ,iSide)                          &
                                           ,BoundingBoxIsEmpty(iSide))
      END DO
    ELSE
      DO iSide=firstSide,LastSide
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
    offsetSideID = ElemInfo_Shared(SideIf
    DO iSide=offsetMPISides_YOUR,LastSide
      dx=ABS(SideSlabIntervals(2)-SideSlabIntervals(1))
      dy=ABS(SideSlabIntervals(4)-SideSlabIntervals(3))
      dz=ABS(SideSlabIntervals(6)-SideSlabIntervals(5))
      SideID = SideInfo
      SideBoundingBoxVolume(SideID)=dx*dy*dz
    END DO
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

! BezierAreaSample stuff:
WRITE(tmpStr,'(L1)') (TrackingMethod.EQ.TRIATRACKING)
TriaSurfaceFlux = GETLOGICAL('TriaSurfaceFlux',TRIM(tmpStr))
IF (Symmetry2D) THEN
  SWRITE(UNIT_stdOut,'(A)') "Surface Flux set to Triangle-aproximation due to Symmetry2D."
  TriaSurfaceFlux = .TRUE.
END IF
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
!> calculates XCL_Ngeo and dXCL_Ngeo for compute node mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc                ,ONLY: N
USE MOD_Basis                  ,ONLY: BuildBezierVdm,BuildBezierDMat
USE MOD_Basis                  ,ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Interpolation          ,ONLY: GetDerivativeMatrix
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars     ,ONLY: NodeType,NodeTypeCL,NodeTypeVISU
USE MOD_Mesh_Vars              ,ONLY: useCurveds
USE MOD_Mesh_Vars              ,ONLY: Elem_xGP
USE MOD_Particle_Mesh_Tools    ,ONLY: GetCNElemID
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeCL
USE MOD_Mesh_Vars              ,ONLY: NGeo,XCL_NGeo,wBaryCL_NGeo,XiCL_NGeo,dXCL_NGeo,InterpolateFromTree,Xi_NGeo
USE MOD_Mesh_Vars              ,ONLY: wBaryCL_NGeo1,Vdm_CLNGeo1_CLNGeo,XiCL_NGeo1,nElems,offsetElem
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: Vdm_Bezier,sVdm_Bezier,D_Bezier
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID
REAL                           :: Vdm_NGeo_CLNGeo(0:NGeo,0:NGeo)
REAL                           :: Vdm_EQNGeo_CLN (0:PP_N ,0:NGeo)
REAL                           :: Vdm_CLNloc_N   (0:PP_N ,0:PP_N)
REAL                           :: DCL_NGeo(0:Ngeo,0:Ngeo)
!INTEGER                        :: firstElem,lastElem
INTEGER                        :: firstHaloElem,lastHaloElem,nComputeNodeHaloElems
INTEGER                        :: CornerNodeIDswitch(8)
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: firstNodeID,nodeID,i,j,k,ll
!INTEGER                        :: nNodeIDs
REAL                           :: NodeCoordstmp(1:3,0:NGeo,0:NGeo,0:NGeo)
!===================================================================================================================================
! small wBaryCL_NGEO
ALLOCATE(wBaryCL_NGeo1(0:1),&
         XiCL_NGeo1   (0:1))
CALL ChebyGaussLobNodesAndWeights(1,XiCL_NGeo1)
CALL BarycentricWeights(1,XiCL_NGeo1,wBaryCL_NGeo1)
ALLOCATE(Vdm_CLNGeo1_CLNGeo(0:NGeo,0:1) )
CALL InitializeVandermonde(1, NGeo,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo ,Vdm_CLNGeo1_CLNGeo)
! new for curved particle sides
ALLOCATE(Vdm_Bezier(0:NGeo,0:NGeo),sVdm_Bezier(0:NGeo,0:NGeo))
! initialize vandermonde for bezier basis surface representation (particle tracking with curved elements)
CALL BuildBezierVdm(NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier) !CHANGETAG
ALLOCATE(D_Bezier(0:NGeo,0:NGeo))
CALL BuildBezierDMat(NGeo,Xi_NGeo,D_Bezier)

CALL InitializeVandermonde(NGeo,NGeo,wBaryCL_NGeo,Xi_NGeo,XiCL_NGeo,Vdm_NGeo_CLNGeo)

#if USE_MPI
! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
!CornerNodeIDswitch(1)=1
!CornerNodeIDswitch(2)=(Ngeo+1)
!CornerNodeIDswitch(3)=(Ngeo+1)**2
!CornerNodeIDswitch(4)=(Ngeo+1)*Ngeo+1
!CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
!CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
!CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2
!CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1
CornerNodeIDswitch(1)=1
CornerNodeIDswitch(2)=2
CornerNodeIDswitch(3)=3
CornerNodeIDswitch(4)=4
CornerNodeIDswitch(5)=5
CornerNodeIDswitch(6)=6
CornerNodeIDswitch(7)=7
CornerNodeIDswitch(8)=8

ASSOCIATE(CNS => CornerNodeIDswitch )

! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
MPISharedSize = INT((3*(NGeo+1)**3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (NGeo+1)*(NGeo+1)*(NGeo+1)*nComputeNodeTotalElems/), XCL_NGeo_Shared_Win,XCL_NGeo_Array)
MPISharedSize = INT((3*(PP_N+1)**3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems/), Elem_xGP_Shared_Win,Elem_xGP_Array)
MPISharedSize = INT((3*3*(NGeo+1)**3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nComputeNodeTotalElems/),dXCL_NGeo_Shared_Win,dXCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,dXCL_NGeo_Shared_Win,IERROR)
XCL_NGeo_Shared (1:3    ,0:NGeo,0:NGeo,0:NGeo,1:nComputeNodeTotalElems) => XCL_NGeo_Array
Elem_xGP_Shared (1:3    ,0:PP_N,0:PP_N,0:PP_N,1:nComputeNodeTotalElems) => Elem_xGP_Array
dXCL_NGeo_Shared(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nComputeNodeTotalElems) => dXCL_NGeo_Array
! Copy local XCL and dXCL into shared
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  DO iElem = 1, nElems
    XCL_NGeo_Shared (:  ,:,:,:,offsetElem+iElem) = XCL_NGeo (:  ,:,:,:,iElem)
    Elem_xGP_Shared (:  ,:,:,:,offsetElem+iElem) = Elem_xGP (:  ,:,:,:,iElem)
    dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
  END DO ! iElem = 1, nElems
ELSE
  DO iElem = 1, nElems
    XCL_NGeo_Shared (:,  :,:,:,GetCNElemID(offsetElem+iElem)) = XCL_NGeo (:,  :,:,:,iElem)
    Elem_xGP_Shared (:,  :,:,:,GetCNElemID(offsetElem+iElem)) = Elem_xGP (:,  :,:,:,iElem)
    dXCL_NGeo_Shared(:,:,:,:,:,GetCNElemID(offsetElem+iElem)) = dXCL_NGeo(:,:,:,:,:,iElem)
  END DO ! iElem = 1, nElems
END IF
nComputeNodeHaloElems = nComputeNodeTotalElems - nComputeNodeElems
IF (nComputeNodeHaloElems.GT.nComputeNodeProcessors) THEN
  firstHaloElem = INT(REAL( myComputeNodeRank   *nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))+1
  lastHaloElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))
ELSE
  firstHaloElem = myComputeNodeRank + 1
  IF (myComputeNodeRank.LT.nComputeNodeHaloElems) THEN
    lastHaloElem = myComputeNodeRank + 1
  ELSE
    lastHaloElem = 0
  END IF
END IF

! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!       Important for curved meshes if NGeo<N, no effect for N>=NGeo
CALL GetVandermonde(    NGeo, NodeTypeVISU, PP_N, NodeTypeCL, Vdm_EQNGeo_CLN,  modal=.FALSE.)
CALL GetVandermonde(    PP_N, NodeTypeCL  , PP_N, NodeType  , Vdm_CLNloc_N,    modal=.FALSE.)

!Transform from EQUI_NGeo to solution points on Nloc
Vdm_EQNGeo_CLN = MATMUL(Vdm_CLNloc_N,Vdm_EQNGeo_CLN)

! Build XCL and dXCL for compute node halo region (each proc of compute-node build only its fair share)
IF(interpolateFromTree) THEN
  CALL abort(__STAMP__,'ERROR: InterpolateFromTree not yet implemented for new halo region!')
ELSE
  CALL GetDerivativeMatrix(NGeo  , NodeTypeCL  , DCL_Ngeo)

  DO iElem = firstHaloElem, lastHaloElem
    ElemID = GetGlobalElemID(nComputeNodeElems+iElem)
!    firstNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+1
    firstNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)
!    nodeID = 0
    nodeID = 1
    IF (useCurveds) THEN
      DO k = 0, NGeo; DO j = 0, NGeo; DO i = 0, NGeo
        NodeCoordstmp(:,i,j,k) = NodeCoords_Shared(:,firstNodeID+NodeID)
        nodeID = nodeID + 1
      END DO; END DO; END DO ! i = 0, NGeo
    ELSE
      NodeCoordstmp(:,0,0,0) = NodeCoords_Shared(:,firstNodeID+CNS(1))
      NodeCoordstmp(:,1,0,0) = NodeCoords_Shared(:,firstNodeID+CNS(2))
      NodeCoordstmp(:,0,1,0) = NodeCoords_Shared(:,firstNodeID+CNS(3))
      NodeCoordstmp(:,1,1,0) = NodeCoords_Shared(:,firstNodeID+CNS(4))
      NodeCoordstmp(:,0,0,1) = NodeCoords_Shared(:,firstNodeID+CNS(5))
      NodeCoordstmp(:,1,0,1) = NodeCoords_Shared(:,firstNodeID+CNS(6))
      NodeCoordstmp(:,0,1,1) = NodeCoords_Shared(:,firstNodeID+CNS(7))
      NodeCoordstmp(:,1,1,1) = NodeCoords_Shared(:,firstNodeID+CNS(8))
!      DO i = 0, NGeo
!        DO j = 0, NGeo
!          DO k = 0, NGeo
!            NodeCoordstmp(:,i,j,k) = NodeCoords_Shared(:,firstNodeID+CNS(NodeID))
!            nodeID = nodeID + 1
!          END DO
!        END DO
!      END DO ! i = 0, NGeo
    END IF
    CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_NGeo_CLNGeo,NodeCoordstmp,XCL_NGeo_Shared(:,:,:,:,nComputeNodeElems+iElem))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_EQNGeo_CLN ,NodeCoordstmp,Elem_xGP_Shared(:,:,:,:,nComputeNodeElems+iElem))

    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      ! Matrix-vector multiplication
      DO ll=0,Ngeo
        dXCL_NGeo_Shared(1,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(1,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(i,ll) &
                                                            *  XCL_NGeo_Shared(: ,ll,j,k,nComputeNodeElems+iElem)
        dXCL_NGeo_Shared(2,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(2,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(j,ll) &
                                                            *  XCL_NGeo_Shared(: ,i,ll,k,nComputeNodeElems+iElem)
        dXCL_NGeo_Shared(3,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(3,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(k,ll) &
                                                            *  XCL_NGeo_Shared(: ,i,j,ll,nComputeNodeElems+iElem)
      END DO
    END DO; END DO; END DO !i,j,k=0,Ngeo
    END DO ! iElem = firstHaloElem, lastHaloElem
END IF

END ASSOCIATE

CALL MPI_WIN_SYNC(XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(XCL_NGeo_Shared (3,  0:NGeo,0:NGeo,0:NGeo,nElems))
ALLOCATE(dXCL_NGeo_Shared(3,3,0:NGeo,0:NGeo,0:NGeo,nElems))
XCL_NGeo_Shared  = XCL_NGeo
dXCL_NGeo_Shared = dXCL_NGeo
#endif /*USE_MPI*/

END SUBROUTINE CalcParticleMeshMetrics


SUBROUTINE CalcBezierControlPoints()
!===================================================================================================================================
!> calculate the bezier control point (+elevated) for shared compute-node mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Mappings               ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID, GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: GetBezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3DElevated,BezierElevation
#if USE_MPI
USE MOD_Particle_Mesh_Vars       ,ONLY: BezierControlPoints3D_Shared,BezierControlPoints3D_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: BezierControlPoints3DElevated_Shared,BezierControlPoints3DElevated_Shared_Win
USE MOD_MPI_Shared               ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars          ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars          ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars          ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
#endif /*USE_MPI*/
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iSide,ilocSide
INTEGER                        :: firstElem,lastElem,firstSide,lastSide
INTEGER                        :: p,q,pq(2), SideID
REAL                           :: tmp(3,0:NGeo,0:NGeo)
REAL                           :: tmp2(3,0:NGeo,0:NGeo)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER                        :: ALLOCSTAT
#endif /*USE_MPI*/
!===================================================================================================================================

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
  CALL Allocate_Shared(MPISharedSize,(/3*(NGeoElevated+1)*(NGeoElevated+1)*nNonUniqueGlobalSides/), &
                                      BezierControlPoints3DElevated_Shared_Win,BezierControlPoints3DElevated_Shared)
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
#endif

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
IF (BezierElevation.GT.0) THEN
  CALL MPI_WIN_SYNC(BezierControlPoints3DElevated_Shared_Win,IERROR)
END IF
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
!firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
!lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
firstSide = INT(REAL (myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif

DO iElem = firstElem, lastElem
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
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)
    DO q=0,NGeo; DO p=0,NGeo
      ! turn into right hand system of side
      pq=CGNS_SideToVol2(NGeo,p,q,iLocSide)
      BezierControlPoints3D(1:3,p,q,SideID)=tmp2(1:3,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! ilocSide=1,6
END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

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
#endif
END IF

END SUBROUTINE CalcBezierControlPoints


SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
! Subroutine for particle geometry initialization (GEO container)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
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
!INTEGER            :: ALLOCSTAT
INTEGER            :: iElem,ElemID,iLocSide,iNode
!INTEGER            :: jNode
INTEGER            :: nStart, NodeNum
INTEGER            :: NodeMap(4,6)
!INTEGER            :: nSides
REAL               :: A(3,3),detcon
!REAL,ALLOCATABLE   :: Coords(:,:,:,:)
!CHARACTER(32)      :: tmpStr
!CHARACTER(LEN=255) :: FileString
INTEGER            :: FirstElem, LastElem, GlobalSideID, nlocSides, localSideID
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

!ALLOCATE(GEO%ElemToNodeID(1:8,1:nElems),       &
!         GEO%ElemSideNodeID(1:4,1:6,1:nElems), &
!         GEO%NodeCoords(1:3,1:nNodes),         &
!         GEO%ConcaveElemSide(1:6,1:nElems), &
!         GEO%ElemMidPoint(1:3,nElems), STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) THEN
! CALL abort(__STAMP__&
! ,'ERROR in InitParticleGeometry: Cannot allocate GEO%... stuff!')
!END IF

#if USE_MPI
MPISharedSize = INT(6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),ConcaveElemSide_Shared_Win,ConcaveElemSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,ConcaveElemSide_Shared_Win,IERROR)
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

MPISharedSize = INT(8*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/8,nComputeNodeTotalElems/),ElemNodeID_Shared_Win,ElemNodeID_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemNodeID_Shared_Win,IERROR)

MPISharedSize = INT(4*6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/4,6,nComputeNodeTotalElems/),ElemSideNodeID_Shared_Win,ElemSideNodeID_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemSideNodeID_Shared_Win,IERROR)

MPISharedSize = INT(3*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemMidPoint_Shared_Win,ElemMidPoint_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemMidPoint_Shared_Win,IERROR)
#else
ALLOCATE(ConcaveElemSide_Shared(   1:6,1:nElems))
ALLOCATE(ElemNodeID_Shared(        1:8,1:nElems))
ALLOCATE(ElemSideNodeID_Shared(1:4,1:6,1:nElems))
ALLOCATE(ElemMidPoint_Shared(      1:3,1:nElems))
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

!ALLOCATE(GEO%ElemToNodeIDGlobal(1:8,1:nElems))

!GEO%ElemToNodeID(:,:)=0
!GEO%ElemSideNodeID(:,:,:)=0
!GEO%NodeCoords(:,:)=0.
!GEO%ConcaveElemSide(:,:)=.FALSE.
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif
  ElemNodeID_Shared      = 0
  ElemSideNodeID_Shared  = 0
  ConcaveElemSide_Shared = .FALSE.
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(ElemNodeID_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemSideNodeID_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ConcaveElemSide_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

!iNode=0
!DO iElem=1,nElems
!  DO jNode=1,8
!    Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=0
!  END DO
!END DO
!DO iElem=1,nElems
!  !--- Save corners of sides
!  DO jNode=1,8
!    IF (Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID.EQ.0) THEN
!      iNode=iNode+1
!      Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=iNode
!      GEO%NodeCoords(1:3,iNode)=Elems(iElem+offsetElem)%ep%node(jNode)%np%x(1:3)
!    END IF
!    GEO%ElemToNodeID(jNode,iElem)=Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID
!    !GEO%ElemToNodeIDGlobal(jNode,iElem) = Elems(iElem+offsetElem)%ep%node(jNode)%np%ind
!  END DO
!END DO
!
!DO iElem=1,nElems
!  DO iLocSide=1,6
!    nStart=MAX(0,ElemToSide(E2S_FLIP,iLocSide,iElem)-1)
!    GEO%ElemSideNodeID(1:4,iLocSide,iElem)=(/Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart  ,4)+1,iLocSide))%np%NodeID,&
!                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+1,4)+1,iLocSide))%np%NodeID,&
!                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+2,4)+1,iLocSide))%np%NodeID,&
!                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+3,4)+1,iLocSide))%np%NodeID/)
!  END DO
!END DO

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  DO iNode = 1,8
    ElemNodeID_Shared(iNode,iElem) = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID) + CNS(iNode)
  END DO

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  DO iLocSide = 1,nlocSides
    ! Get global SideID
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    ! Find start of CGNS mapping from flip
    nStart = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)-1)
    ! Shared memory array starts at 1, but NodeID at 0
    ElemSideNodeID_Shared(1:4,localSideID,iElem) = (/ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+NodeMap(MOD(nStart  ,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+NodeMap(MOD(nStart+1,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+NodeMap(MOD(nStart+2,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+NodeMap(MOD(nStart+3,4)+1,localSideID)-1/)
  END DO
END DO
END ASSOCIATE
#if USE_MPI
CALL MPI_WIN_SYNC(ElemNodeID_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemSideNodeID_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif
!--- Save whether Side is concave or convex
DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  DO iLocSide = 1,nlocSides
    !--- Check whether the bilinear side is concave
    !--- Node Number 4 and triangle 1-2-3
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    DO NodeNum = 1,3               ! for all 3 nodes of triangle
!      A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,iElem)) &
!                   - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,iLocSide,iElem))
       A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,iElem)+1) &
                    - NodeCoords_Shared(:,ElemSideNodeID_Shared(4      ,localSideID,iElem)+1)
    END DO
    !--- concave if detcon < 0:
    detcon = ((A(2,1) * A(3,2) - A(3,1) * A(2,2)) * A(1,3) +     &
              (A(3,1) * A(1,2) - A(1,1) * A(3,2)) * A(2,3) +     &
              (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * A(3,3))
!    IF (detcon.LT.0) GEO%ConcaveElemSide(iLocSide,iElem)=.TRUE.
    IF (detcon.LT.0) ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
  END DO
END DO
!SWRITE(*,*) ConcaveElemSide_Shared(:,25)
!SideInfo_Shared(SIDE_NBELEMID,ElemInfo_Shared(ELEM_FIRSTSIDEIND,1)+1)

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

!DO iElem =1, nElems
!  GEO%ElemMidPoint(:,iElem) = 0.0
!  DO iNode = 1,8
!    GEO%ElemMidPoint(1:3,iElem) = GEO%ElemMidPoint(1:3,iElem) + GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
!  END DO
!  GEO%ElemMidPoint(1:3,iElem) = GEO%ElemMidPoint(1:3,iElem) / 8.
!END DO

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ElemMidPoint_Shared(:,iElem) = 0.
  DO iNode = 1,8
    ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode)
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


SUBROUTINE WriteTriaDataToVTK(nSides,nElems,Coord,FileString)
!===================================================================================================================================
!> Routine writing data to VTK Triangles (cell type = 5)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nSides               !< Number of sides per element
INTEGER,INTENT(IN)          :: nElems               !< Number of elements
REAL   ,INTENT(IN)          :: Coord(1:3,1:4,1:nSides,1:nElems)
CHARACTER(LEN=*),INTENT(IN) :: FileString           ! < Output file name
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER            :: iElem,nVTKElems,nVTKCells,ivtk=44,iSide!,iVal,iVar,str_len
INTEGER(KIND=8)    :: Offset, nBytes
INTEGER            :: IntegerType
INTEGER            :: Vertex(3,nSides*nElems*2)
INTEGER            :: NodeID,CellID,CellType
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf!,components_string
!CHARACTER(LEN=255) :: VarNameString
REAL(KIND=4)       :: FloatType
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE TRIA DATA TO VTX XML BINARY (VTU) FILE..."
IF(nSides.LT.1)THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
  RETURN
END IF

! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
nVTKElems=nSides*nElems*4 ! number of Nodes
nVTKCells=nSides*2*nElems ! number of Triangles

Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKElems
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//&
'" NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
!IF (nVal .GT.0)THEN
!  DO iVar=1,nVal
!    IF (VarNamePartCombine(iVar).EQ.0) THEN
!      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNamePartVisu(iVar))//&
!      '" NumberOfComponents="1" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
!      Offset=Offset+SIZEOF(IntegerType)+nVTKElems*SIZEOF(FloatType)
!      WRITE(StrOffset,'(I16)')Offset
!    ELSE IF (VarNamePartCombine(iVar).EQ.1) THEN
!      str_len = LEN_TRIM(VarNamePartVisu(iVar))
!      write(components_string,'(I1)') VarNamePartCombineLen(iVar)
!      !IF(FileType.EQ.'DSMCHOState')THEN
!      !  VarNameString = VarNamePartVisu(iVar)(1:str_len-4)//VarNamePartVisu(iVar)(str_len-2:str_len)
!      !ELSE
!        VarNameString = VarNamePartVisu(iVar)(1:str_len-1)
!      !END IF
!      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameString)//&
!      '" NumberOfComponents="'//components_string//'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf
!      WRITE(ivtk) TRIM(Buffer)
!      Offset=Offset+SIZEOF(IntegerType)+nVTKElems*SIZEOF(FloatType)*VarNamePartCombineLen(iVar)
!      WRITE(StrOffset,'(I16)')Offset
!    END IF
!  END DO
!END IF
Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(IntegerType)+3*nVTKElems*SIZEOF(FloatType)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(IntegerType)+nVTKCells*3*SIZEOF(IntegerType)
WRITE(StrOffset,'(I16)')Offset
! Offsets
Buffer='        <DataArray type="Int32" Name="offsets" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(IntegerType)+nVTKCells*SIZEOF(IntegerType)
WRITE(StrOffset,'(I16)')Offset
! Elem types
Buffer='        <DataArray type="Int32" Name="types" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
! Prepare append section
Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
! Write leading data underscore
Buffer='_';WRITE(ivtk) TRIM(Buffer)

! Write binary raw data into append section
! Point data
nBytes = nVTKElems*SIZEOF(FloatType)
!DO iVal=1,nVal
!  IF (VarNamePartCombine(iVal).EQ.0) THEN
!    WRITE(ivtk) nBytes,REAL(Value(1:nParts,iVal),4)
!  ELSEIF(VarNamePartCombine(iVal).EQ.1) THEN
!    WRITE(ivtk) nBytes*VarNamePartCombineLen(iVal),REAL(Value(1:nParts,iVal:iVal+VarNamePartCombineLen(iVal)-1),4)
!  ENDIF
!END DO
! Points
nBytes = nBytes * 3
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coord,4)
! Connectivity
NodeID = -1
CellID = 1
DO iElem=1,nElems
  DO iSide=1,6
    ! nodes 1,2,3 and nodes 1,3,4 forming one triangle
    ! nodes indexes start with 0 in vtk
    Vertex(:,CellID) = (/NodeID+1,NodeID+2,NodeID+3/)
    Vertex(:,CellID+1) = (/NodeID+1,NodeID+3,NodeID+4/)
    CellID=CellID+2
    NodeID=NodeID+4
  END DO
END DO
nBytes = 3*nVTKCells*SIZEOF(IntegerType)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset
nBytes = nVTKCells*SIZEOF(IntegerType)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=3,3*nVTKCells,3)
! Cell type
CellType = 5  ! VTK_TRIANGLE
!CellType = 6  ! VTK_TRIANGLE_STRIP
WRITE(ivtk) nBytes
WRITE(ivtk) (CellType,iElem=1,nVTKCells)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"

END SUBROUTINE WriteTriaDataToVTK


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
USE MOD_Particle_Mesh_Tools       ,ONLY: GetGlobalElemID
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
DO iElem = firstElem,lastElem ! go through all elements
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

SWRITE(UNIT_stdOut,'(A)')' CHECKING FOR WEIRD ELEMENTS DONE!'
IF(WeirdElems.GT.0) THEN
  IPWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
  IPWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
  DO iElem = 1,WeirdElems
    IPWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
  END DO
END IF

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
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID, GetGlobalNonUniqueSideID
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadiusNGEO_Shared,ElemRadiusNGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo_Shared,ElemRadius2NGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
#else
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo
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
!INTEGER                        :: ALLOCSTAT
INTEGER                        :: iElem,SideID
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
#endif

#if USE_MPI
IF(myComputeNodeRank.EQ.0) THEN
#endif
  ElemRadiusNGeo =0.
  ElemRadius2NGeo=0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

#if USE_MPI
firstElem=INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

DO iElem=firstElem,lastElem
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
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

IF (TRIM(DepositionType).EQ.'shape_function_simple')THEN
  CALL abort(&
  __STAMP__&
  ,'not implemented yet')
  !ALLOCATE(ElemRadius2_sf(1:nElems))
  !DO iElem=1,nElems
  !  ElemRadius2_sf(iElem)=(ElemRadiusNGeo(iElem)+r_sf)*(ElemRadiusNGeo(iElem)+r_sf)
  !END DO ! iElem=1,nElems
END IF

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(XiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(slenXiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

END SUBROUTINE BuildElementBasisAndRadius


SUBROUTINE BuildElementRadiusTria()
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo_Shared,ElemRadius2NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared_Win,ElemRadius2NGeo_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: nELems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                        :: ALLOCSTAT
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
#else
ALLOCATE(ElemBaryNGeo(1:3,nElems) &
        ,ElemRadius2NGeo( nElems))
#endif

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemRadius2NGeo = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

#if USE_MPI
firstElem=INT(REAL(myComputeNodeRank*    nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif

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
END DO ! iElem

#if USE_MPI
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

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
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Mesh_Vars               ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
USE MOD_MPI_Shared              ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo
USE MOD_Mesh_Vars               ,ONLY: nElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iDir
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
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos = 0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos = xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
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
USE MOD_Particle_Mesh_Vars ,ONLY: ElemBaryNGeo
#if USE_MPI
USE MOD_MPI_Shared        ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars   ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars   ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars   ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars, ONLY: ElemBaryNGeo_Shared,ElemBaryNGeo_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Mesh_Vars ,ONLY: nComputeNodeElems
USE MOD_Mesh_Vars          ,ONLY: XCL_NGeo
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,i,j,k
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

DO iElem=firstElem,lastElem
  ! evaluate the polynomial at origin: Xi=(/0.0,0.0,0.0/)
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  xPos=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buf=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        xPos=xPos+XCL_NGeo(1:3,i,j,k,iElem)*Lag(1,i)*buf
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
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID,GetCNElemID, GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,SideType,SideNormVec,SideDistance
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideDistance_Shared,SideDistance_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideType_Shared,SideType_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideNormVec_Shared,SideNormVec_Shared_Win
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
!IF (.NOT.DoRefMapping) THEN
!  ALLOCATE(ElemType(1:nTotalElems))
!  ElemType=-1
!END IF

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
  XCL_NGeoLoc = XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
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
  ! check the coordinates of all Chebychev-Lobatto geometry points between the bi-linear and used
  ! mapping
  CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo),ElemCurved(iElem))

  ! 2) check sides
  ! loop over all 6 sides of element
  ! a) check if the sides are straight
  ! b) use curved information to decide side type
  DO ilocSide=1,6
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)

    ! Why were only flips LT. 0 considered? All flips are equal!
!    IF (SideInfo_Shared(SIDE_ID,SideID).GT.0) THEN
!      flip = 0
!    ELSE
!      flip = MOD(Sideinfo_Shared(SIDE_FLIP,SideID),10)
!    END IF
    flip = Sideinfo_Shared(SIDE_FLIP,SideID)

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

! TODO with PERIODIC SIDES
! sanity check for side periodic type
!DO iSide=1,nPartSides
!  IF(DoRefmapping)THEN
!    BCSideID  =PartBCSideList(iSide)
!    IF(BCSideID.LE.0) CYCLE
!  ELSE
!    BCSideID  =iSide
!  END IF
!  PVID=SidePeriodicType(iSide)
!  IF(PVID.EQ.0) CYCLE
!  IF(.NOT.PartMeshHasPeriodicBCs) CYCLE
!  Vec1=SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
!  ScalarProduct=DOT_PRODUCT(SideNormVec(1:3,BCSideID),Vec1)
!  IF(ALMOSTEQUAL(ScalarProduct,GEO%PeriodicVectorsLength(ABS(PVID))))THEN
!    SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!  ELSEIF(.NOT.ALMOSTEQUAL(ScalarProduct,-GEO%PeriodicVectorsLength(ABS(PVID))))THEN
!    WRITE (*,*) "BCSideID                  : ", BCSideID
!    WRITE (*,*) "SideNormVec(1:3,BCSideID) : ", SideNormVec(1:3,BCSideID)
!    WRITE (*,*) "Vec1                      : ", Vec1
!    CALL abort(&
!__STAMP__&
!        , ' Missalignment between SideNormVec and PeriodicVector!',ABS(PVID),ScalarProduct)
!  END IF
!END DO ! iSide=1,nPartSides

!! fill Element type checking sides
!IF (.NOT.DoRefMapping) THEN
!  DO iElem=1,nTotalElems
!    DO ilocSide=1,6
!      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!      SELECT CASE(SideType(SideID))
!      CASE(PLANAR_RECT,PLANAR_NONRECT)
!        IF (ElemType(iElem).GE.1) THEN
!          CYCLE
!        ELSE
!          ElemType(iElem) = 1
!        END IF
!      CASE(BILINEAR)
!        IF (ElemType(iElem).GE.2) THEN
!          CYCLE
!        ELSE
!          ElemType(iElem) = 2
!        END IF
!      CASE(PLANAR_CURVED,CURVED)
!        ElemType(iElem) = 3
!        EXIT
!      END SELECT
!    END DO ! ilocSide=1,6
!  END DO ! iElem=1,nTotalElems
!END IF

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
! This only works if the full mesh is built on one node. Make it failproof using global communication
!IF (myComputeNodeRank.EQ.0) THEN
!  CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nLinearElems         ,nLinearElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nCurvedElems         ,nCurvedElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!ELSE
!  CALL MPI_REDUCE(nPlanarRectangular   ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nPlanarNonRectangular,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nBilinear            ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nPlanarCurved        ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nCurved              ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nLinearElems         ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!  CALL MPI_REDUCE(nCurvedElems         ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
!END IF
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
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,nUniqueBCSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID_Shared,SideID2BCSide_Shared,BCSideMetrics_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID_Shared_Win,SideID2BCSide_Shared_Win,BCSideMetrics_Shared_Win
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
INTEGER                        :: iSide,firstSide,lastSide,BCSideID
INTEGER                        :: nUniqueBCSidesProc,offsetUniqueBCSidesProc
REAL                           :: origin(1:3),xi(1:2),radius,radiusMax,vec(1:3)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! Count number of BC sides in range
nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
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
ALLOCATE(BCSideMetrics(1:4,1:nUniqueBCSides)
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
  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
  BCSideID           = offsetUniqueBCSidesProc + nUniqueBCSidesProc
  BCSide2SideID(BCSideID) = iSide
  SideID2BCSide(iSide)    = BCSideID

  ! calculate origin, radius for all BC sides
  !> build side origin
  xi     = 0.
  ! TODO: BezierControlPoints are alloced with global side ID, so this SHOULD work. Breaks if we reduce the halo region
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
USE MOD_Equation_Vars          ,ONLY: c
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides,SideBCMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo,ElemRadiusNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,nUniqueBCSides
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID, GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Vars          ,ONLY: ManualTimeStep
USE MOD_Utils                  ,ONLY: InsertionSort
#if USE_MPI
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
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
USE MOD_Particle_Mesh_Vars     ,ONLY: nComputeNodeElems,nComputeNodeSides
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
INTEGER,ALLOCATABLE            :: ElemToBCSidesProc(:,:)
REAL                           :: dX,dY,dZ
REAL                           :: origin(1:3),vec(1:3)
REAL                           :: BC_halo_eps,BC_halo_eps_velo,BC_halo_diag,deltaT
LOGICAL                        :: fullMesh
REAL,ALLOCATABLE               :: tmpSideBCMetrics(:,:)
REAL,ALLOCATABLE               :: tmpSideBCDistance(:)
INTEGER,ALLOCATABLE            :: intSideBCMetrics(:)
#if USE_MPI
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
      ElemToBCSidesProc(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides + 1
    END IF
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

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

    ! loop over all sides. Check distance from every local side to total sides.
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,ElemID) - ElemBaryNGeo(:,BCElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(ElemID) + ElemRadiusNGeo(BCElemID))) CYCLE

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
      ElemToBCSidesProc(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides + 1
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

ElemToBCSides(1,firstElem:lastElem) = ElemToBCSidesProc(1,firstElem:lastElem)
ElemToBCSides(2,firstElem:lastElem) = ElemToBCSidesProc(2,firstElem:lastElem) + offsetBCSidesProc
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
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(iElem)
    END DO

    ! loop over all sides. Check distance from every local side to total sides. Once a side has been flagged,
    ! it must not be counted again
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,ElemID) - ElemBaryNGeo(:,BCElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(ElemID) + ElemRadiusNGeo(BCElemID))) CYCLE

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
  ElemID   = INT(SideBCMetrics(BCSIDE_ELEMID,iSide))

  !> get origin and radius from BC Side
  SideBCMetrics(5:7          ,iSide) = BCSideMetrics(1:3,BCSideID)
  SideBCMetrics(BCSIDE_RADIUS,iSide) = BCSideMetrics(4  ,BCSideID)

  !> build side distance
  origin(1:3) = ElemBaryNGeo(1:3,ElemID)
  vec(1:3)    = origin(1:3) - SideBCMetrics(5:7,iSide)
  SideBCMetrics(BCSIDE_DISTANCE,iSide) = SQRT(DOT_PRODUCT(vec,vec))-ElemRadiusNGeo(ElemID)-SideBCMetrics(BCSIDE_RADIUS,iSide)
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
  firstSide    = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem)
  lastSide     = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem)+ElemToBCSides(ELEM_NBR_BCSIDES,iElem)-1
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
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoRef
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
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
ALLOCATE(ElemsJ(0:PP_N,0:PP_N,0:PP_N,1:nElems))
ElemsJ => sJ
#endif /* USE_MPI*/

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
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: NodeCoords_Shared,ElemBaryNGeo_Shared
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
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nElems
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,NodeCoords_Shared
#if USE_MPI
USE MOD_Particle_Mesh_Vars ,ONLY: offsetComputeNodeNode,nComputeNodeNodes
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: offsetLocalNode,nLocalNodes
!===================================================================================================================================
! calculate all offsets
offsetLocalNode = ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)
nLocalNodes     = ElemInfo_Shared(ELEM_LASTNODEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)

GEO%xmin     = MINVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
GEO%xmax     = MAXVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
GEO%ymin     = MINVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
GEO%ymax     = MAXVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
GEO%zmin     = MINVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
GEO%zmax     = MAXVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))

#if USE_MPI
GEO%CNxmin   = MINVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNxmax   = MAXVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNymin   = MINVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNymax   = MAXVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNzmin   = MINVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNzmax   = MAXVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%xminglob = MINVAL(NodeCoords_Shared(1,:))
GEO%xmaxglob = MAXVAL(NodeCoords_Shared(1,:))
GEO%yminglob = MINVAL(NodeCoords_Shared(2,:))
GEO%ymaxglob = MAXVAL(NodeCoords_Shared(2,:))
GEO%zminglob = MINVAL(NodeCoords_Shared(3,:))
GEO%zmaxglob = MAXVAL(NodeCoords_Shared(3,:))
#else
GEO%CNxmin   = GEO%xmin
GEO%CNxmax   = GEO%xmax
GEO%CNymin   = GEO%ymin
GEO%CNymax   = GEO%ymax
GEO%CNzmin   = GEO%zmin
GEO%CNzmax   = GEO%zmax
GEO%xminglob = GEO%xmin
GEO%xmaxglob = GEO%xmax
GEO%yminglob = GEO%ymin
GEO%ymaxglob = GEO%ymax
GEO%zminglob = GEO%zmin
GEO%zmaxglob = GEO%zmax
#endif /*USE_MPI*/

END SUBROUTINE GetMeshMinMax


SUBROUTINE FinalizeParticleMesh()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nElems,nNodes
USE MOD_Particle_BGM           ,ONLY: FinalizeBGM
USE MOD_Particle_Mesh_Readin   ,ONLY: FinalizeMeshReadin
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod,Distance,ListDistance
#if USE_MPI
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
#endif
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
INTEGER                             :: iELem,iNode
!===================================================================================================================================

! Particle mesh readin happens during mesh readin, finalize with gathered routine here
CALL FinalizeMeshReadin()

CALL FinalizeBGM()

SELECT CASE (TrackingMethod)

  ! RefMapping, Tracing
  CASE(1,2)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! GetBCSidesAndOrgin
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL MPI_WIN_UNLOCK_ALL(BCSide2SideID_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      BCSide2SideID_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(SideID2BCSide_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      SideID2BCSide_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(BCSideMetrics_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      BCSideMetrics_Shared_Win        ,iError)
    END IF

    ! CalcParticleMeshMetrics
    CALL MPI_WIN_UNLOCK_ALL(XCL_NGeo_Shared_Win             ,iError)
    CALL MPI_WIN_FREE(      XCL_NGeo_Shared_Win             ,iError)
    CALL MPI_WIN_UNLOCK_ALL(Elem_xGP_Shared_Win             ,iError)
    CALL MPI_WIN_FREE(      Elem_xGP_Shared_Win             ,iError)
    CALL MPI_WIN_UNLOCK_ALL(dXCL_NGeo_Shared_Win            ,iError)
    CALL MPI_WIN_FREE(      dXCL_NGeo_Shared_Win            ,iError)

    ! CalcBezierControlPoints
    CALL MPI_WIN_UNLOCK_ALL(BezierControlPoints3D_Shared_Win,iError)
    CALL MPI_WIN_FREE(      BezierControlPoints3D_Shared_Win,iError)
    IF (BezierElevation.GT.0) THEN
      CALL MPI_WIN_UNLOCK_ALL(BezierControlPoints3DElevated_Shared_Win,iError)
      CALL MPI_WIN_FREE(      BezierControlPoints3DElevated_Shared_Win,iError)
    END IF

    ! GetSideSlabNormalsAndIntervals (allocated in particle_mesh.f90)
    CALL MPI_WIN_UNLOCK_ALL(SideSlabNormals_Shared_Win      ,iError)
    CALL MPI_WIN_FREE(      SideSlabNormals_Shared_Win      ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideSlabIntervals_Shared_Win    ,iError)
    CALL MPI_WIN_FREE(      SideSlabIntervals_Shared_Win    ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BoundingBoxIsEmpty_Shared_Win   ,iError)
    CALL MPI_WIN_FREE(      BoundingBoxIsEmpty_Shared_Win   ,iError)

    ! BuildElementOriginShared
    CALL MPI_WIN_UNLOCK_ALL(ElemBaryNGeo_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      ElemBaryNGeo_Shared_Win         ,iError)

    ! IdentifyElemAndSideType
    CALL MPI_WIN_UNLOCK_ALL(ElemCurved_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      ElemCurved_Shared_Win           ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideType_Shared_Win             ,iError)
    CALL MPI_WIN_FREE(      SideType_Shared_Win             ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideDistance_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      SideDistance_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideNormVec_Shared_Win          ,iError)
    CALL MPI_WIN_FREE(      SideNormVec_Shared_Win          ,iError)

    ! BuildElementBasisAndRadius
    CALL MPI_WIN_UNLOCK_ALL(ElemRadiusNGeo_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      ElemRadiusNGeo_Shared_Win       ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemRadius2NGeo_Shared_Win      ,iError)
    CALL MPI_WIN_FREE(      ElemRadius2NGeo_Shared_Win      ,iError)
    CALL MPI_WIN_UNLOCK_ALL(XiEtaZetaBasis_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      XiEtaZetaBasis_Shared_Win       ,iError)
    CALL MPI_WIN_UNLOCK_ALL(slenXiEtaZetaBasis_Shared_Win   ,iError)
    CALL MPI_WIN_FREE(      slenXiEtaZetaBasis_Shared_Win   ,iError)

    ! GetLinearSideBaseVectors
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors0_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors0_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors1_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors1_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors2_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors2_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors3_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors3_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectorsScale_Shared_Win     ,iError)
    CALL MPI_WIN_FREE(      BaseVectorsScale_Shared_Win     ,iError)

    ! BuildBCElemDistance
    IF (TrackingMethod.EQ.1) THEN
      CALL MPI_WIN_UNLOCK_ALL(ElemToBCSides_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      ElemToBCSides_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(SideBCMetrics_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      SideBCMetrics_Shared_Win        ,iError)
    END IF

    ! BuildEpsOneCell
    CALL MPI_WIN_UNLOCK_ALL(ElemsJ_Shared_Win               ,iError)
    CALL MPI_WIN_FREE(      ElemsJ_Shared_Win               ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemEpsOneCell_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      ElemEpsOneCell_Shared_Win       ,iError)

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

    ! CalcParticleMeshMetrics
    ADEALLOCATE(XCL_NGeo_Array)
    ADEALLOCATE(Elem_xGP_Array)
    ADEALLOCATE(dXCL_NGeo_Array)

    ! CalcBezierControlPoints
    ADEALLOCATE(BezierControlPoints3D_Shared)
    ADEALLOCATE(BezierControlPoints3DElevated_Shared)

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
!  CASE(2)

  ! TriaTracking
  CASE(3)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! CalcParticleMeshMetrics
    CALL MPI_WIN_UNLOCK_ALL(XCL_NGeo_Shared_Win             ,iError)
    CALL MPI_WIN_FREE(      XCL_NGeo_Shared_Win             ,iError)
    CALL MPI_WIN_UNLOCK_ALL(Elem_xGP_Shared_Win             ,iError)
    CALL MPI_WIN_FREE(      Elem_xGP_Shared_Win             ,iError)
    CALL MPI_WIN_UNLOCK_ALL(dXCL_NGeo_Shared_Win            ,iError)
    CALL MPI_WIN_FREE(      dXCL_NGeo_Shared_Win            ,iError)

    ! InitParticleGeometry
    CALL MPI_WIN_UNLOCK_ALL(ConcaveElemSide_Shared_Win,iError)
    CALL MPI_WIN_FREE(ConcaveElemSide_Shared_Win,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemNodeID_Shared_Win,iError)
    CALL MPI_WIN_FREE(ElemNodeID_Shared_Win,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemSideNodeID_Shared_Win,iError)
    CALL MPI_WIN_FREE(ElemSideNodeID_Shared_Win,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemMidPoint_Shared_Win,iError)
    CALL MPI_WIN_FREE(ElemMidPoint_Shared_Win,iError)

    ! BuildElementRadiusTria
    CALL MPI_WIN_UNLOCK_ALL(ElemBaryNGeo_Shared_Win,iError)
    CALL MPI_WIN_FREE(ElemBaryNGeo_Shared_Win,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemRadius2NGeo_Shared_Win,iError)
    CALL MPI_WIN_FREE(ElemRadius2NGeo_Shared_Win,iError)
    CALL MPI_WIN_UNLOCK_ALL(XiEtaZetaBasis_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      XiEtaZetaBasis_Shared_Win       ,iError)
    CALL MPI_WIN_UNLOCK_ALL(slenXiEtaZetaBasis_Shared_Win   ,iError)
    CALL MPI_WIN_FREE(      slenXiEtaZetaBasis_Shared_Win   ,iError)

    ! BuildElemTypeTria
    CALL MPI_WIN_UNLOCK_ALL(ElemCurved_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      ElemCurved_Shared_Win           ,iError)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
    ! CalcParticleMeshMetrics
    ADEALLOCATE(XCL_NGeo_Array)
    ADEALLOCATE(Elem_xGP_Array)
    ADEALLOCATE(dXCL_NGeo_Array)

    ! InitParticleGeometry
    ADEALLOCATE(ConcaveElemSide_Shared)
    ADEALLOCATE(ElemNodeID_Shared)
    ADEALLOCATE(ElemSideNodeID_Shared)
    ADEALLOCATE(ElemMidPoint_Shared)

    ! BuildElementRadiusTria
    ADEALLOCATE(ElemBaryNGeo_Shared)
    ADEALLOCATE(ElemRadius2NGEO_Shared)
    ADEALLOCATE(XiEtaZetaBasis_Shared)
    ADEALLOCATE(slenXiEtaZetaBasis_Shared)

    !  BuildElemTypeTria
    ADEALLOCATE(ElemCurved)
    ADEALLOCATE(ElemCurved_Shared)

END SELECT

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
SDEALLOCATE(ElemTolerance)
SDEALLOCATE(ElemToGlobalElemID)

! Old stuff
!!SDEALLOCATE(PartElemToSide)
!!SDEALLOCATE(PartSideToElem)
!!SDEALLOCATE(PartElemIsMortar)
!!SDEALLOCATE(PartElemToElemGlob)
!!SDEALLOCATE(PartElemToElemAndSide)
!SDEALLOCATE(PartBCSideList)
!SDEALLOCATE(SidePeriodicType)
!SDEALLOCATE(SidePeriodicDisplacement)
!!SDEALLOCATE(IsTracingBCElem)
!!SDEALLOCATE(TracingBCInnerSides)
!!SDEALLOCATE(TracingBCTotalSides)
!SDEALLOCATE(ElemType)
!SDEALLOCATE(GEO%PeriodicVectors)
!SDEALLOCATE(GEO%PeriodicVectorsLength)
!SDEALLOCATE(GEO%FIBGM)
!!SDEALLOCATE(GEO%Volume)
!!SDEALLOCATE(GEO%MPVolumePortion)
!!SDEALLOCATE(GEO%CharLength)
!SDEALLOCATE(GEO%ElemToFIBGM)
!SDEALLOCATE(GEO%TFIBGM)
!
!!SDEALLOCATE(GEO%ElemToNodeID)
!!SDEALLOCATE(GEO%ElemSideNodeID)
!!SDEALLOCATE(GEO%ElemToNodeIDGlobal)
!!SDEALLOCATE(GEO%NodeCoords)
!SDEALLOCATE(GEO%ElemsOnNode)
!SDEALLOCATE(GEO%NeighNodesOnNode)
!SDEALLOCATE(GEO%NumNeighborElems)
!IF (ALLOCATED(GEO%ElemToNeighElems)) THEN
!  DO iElem=1,nElems
!    SDEALLOCATE(GEO%ElemToNeighElems(iElem)%ElemID)
!  END DO
!END IF
!SDEALLOCATE(GEO%ElemToNeighElems)
!IF (ALLOCATED(GEO%NodeToElem)) THEN
!  DO iNode=1,nNodes
!    SDEALLOCATE(GEO%NodeToElem(iNode)%ElemID)
!  END DO
!END IF
!SDEALLOCATE(GEO%NodeToElem)
!IF (ALLOCATED(GEO%NodeToNeighNode)) THEN
!  DO iNode=1,nNodes
!    SDEALLOCATE(GEO%NodeToNeighNode(iNode)%ElemID)
!  END DO
!END IF
!SDEALLOCATE(GEO%NodeToNeighNode)
!!SDEALLOCATE(GEO%ConcaveElemSide)
!!SDEALLOCATE(GEO%ElemMidPoint)
!!SDEALLOCATE(GEO%BoundsOfElem)
!
!SDEALLOCATE(BCElem)
!ADEALLOCATE(XiEtaZetaBasis)
!ADEALLOCATE(slenXiEtaZetaBasis)
!ADEALLOCATE(ElemRadiusNGeo)
!ADEALLOCATE(ElemRadius2NGeo)
!ADEALLOCATE(ElemEpsOneCell)
!SDEALLOCATE(Distance)
!SDEALLOCATE(ListDistance)
!!SDEALLOCATE(isTracingTrouble)
!SDEALLOCATE(ElemTolerance)
!SDEALLOCATE(ElemToGlobalElemID)
!SDEALLOCATE(ElemHaloInfoProc)

!ADEALLOCATE(ElemCurved)

ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh


!SUBROUTINE InitFIBGM()
!!===================================================================================================================================
!! Build Fast-Init-Background-Mesh.
!! The BGM is a cartesian mesh for easier locating of particles
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_ReadInTools            ,ONLY: GetRealArray,GetLogical
!USE MOD_Mesh_Vars              ,ONLY: nElems
!USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,nTotalElems,nTotalBCSides, FindNeighbourElems
!USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
!#if USE_MPI
!USE MOD_Particle_MPI           ,ONLY: InitHALOMesh, AddHaloNodeData
!USE MOD_Particle_MPI_Vars      ,ONLY: printMPINeighborWarnings,printBezierControlPointsWarnings
!USE MOD_PICDepo_Vars           ,ONLY: CellLocNodes_Volumes, DepositionType
!#endif /*USE_MPI*/
!USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
!USE MOD_PICDepo_Vars           ,ONLY: ElemRadius2_sf,DepositionType,DoSFLocalDepoAtBounds
!USE MOD_Particle_Mesh_Tools    ,ONLY: BoundsOfElement
!#if USE_MPI
!USE MOD_Analyze_Vars           ,ONLY: CalcHaloInfo
!#endif /*USE_MPI*/
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL                     :: StartT,EndT
!INTEGER                  :: iElem,ElemToBGM(1:6,1:PP_nElems)
!INTEGER,ALLOCATABLE      :: HaloElemToBGM(:,:)
!REAL,ALLOCATABLE         :: SideOrigin(:,:), SideRadius(:)
!!=================================================================================================================================

!SWRITE(UNIT_StdOut,'(66("-"))')
!SWRITE(UNIT_stdOut,'(A)')' INIT FIBGM...'
!StartT=PICLASTIME()
!!! Read parameter for FastInitBackgroundMesh (FIBGM)
!GEO%FIBGMdeltas(1:3) = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
!GEO%FactorFIBGM(1:3) = GETREALARRAY('Part-FactorFIBGM',3,'1. , 1. , 1.')
!GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

!! build elem basis before halo region build
!ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:nTotalElems) &
!        ,slenXiEtaZetaBasis(1:6,1:nTotalElems) &
!        ,ElemRadiusNGeo(1:nTotalElems)         &
!        ,ElemRadius2NGeo(1:nTotalElems)        )
!CALL BuildElementBasis()

!! sort elem in bgm cells
!SWRITE(UNIT_stdOut,'(A)')' Getting element range in FIBGM ...'
!DO iElem=1,PP_nElems
!  CALL BGMIndexOfElement(iElem,ElemToBGM(1:6,iElem))
!END DO ! iElem = nElems+1,nTotalElems
!SWRITE(UNIT_stdOut,'(A)')' Building FIBGM ...'
!CALL GetFIBGM(ElemToBGM)
!EndT=PICLASTIME()
!IF(PartMPI%MPIROOT)THEN
!  WRITE(UNIT_stdOut,'(A,F12.3,A)',ADVANCE='YES')' Init FIBGM took                  [',EndT-StartT,'s]'
!END IF
!SWRITE(UNIT_StdOut,'(66("-"))')

!CALL DuplicateSlavePeriodicSides()
!! CAUTION: in MarkAllBCSides, a counter is reset for refmapping
!IF(DoRefMapping)THEN
!  CALL MarkAllBCSides()
!END IF
!! get elem and side types
!CALL GetElemAndSideType()

!StartT=PICLASTIME()
!#if USE_MPI
!SWRITE(UNIT_stdOut,'(A)')' INIT HALO REGION...'
!!CALL Initialize()  ! Initialize parallel environment for particle exchange between MPI domains
!printMPINeighborWarnings = GETLOGICAL('printMPINeighborWarnings','.FALSE.')
!printBezierControlPointsWarnings = GETLOGICAL('printBezierControlPointsWarnings','.FALSE.')
!CALL InitHaloMesh()
!! HALO mesh and region build. Unfortunately, the local FIBGM has to be extended to include the HALO elements :(
!! rebuild is a local operation
!IF(.NOT.DoRefMapping)THEN
!  CALL MarkAllBCSides()
!END IF
!#endif /*USE_MPI*/

!IF(nTotalElems.GT.PP_nElems)THEN
!  ALLOCATE(HaloElemToBGM(1:6,PP_nElems+1:nTotalElems))
!  DO iElem=PP_nElems+1,nTotalElems
!    CALL BGMIndexOfElement(iElem,HaloElemToBGM(1:6,iElem))
!  END DO ! iElem = nElems+1,nTotalElems
!  CALL AddHALOCellsToFIBGM(ElemToBGM,HaloElemToBGM)
!  DEALLOCATE(HaloElemToBGM)
!ELSE
!  CALL AddHALOCellsToFIBGM(ElemToBGM)
!END IF

!EndT=PICLASTIME()
!IF(PartMPI%MPIROOT)THEN
!   WRITE(UNIT_stdOut,'(A,F8.3,A)',ADVANCE='YES')' Construction of halo region took [',EndT-StartT,'s]'
!END IF

!! Compute the element bounding boxes before the arrays might be reduced in RefMapping
!ALLOCATE(GEO%BoundsOfElem(1:2,1:3,1:nElems))
!DO iElem = 1, nElems
!  CALL BoundsOfElement(iElem,GEO%BoundsOfElem(1:2,1:3,iElem))
!END DO ! iElem = 1, nElems

!! Reduce the Bezier control point arrays for RefMapping, as only the boundary faces are required for this type of tracking
!IF(DoRefMapping)THEN
!  IF(PartMPI%MPIROOT)THEN
!     WRITE(UNIT_stdOut,'(A)') ' Reshaping arrays to reduced list...'
!  END IF
!  ! remove inner BezierControlPoints3D and SlabNormals, etc.
!  CALL ReShapeBezierSides()
!  ! compute side origin and radius for all sides in PartBCSideList
!  IF(PartMPI%MPIROOT)THEN
!     WRITE(UNIT_stdOut,'(A)') ' GetSideOrigin and Radius..'
!  END IF
!  ! remove inner BezierControlPoints3D and SlabNormals, etc.
!  ALLOCATE( SideOrigin(1:3,1:nTotalBCSides) &
!          , SideRadius(    1:nTotalBCSides) )
!  CALL GetSideOriginAndRadius(nTotalBCSides,SideOrigin,SideRadius)
!END IF

!! Get BCElem mapping, epsOnCell and calculate number of different elements and sides
!IF(PartMPI%MPIROOT)THEN
!   WRITE(UNIT_stdOut,'(A)') ' GetBCElemMap ...'
!END IF
!CALL GetBCElemMap()

!! Identify all elements that are close to boundaries, where the deposition via shape function would cause the shape function sphere
!! to be truncated by the boundary: e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical'
!IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
!  IF(DoSFLocalDepoAtBounds)THEN
!    IF(PartMPI%MPIROOT)THEN
!       WRITE(UNIT_stdOut,'(A)') ' GetShapeFunctionBCElems ...'
!    END IF
!    IF(.NOT.DoRefMapping) CALL abort(&
!    __STAMP__&
!    ,'DoSFLocalDepoAtBounds=T only with DoRefMapping=T because the algorithm requires arrays that are only set for DoRefMapping=T!')
!    CALL GetShapeFunctionBCElems()
!  END IF
!END IF

!! Calculate the number of different element and side types (linear, bi-linear, curved, etc.)
!IF(PartMPI%MPIROOT)THEN
!   WRITE(UNIT_stdOut,'(A)') ' CaclElemAndSideNum ...'
!END IF
!CALL CalcElemAndSideNum()

!! Get basevectors for (bi-)linear sides
!IF(PartMPI%MPIROOT)THEN
!   WRITE(UNIT_stdOut,'(A)') ' LinearSideBaseVectors ...'
!END IF
!CALL GetLinearSideBaseVectors()

!! Set element connectivity
!IF(PartMPI%MPIROOT)THEN
!   WRITE(UNIT_stdOut,'(A)') ' Elem-Connectivity ...'
!END IF
!! check connectivity of particle mesh
!CALL ElemConnectivity()

!IF (FindNeighbourElems) THEN
!  ! build node conectivity of particle mesh
!  IF(PartMPI%MPIROOT)THEN
!     WRITE(UNIT_stdOut,'(A)') ' Node-Neighbourhood ...'
!  END IF
!  CALL NodeNeighbourhood()
!#if USE_MPI
!  IF(PartMPI%MPIROOT)THEN
!     WRITE(UNIT_stdOut,'(A)') ' Node-Communication ...'
!  END IF
!  CALL BuildLocNodeToHaloNodeComm()
!  ! calculate additional volumes and weightings for halo region if deposition type is "cell_volume_mean"
!  IF (TRIM(DepositionType).EQ.'cell_volweight_mean') CALL AddHaloNodeData(CellLocNodes_Volumes)
!#endif /*USE_MPI*/
!END IF

!SDEALLOCATE(XiEtaZetaBasis)
!SDEALLOCATE(slenXiEtaZetaBasis)
!SDEALLOCATE(ElemRadiusNGeo)
!SDEALLOCATE(ElemRadius2NGeo)
!ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:nTotalElems) &
!        ,slenXiEtaZetaBasis(1:6,1:nTotalElems) &
!        ,ElemRadiusNGeo(1:nTotalElems)         &
!        ,ElemRadius2NGeo(1:nTotalElems)        )
!SWRITE(UNIT_stdOut,'(A)')' BUILD ElementBasis ...'
!SDEALLOCATE(ElemRadius2_sf) ! deallocate when using LB (it would be allocated twice because the call is executed twice)
!! second build of elem basis after halo region build
!CALL BuildElementBasis()
!SWRITE(UNIT_stdOut,'(A)')' BUILD ElementBasis DONE!'
!IF(DoRefMapping) THEN
!  ! compute distance between each side associated with  the element and its origin
!  CALL GetElemToSideDistance(nTotalBCSides,SideOrigin,SideRadius)
!  DEALLOCATE( SideOrigin, SideRadius)
!END IF

!#if USE_MPI
!! Output halo element info
!CalcHaloInfo = GETLOGICAL('CalcHaloInfo')
!IF(CalcHaloInfo)THEN
!  CALL SetHaloInfo()
!END IF
!#endif /*USE_MPI*/

!SWRITE(UNIT_stdOut,'(A)')' ... DONE!'
!SWRITE(UNIT_StdOut,'(132("-"))')

!END SUBROUTINE InitFIBGM


!SUBROUTINE GetFIBGM(ElemToBGM)
!!===================================================================================================================================
!! build local FIBGM mesh for process local FIBGM mesh including HALO region
!! mode 1: build local BGM and interconnections with other processes
!! mode 2: rebuild BGM including HALO region
!!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Globals
!USE MOD_Partilce_Periodic_BC ,ONLY: InitPeriodicBC
!USE MOD_Particle_Mesh_Vars   ,ONLY: GEO
!USE MOD_Particle_MPI_Vars    ,ONLY: SafetyFactor,halo_eps_velo,halo_eps,halo_eps2
!#if !(USE_HDG)
!USE MOD_CalcTimeStep         ,ONLY: CalcTimeStep
!#endif /*USE_HDG*/
!USE MOD_Equation_Vars        ,ONLY: c
!USE MOD_Particle_Vars        ,ONLY: manualtimestep
!#if (PP_TimeDiscMethod==201)
!USE MOD_Particle_Vars        ,ONLY: dt_part_ratio
!#endif
!USE MOD_ChangeBasis          ,ONLY: ChangeBasis2D
!#if USE_MPI
!USE MOD_Particle_MPI         ,ONLY: InitHALOMesh
!USE MOD_Particle_Mesh_Vars   ,ONLY: FIBGMCellPadding
!USE MOD_PICDepo_Vars         ,ONLY: DepositionType, r_sf
!USE MOD_Particle_MPI_Vars    ,ONLY: PartMPI
!USE MOD_Particle_Mesh_Vars   ,ONLY: NbrOfCases,casematrix
!#endif /*USE_MPI*/
!#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
!USE MOD_TimeDisc_Vars        ,ONLY: RK_c,nRKStages
!#endif
!USE MOD_ReadInTools          ,ONLY: PrintOption
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!INTEGER,INTENT(IN)    :: mode
!INTEGER,INTENT(IN)     :: ElemToBGM(1:6,1:PP_nElems)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!REAL                  :: localXmin,localXmax,localymin,localymax,localzmin,localzmax
!INTEGER                          :: BGMimin,BGMimax,BGMjmin,BGMjmax,BGMkmin,BGMkmax
!!REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax
!INTEGER                          :: iBGM,jBGM,kBGM,iElem
!INTEGER                          :: BGMCellXmax,BGMCellXmin
!INTEGER                          :: BGMCellYmax,BGMCellYmin
!INTEGER                          :: BGMCellZmax,BGMCellZmin
!INTEGER                          :: ALLOCSTAT
!INTEGER                          :: iProc
!REAL                             :: deltaT
!REAL                             :: globalDiag
!#if USE_MPI
!INTEGER                          :: ii,jj,kk,i,j
!INTEGER                          :: BGMCells,  m, CurrentProc, Cell, Procs
!INTEGER                          :: imin, imax, kmin, kmax, jmin, jmax
!INTEGER                          :: nPaddingCellsX, nPaddingCellsY, nPaddingCellsZ
!INTEGER                          :: nShapePaddingX, nShapePaddingY, nShapePaddingZ
!INTEGER                          :: NbrOfBGMCells(0:PartMPI%nProcs-1)
!INTEGER                          :: Displacement(1:PartMPI%nProcs)
!INTEGER, ALLOCATABLE             :: BGMCellsArray(:),CellProcNum(:,:,:)
!INTEGER, ALLOCATABLE             :: GlobalBGMCellsArray(:), ReducedBGMArray(:)
!INTEGER                          :: ReducedNbrOfBGMCells(0:PartMPI%nProcs-1)
!INTEGER, ALLOCATABLE             :: CellProcList(:,:,:,:)
!INTEGER                          :: tempproclist(0:PartMPI%nProcs-1)
!INTEGER                          :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
!INTEGER                          :: ind, Shift(1:3), iCase
!INTEGER                          :: j_offset
!#endif /*USE_MPI*/
!#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
!INTEGER                          :: iStage
!#endif
!!===================================================================================================================================

!! zeros
!#if USE_MPI
!ii=0
!jj=0
!kk=0
!#endif /*USE_MPI*/


!#if USE_MPI
!! allocate and initialize MPINeighbor
!ALLOCATE(PartMPI%isMPINeighbor(0:PartMPI%nProcs-1))
!PartMPI%isMPINeighbor(:) = .FALSE.
!PartMPI%nMPINeighbors=0
!#endif

!CALL InitPeriodicBC()

!! deallocate stuff // required for dynamic load balance
!#if USE_MPI
!IF (ALLOCATED(GEO%FIBGM)) THEN
!  DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
!    DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
!      DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
!        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
!        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%ShapeProcs)
!        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)
!!           SDEALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs)
!      END DO
!    END DO
!  END DO
!  DEALLOCATE(GEO%FIBGM)
!END IF
!#endif /*USE_MPI*/

!!--- compute number of background cells in each direction
!!BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)+1.00001)
!!BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)+0.99999)
!!BGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)+1.00001)
!!BGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)+0.99999)
!!BGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)+1.00001)
!!BGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)+0.99999)

!! now fail safe, enlarge the BGM grid for safety reasons
!BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
!BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))-1
!BGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
!BGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))-1
!BGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
!BGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))-1

!!--- JN: For MPI communication, information also about the neighboring FIBGM cells is needed
!!--- AS: shouldn't we add up here the nPaddingCells?
!!--- TS: What we need to do is increase the BGM area for shape_function ONLY
!!        Reason: if a particle moves outside the domain, there still needs to be a
!!                BGM with an associated ShapeProc at the particle position
!!        Particle may only move c*dt*Safetyfactor.
!!--- PO: modified for curved and shape-function influence
!!        c*dt*SafetyFactor+r_cutoff
!IF (ManualTimeStep.EQ.0.0) THEN
!#if !(USE_HDG)
!  deltaT=CALCTIMESTEP()
!#else
!   CALL abort(&
!__STAMP__&
!, 'ManualTimeStep.EQ.0.0 -> ManualTimeStep is not defined correctly! Particles-ManualTimeStep = ',RealInfoOpt=ManualTimeStep)
!#endif /*USE_HDG*/
!ELSE
!  deltaT=ManualTimeStep
!END IF
!IF (halo_eps_velo.EQ.0) halo_eps_velo = c
!#if (PP_TimeDiscMethod==4 || PP_TimeDiscMethod==200 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==43)
!IF (halo_eps_velo.EQ.c) THEN
!   CALL abort(&
!__STAMP__&
!, 'halo_eps_velo.EQ.c -> Halo Eps Velocity for MPI not defined')
!END IF
!#endif
!#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
!halo_eps = RK_c(2)
!DO iStage=2,nRKStages-1
!  halo_eps = MAX(halo_eps,RK_c(iStage+1)-RK_c(iStage))
!END DO
!halo_eps = MAX(halo_eps,1.-RK_c(nRKStages))
!CALL PrintOption('max. RKdtFrac','CALCUL.',RealOpt=halo_eps)
!halo_eps = halo_eps*halo_eps_velo*deltaT*SafetyFactor !dt multiplied with maximum RKdtFrac
!#else
!halo_eps = halo_eps_velo*deltaT*SafetyFactor ! for RK too large
!#endif

!#if USE_MPI
!! Check whether halo_eps is smaller than shape function radius
!IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
!  IF(halo_eps.LT.r_sf)THEN
!    SWRITE(UNIT_stdOut,'(A)') ' halo_eps is smaller than shape function radius. Setting halo_eps=r_sf'
!    halo_eps = r_sf
!    CALL PrintOption('max. RKdtFrac','CALCUL.',RealOpt=halo_eps)
!  END IF
!END IF
!#endif /*USE_MPI*/

!! limit halo_eps to diagonal of bounding box
!globalDiag = SQRT( (GEO%xmaxglob-GEO%xminglob)**2 &
!                 + (GEO%ymaxglob-GEO%yminglob)**2 &
!                 + (GEO%zmaxglob-GEO%zminglob)**2 )
!IF(halo_eps.GT.globalDiag)THEN
!  CALL PrintOption('unlimited halo distance','CALCUL.',RealOpt=halo_eps)
!  SWRITE(UNIT_stdOut,'(A38)') ' |   limitation of halo distance  |    '
!  halo_eps=globalDiag
!END IF

!halo_eps2=halo_eps*halo_eps
!CALL PrintOption('halo distance','CALCUL.',RealOpt=halo_eps)


!#if USE_MPI
!! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical', 'shape_function_spherical', 'shape_function_simple'
!IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
!  ! and changed, tooo
!  BGMimax = INT((MIN(GEO%xmax+halo_eps,GEO%xmaxglob)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
!  BGMimin = INT((MAX(GEO%xmin-halo_eps,GEO%xminglob)-GEO%xminglob)/GEO%FIBGMdeltas(1))-1
!  BGMjmax = INT((MIN(GEO%ymax+halo_eps,GEO%ymaxglob)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
!  BGMjmin = INT((MAX(GEO%ymin-halo_eps,GEO%yminglob)-GEO%yminglob)/GEO%FIBGMdeltas(2))-1
!  BGMkmax = INT((MIN(GEO%zmax+halo_eps,GEO%zmaxglob)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
!  BGMkmin = INT((MAX(GEO%zmin-halo_eps,GEO%zminglob)-GEO%zminglob)/GEO%FIBGMdeltas(3))-1
!END IF
!#endif

!GEO%FIBGMimax=BGMimax
!GEO%FIBGMimin=BGMimin
!GEO%FIBGMjmax=BGMjmax
!GEO%FIBGMjmin=BGMjmin
!GEO%FIBGMkmax=BGMkmax
!GEO%FIBGMkmin=BGMkmin


!! allocate space for BGM
!ALLOCATE(GEO%FIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax), STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) THEN
!  WRITE(*,'(A,6(I0,A))')'Problem allocating GEO%FIBGM(',BGMimin,':',BGMimax,',', &
!                                                        BGMjmin,':',BGMjmax,',', &
!                                                        BGMkmin,':',BGMkmax,')'
!#if USE_MPI
!  iProc=PartMPI%MyRank
!#else
!  iProc=0
!#endif /*USE_MPI*/
!  CALL abort(&
!__STAMP__&
!, 'Problem allocating GEO%FIBGM!' )
!END IF

!! null number of element per BGM cell
!DO kBGM = BGMkmin,BGMkmax
!   DO jBGM = BGMjmin,BGMjmax
!     DO iBGM = BGMimin,BGMimax
!         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
!      END DO
!   END DO
!END DO

!!--- compute number of elements in each background cell
!DO iElem=1,PP_nElems
!  ! here fancy stuff, because element could be wide out of element range
!  BGMCellXmin = ElemToBGM(1,iElem)
!  BGMCellXmax = ElemToBGM(2,iElem)
!  BGMCellYmin = ElemToBGM(3,iElem)
!  BGMCellYmax = ElemToBGM(4,iElem)
!  BGMCellZmin = ElemToBGM(5,iElem)
!  BGMCellZmax = ElemToBGM(6,iElem)
!  ! add ecurrent element to number of BGM-elems
!  DO iBGM = BGMCellXmin,BGMCellXmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO kBGM = BGMCellZmin,BGMCellZmax
!         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem

!!--- allocate mapping variable and clean number for mapping (below)
!DO kBGM = BGMkmin,BGMkmax
!  DO jBGM = BGMjmin,BGMjmax
!    DO iBGM = BGMimin,BGMimax
!      IF(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem.EQ.0) CYCLE
!      ALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%FIBGM(iBGM,jBGM,kBGM)%nElem))
!      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
!    END DO ! kBGM
!  END DO ! jBGM
!END DO ! iBGM

!!--- map elements to background cells
!DO iElem=1,PP_nElems
!  ! here fancy stuff, because element could be wide out of element range
!  BGMCellXmin = ElemToBGM(1,iElem)
!  BGMCellXmax = ElemToBGM(2,iElem)
!  BGMCellYmin = ElemToBGM(3,iElem)
!  BGMCellYmax = ElemToBGM(4,iElem)
!  BGMCellZmin = ElemToBGM(5,iElem)
!  BGMCellZmax = ElemToBGM(6,iElem)
!  ! add current Element to BGM-Elem
!  DO kBGM = BGMCellZmin,BGMCellZmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO iBGM = BGMCellXmin,BGMCellXmax
!        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
!        GEO%FIBGM(iBGM,jBGM,kBGM)%Element(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem


!!IF(mode.EQ.2) RETURN
!SWRITE(UNIT_stdOut,'(A)')' Building MPI-FIBGM ...'
!#if USE_MPI
!!--- MPI stuff for background mesh (FastinitBGM)
!BGMCells=0
!ALLOCATE(BGMCellsArray(1:(BGMimax-BGMimin+1)*(BGMjmax-BGMjmin+1)*(BGMkmax-BGMkmin+1)*3))
!!Count BGMCells with Elements inside and save their indices in BGMCellsArray
!DO kBGM=BGMkmin, BGMkmax
!  DO jBGM=BGMjmin, BGMjmax
!    DO iBGM=BGMimin, BGMimax
!      IF (GEO%FIBGM(iBGM,jBGM,kBGM)%nElem .GT. 0) THEN
!        BGMCellsArray(BGMCells*3+1)= iBGM
!        BGMCellsArray(BGMCells*3+2)= jBGM
!        BGMCellsArray(BGMCells*3+3)= kBGM
!        BGMCells=BGMCells+1
!      END IF
!    END DO ! kBGM
!  END DO ! jBGM
!END DO ! iBGM

!!Communicate number of BGMCells
!CALL MPI_ALLGATHER(BGMCells, 1, MPI_INTEGER, NbrOfBGMCells(0:PartMPI%nProcs-1), 1, MPI_INTEGER, PartMPI%COMM, IERROR)
!ALLOCATE(GlobalBGMCellsArray(1:SUM(NbrOfBGMCells)*3))
!Displacement(1)=0
!DO i=2, PartMPI%nProcs
!  Displacement(i) = SUM(NbrOfBGMCells(0:i-2))*3
!END DO
!!Gather indices of every Procs' Cells
!CALL MPI_ALLGATHERV(BGMCellsArray(1:BGMCells*3), BGMCells*3, MPI_INTEGER, GlobalBGMCellsArray, &
!                   & NbrOfBGMCells(0:PartMPI%nProcs-1)*3, Displacement, MPI_INTEGER, PartMPI%COMM, IERROR)

!!--- JN: first: count required array size for ReducedBGMArray
!!--- TS: Define padding stencil (max of halo and shape padding)
!!        Reason: This padding is used to build the ReducedBGM, so any information
!!                outside this region is lost
!IF (GEO%nPeriodicVectors.GT.0) THEN  !Periodic (can't be done below because ReducedBGMArray is sorted by proc)
!  FIBGMCellPadding(1:3)=1
!  IF(.NOT.GEO%directions(1)) FIBGMCellPadding(1) = INT(halo_eps/GEO%FIBGMdeltas(1))+1
!  IF(.NOT.GEO%directions(2)) FIBGMCellPadding(2) = INT(halo_eps/GEO%FIBGMdeltas(2))+1
!  IF(.NOT.GEO%directions(3)) FIBGMCellPadding(3) = INT(halo_eps/GEO%FIBGMdeltas(3))+1
!ELSE
!  FIBGMCellPadding(1:3) = INT(halo_eps/GEO%FIBGMdeltas(1:3))+1
!END IF
!! halo region already included in BGM
!!FIBGMCellPadding(1:3) = 0
!nShapePaddingX = 0
!nShapePaddingY = 0
!nShapePaddingZ = 0
!! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical', 'shape_function_spherical', 'shape_function_simple'
!IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
!  nShapePaddingX = INT(r_sf/GEO%FIBGMdeltas(1)+0.9999999)
!  nShapePaddingY = INT(r_sf/GEO%FIBGMdeltas(2)+0.9999999)
!  nShapePaddingZ = INT(r_sf/GEO%FIBGMdeltas(3)+0.9999999)
!  !IPWRITE(*,*) 'nShapePaddingX',nShapePaddingX
!  !IPWRITE(*,*) 'nShapePaddingY',nShapePaddingY
!  !IPWRITE(*,*) 'nShapePaddingZ',nShapePaddingZ
! ! IF(mode.EQ.2) THEN
! !   IF((nShapePaddingX.EQ.0)    &
! !     .OR.(nShapePaddingY.EQ.0) &
! !     .OR.(nShapePaddingZ.EQ.0))THEN
! !       CALL abort(__STAMP__&
! !         'Error in stencil calculation for FIBGM and shape function')
! !   END IF
! ! END IF
!! 0.999999 in order to prevent stencil to get too big in case of r_sf==c_int*deltas
!!  -> worst case: last 0.000001 gets cut off -> insignificant
!END IF
!nPaddingCellsX = MAX(nShapePaddingX,FIBGMCellPadding(1))
!nPaddingCellsY = MAX(nShapePaddingY,FIBGMCellPadding(2))
!nPaddingCellsZ = MAX(nShapePaddingZ,FIBGMCellPadding(3))

!j=0
!CurrentProc=0
!DO i=1, SUM(NbrOfBGMCells)*3, 3
!  IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
!    CurrentProc=CurrentProc+1
!  END IF
!  IF  (.NOT.(GlobalBGMCellsArray(i) .LT. BGMimin-nPaddingCellsX .OR. GlobalBGMCellsArray(i).GT. BGMimax+nPaddingCellsX &
!      & .OR. GlobalBGMCellsArray(i+1) .LT. BGMjmin-nPaddingCellsY .OR. GlobalBGMCellsArray(i+1) .GT. BGMjmax+nPaddingCellsY &
!      & .OR. GlobalBGMCellsArray(i+2) .LT. BGMkmin-nPaddingCellsZ .OR. GlobalBGMCellsArray(i+2) .GT. BGMkmax+nPaddingCellsZ &
!      & .OR. CurrentProc .EQ. PartMPI%Myrank)) THEN
!    j=j+3
!  END IF
!END DO !i

!! Periodic: ReducedBGMArray needs to include cells on the other side of periodic vectors
!! --- PO: CAUTION: changes throuogh curved
!Vec1(1:3) = 0
!Vec2(1:3) = 0
!Vec3(1:3) = 0
!IF (GEO%nPeriodicVectors.GT.0) THEN
!  ! build case matrix
!  IF (GEO%nPeriodicVectors.EQ.1) THEN
!    DO ind = 1,3
!      Vec1(ind) = INT(GEO%PeriodicVectors(ind,1)/GEO%FIBGMdeltas(ind)+0.1)
!    END DO
!  END IF
!  IF (GEO%nPeriodicVectors.EQ.2) THEN
!    DO ind = 1,3
!      Vec1(ind) = INT(GEO%PeriodicVectors(ind,1)/GEO%FIBGMdeltas(ind)+0.1)
!      Vec2(ind) = INT(GEO%PeriodicVectors(ind,2)/GEO%FIBGMdeltas(ind)+0.1)
!    END DO
!  END IF
!  IF (GEO%nPeriodicVectors.EQ.3) THEN
!    DO ind = 1,3
!      Vec1(ind) = INT(GEO%PeriodicVectors(ind,1)/GEO%FIBGMdeltas(ind)+0.1)
!      Vec2(ind) = INT(GEO%PeriodicVectors(ind,2)/GEO%FIBGMdeltas(ind)+0.1)
!      Vec3(ind) = INT(GEO%PeriodicVectors(ind,3)/GEO%FIBGMdeltas(ind)+0.1)
!    END DO
!  END IF
!  CurrentProc=0
!  DO i=1, SUM(NbrOfBGMCells)*3, 3
!    DO iCase = 1, NbrOfCases
!      IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
!          (casematrix(iCase,2).EQ.0) .AND. &
!          (casematrix(iCase,3).EQ.0)) CYCLE
!      Shift(1:3) = casematrix(iCase,1)*Vec1(1:3) + &
!                   casematrix(iCase,2)*Vec2(1:3) + &
!                   casematrix(iCase,3)*Vec3(1:3)
!      IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
!        CurrentProc=CurrentProc+1
!      END IF
!      IF  (.NOT.(GlobalBGMCellsArray(i)  +Shift(1) .LT. BGMimin-nPaddingCellsX &
!           .OR.  GlobalBGMCellsArray(i)  +Shift(1) .GT. BGMimax+nPaddingCellsX &
!           .OR.  GlobalBGMCellsArray(i+1)+Shift(2) .LT. BGMjmin-nPaddingCellsY &
!           .OR.  GlobalBGMCellsArray(i+1)+Shift(2) .GT. BGMjmax+nPaddingCellsY &
!           .OR.  GlobalBGMCellsArray(i+2)+Shift(3) .LT. BGMkmin-nPaddingCellsZ &
!           .OR.  GlobalBGMCellsArray(i+2)+Shift(3) .GT. BGMkmax+nPaddingCellsZ &
!           .OR. CurrentProc .EQ. PartMPI%MyRank)) THEN
!        j=j+3
!      END IF
!    END DO !iCase
!  END DO !i
!END IF !nPeriodic>0

!ALLOCATE(ReducedBGMArray(1:j))
!!Reduce GlobalBGMCellsArray: erase cells far away from iprocs domain
!!--- JN: ReducedBGMArray contains data only from other MPI procs!

!IF (GEO%nPeriodicVectors.GT.0) THEN  !Periodic (can't be done below because ReducedBGMArray is sorted by proc)
!  j=1
!  CurrentProc=0
!  ReducedBGMArray=0
!  ReducedNbrOfBGMCells=0
!  DO i=1, SUM(NbrOfBGMCells)*3, 3
!    DO iCase = 1, NbrOfCases         ! This time INCLUDING non-moved
!      Shift(1:3) = casematrix(iCase,1)*Vec1(1:3) + &
!                   casematrix(iCase,2)*Vec2(1:3) + &
!                   casematrix(iCase,3)*Vec3(1:3)
!      IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
!        CurrentProc=CurrentProc+1
!      END IF
!      IF  (.NOT.(GlobalBGMCellsArray(i)   +Shift(1) .LT. BGMimin-nPaddingCellsX &
!           .OR.  GlobalBGMCellsArray(i)   +Shift(1) .GT. BGMimax+nPaddingCellsX &
!           .OR.  GlobalBGMCellsArray(i+1) +Shift(2) .LT. BGMjmin-nPaddingCellsY &
!           .OR.  GlobalBGMCellsArray(i+1) +Shift(2) .GT. BGMjmax+nPaddingCellsY &
!           .OR.  GlobalBGMCellsArray(i+2) +Shift(3) .LT. BGMkmin-nPaddingCellsZ &
!           .OR.  GlobalBGMCellsArray(i+2) +Shift(3) .GT. BGMkmax+nPaddingCellsZ &
!           .OR.  CurrentProc .EQ. PartMPI%MyRank)) THEN
!        ReducedBGMArray(j)=GlobalBGMCellsArray(i)     +Shift(1)
!        ReducedBGMArray(j+1)=GlobalBGMCellsArray(i+1) +Shift(2)
!        ReducedBGMArray(j+2)=GlobalBGMCellsArray(i+2) +Shift(3)
!        j=j+3
!        ReducedNbrOfBGMCells(CurrentProc)=ReducedNbrOfBGMCells(CurrentProc)+1
!      END IF
!    END DO ! iCase
!  END DO !i
!ELSE ! non periodic case
!  j=1
!  CurrentProc=0
!  ReducedBGMArray=0
!  ReducedNbrOfBGMCells=0
!  DO i=1, SUM(NbrOfBGMCells)*3, 3
!    IF  (i .GT. SUM(NbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PartMPI%nProcs-1) THEN
!      CurrentProc=CurrentProc+1
!    END IF
!    IF  (.NOT.(GlobalBGMCellsArray(i)   .LT. BGMimin-nPaddingCellsX .OR. GlobalBGMCellsArray(i).GT.    BGMimax+nPaddingCellsX &
!        & .OR. GlobalBGMCellsArray(i+1) .LT. BGMjmin-nPaddingCellsY .OR. GlobalBGMCellsArray(i+1) .GT. BGMjmax+nPaddingCellsY &
!        & .OR. GlobalBGMCellsArray(i+2) .LT. BGMkmin-nPaddingCellsZ .OR. GlobalBGMCellsArray(i+2) .GT. BGMkmax+nPaddingCellsZ &
!         & .OR. CurrentProc .EQ. PartMPI%MyRank)) THEN
!      ReducedBGMArray(j  )=GlobalBGMCellsArray(i  )
!      ReducedBGMArray(j+1)=GlobalBGMCellsArray(i+1)
!      ReducedBGMArray(j+2)=GlobalBGMCellsArray(i+2)
!      j=j+3
!      ReducedNbrOfBGMCells(CurrentProc)=ReducedNbrOfBGMCells(CurrentProc)+1
!    END IF
!  END DO !i
!END IF !periodic


!!--- JN: Determine required size of CellProcList array (hope this works, everytime I try to again understand this
!!        shape function parallelization stuff, I get confused...)
!!--- JN: But therefore we first have to refill BGMCellsArray to not only contain
!!        cells with PIC%FastInitBGM%nElem.GT.0 but also those adjacent to them!
!!--- TS: Actually, not the adjacent cell needs to be considered but a shape_proc stencil
!!        Usually, the shape function radius is chosen to be the size of one BGM, but this
!!        is not necessarily always true. Hence new shape_proc padding:

!BGMCells=0
!DO iBGM=BGMimin, BGMimax  !Count BGMCells with Elements inside or adjacent and save their indices in BGMCellsArray
!  DO jBGM=BGMjmin, BGMjmax
!    DO kBGM=BGMkmin, BGMkmax
!      iMin=MAX(iBGM-nShapePaddingX,BGMimin); iMax=MIN(iBGM+nShapePaddingX,BGMimax)
!      jMin=MAX(jBGM-nShapePaddingY,BGMjmin); jMax=MIN(jBGM+nShapePaddingY,BGMjmax)
!      kMin=MAX(kBGM-nShapePaddingZ,BGMkmin); kMax=MIN(kBGM+nShapePaddingZ,BGMkmax)
!      IF (SUM(GEO%FIBGM(iMin:iMax,jMin:jMax,kMin:kMax)%nElem) .GT. 0) THEN
!        ! debug here changed i,j,k to ibgm,jbgm,kbgm
!        BGMCellsArray(BGMCells*3+1)= iBGM
!        BGMCellsArray(BGMCells*3+2)= jBGM
!        BGMCellsArray(BGMCells*3+3)= kBGM
!        BGMCells=BGMCells+1
!      END IF
!    END DO !iBGM
!  END DO !jBGM
!END DO !kBGM

!! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved
!! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

!! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
!ALLOCATE(CellProcNum(BGMimin-nShapePaddingX:BGMimax+nShapePaddingX, &
!                     BGMjmin-nShapePaddingY:BGMjmax+nShapePaddingY, &
!                     BGMkmin-nShapePaddingZ:BGMkmax+nShapePaddingZ))
!CellProcNum = 0
!Procs = 0 ! = maximum number of procs in one BGM cell
!DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
!  IF((ReducedBGMArray(j).GE.BGMimin-nShapePaddingX).AND.(ReducedBGMArray(j).LE.BGMimax+nShapePaddingX))THEN
!    IF((ReducedBGMArray(j+1).GE.BGMjmin-nShapePaddingY).AND.(ReducedBGMArray(j+1).LE.BGMjmax+nShapePaddingY))THEN
!      IF((ReducedBGMArray(j+2).GE.BGMkmin-nShapePaddingZ).AND.(ReducedBGMArray(j+2).LE.BGMkmax+nShapePaddingZ))THEN !inside
!        CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
!             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
!        Procs = MAX(Procs, CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)))
!      END IF
!    END IF
!  END IF
!END DO
!! allocate the temporary array
!ALLOCATE(CellProcList(BGMimin-nShapePaddingX:BGMimax+nShapePaddingX, &
!                      BGMjmin-nShapePaddingY:BGMjmax+nShapePaddingY, &
!                      BGMkmin-nShapePaddingZ:BGMkmax+nShapePaddingZ, &
!                      1:Procs))
!CellProcList = -1

!! fill array with proc numbers

!CellProcNum = 0
!j_offset = 0
!DO CurrentProc = 0,PartMPI%nProcs-1
!  DO j = 1+j_offset, ReducedNbrOfBGMCells(CurrentProc)*3-2+j_offset,3
!    IF((ReducedBGMArray(j).GE.BGMimin-nShapePaddingX).AND.(ReducedBGMArray(j).LE.BGMimax+nShapePaddingX))THEN
!      IF((ReducedBGMArray(j+1).GE.BGMjmin-nShapePaddingY).AND.(ReducedBGMArray(j+1).LE.BGMjmax+nShapePaddingY))THEN
!        IF((ReducedBGMArray(j+2).GE.BGMkmin-nShapePaddingZ).AND.(ReducedBGMArray(j+2).LE.BGMkmax+nShapePaddingZ))THEN
!          CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
!             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
!          CellProcList(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2), &
!             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2))) = CurrentProc
!        END IF
!      END IF
!    END IF
!  END DO
!  j_offset = j_offset + ReducedNbrOfBGMCells(CurrentProc)*3
!END DO
!! fill real array
!DO Cell=0, BGMCells-1
!  TempProcList=0
!  DO iBGM = BGMCellsArray(Cell*3+1)-nShapePaddingX, BGMCellsArray(Cell*3+1)+nShapePaddingX
!    DO jBGM = BGMCellsArray(Cell*3+2)-nShapePaddingY, BGMCellsArray(Cell*3+2)+nShapePaddingY
!      DO kBGM = BGMCellsArray(Cell*3+3)-nShapePaddingZ, BGMCellsArray(Cell*3+3)+nShapePaddingZ
!        DO m = 1,CellProcNum(iBGM,jBGM,kBGM)
!          TempProcList(CellProcList(iBGM,jBGM,kBGM,m))=1       ! every proc that is within the stencil gets a 1
!        END DO ! m
!        kk = kBGM
!      END DO !kBGM
!      jj = jBGM
!    END DO !jBGM
!    ii = iBGM
!  END DO !iBGM
!  Procs=SUM(TempProcList)
!  IF (Procs.NE.0) THEN
!    ALLOCATE(GEO%FIBGM(ii-nShapePaddingX,jj-nShapePaddingY,kk-nShapePaddingZ)%ShapeProcs(1:Procs+1))
!    GEO%FIBGM(ii-nShapePaddingX,jj-nShapePaddingY,kk-nShapePaddingZ)%ShapeProcs(1) = Procs
!    j=2
!    DO m=0,PartMPI%nProcs-1
!      IF (TempProcList(m) .EQ. 1) THEN
!        IF(.NOT.PartMPI%isMPINeighbor(m))THEN
!          !IF(mode.EQ.2)THEN
!          !  IPWRITE(UNIT_stdOut,*) ' Warning, something wrong with halo region'
!          !  CALL abort(__STAMP__&
!          !      , ' Something wrong with Halo region' )
!          !END IF
!          PartMPI%isMPINeighbor(m) = .true.
!          PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
!        END IF
!        GEO%FIBGM(ii-nShapePaddingX,jj-nShapePaddingY,kk-nShapePaddingZ)%ShapeProcs(j)=m
!        j=j+1
!      END IF
!    END DO !m
!  END IF
!END DO !Cell

!   !Compare own BGMCells and their Neighbors with ReducedBGMArray and save other Processes in BGM-Cells
!   !--- JN: ReducedBGMArray contains data only from other MPI procs!
!   !--- JN: BGMCellsArray contains in index triplets (i,k,l) all BGM cells containing elements from the local MPI proc
!   !        plus the index triplets of BGM cells adjacent to cells containing elements from the local MPI proc

!!   !--- JN: First identify only procs that share the exact same BGM cell as I (SharedProcs)
!!   Procs = 0
!!   CellProcList=-1
!!   DO Cell=0, BGMCells-1
!!     TempProcList=0
!!     i = BGMCellsArray(Cell*3+1)
!!     k = BGMCellsArray(Cell*3+2)
!!     l = BGMCellsArray(Cell*3+3)
!!     IF (GEO%FIBGM(i,k,l)%nElem.EQ.0) CYCLE
!!     CurrentProc=0
!!     m=2
!!     DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
!!       !--- JN: Slide CurrentProc to the MPI Proc that the currently checked BGMCell belongs to
!!       DO WHILE (j .GT. SUM(ReducedNbrOfBGMCells(0: CurrentProc))*3 .AND. CurrentProc .LT. PMPIVAR%nProcs-1)
!!         CurrentProc=CurrentProc+1
!!       END DO
!!       IF (i .EQ. ReducedBGMArray(j) .AND. k .EQ. ReducedBGMArray(j+1) .AND. l .EQ. ReducedBGMArray(j+2)) THEN
!!         IF (m .GT. MaxShapeProcs) THEN
!!           CALL abort(__STAMP__&
!!                                'ERROR in Boundary_PIC.f90: Cellproclist can contain only MaxShapeProcs=',MaxShapeProcs,999.)
!!         END IF
!!         CellProcList(i,k,l,m)=CurrentProc
!!         m=m+1
!!         TempProcList(CurrentProc)=1
!!       END IF
!!     END DO !j
!!     Procs=SUM(TempProcList)
!!     ALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs(1:Procs+1))
!!     GEO%FIBGM(i,k,l)%SharedProcs(1) = Procs
!!     j=2
!!     DO m=0,PMPIVAR%nProcs-1
!!       IF (TempProcList(m) .EQ. 1) THEN
!!         GEO%FIBGM(i,k,l)%SharedProcs(j)=m
!!         j=j+1
!!       END IF
!!     END DO !m
!!   END DO !Cell


!! ----------------------------------------------------------------!
!!--- AS: Do it again for Paddingcells
!DEALLOCATE(CellProcList)
!DEALLOCATE(CellProcNum)
!!--- JN: Determine required size of CellProcList array (hope this works, everytime I try to again understand this
!!        shape function parallelization stuff, I get confused...)
!!--- JN: But therefore we first have to refill BGMCellsArray to not only contain
!!        cells with PIC%FastInitBGM%nElem.GT.0 but also those adjacent and the paddingcells to them!
!BGMCells=0
!DO iBGM=BGMimin, BGMimax  !Count BGMCells with Elements inside or adjacent and save their indices in BGMCellsArray
!  DO jBGM=BGMjmin, BGMjmax
!    DO kBGM=BGMkmin, BGMkmax
!      iMin=MAX(iBGM-nPaddingCellsX,BGMimin); iMax=MIN(iBGM+nPaddingCellsX,BGMimax)
!      jMin=MAX(jBGM-nPaddingCellsY,BGMjmin); jMax=MIN(jBGM+nPaddingCellsY,BGMjmax)
!      kMin=MAX(kBGM-nPaddingCellsZ,BGMkmin); kMax=MIN(kBGM+nPaddingCellsZ,BGMkmax)
!      IF (SUM(GEO%FIBGM(iMin:iMax,jMin:jMax,kMin:kMax)%nElem) .GT. 0) THEN
!        BGMCellsArray(BGMCells*3+1)= iBGM
!        BGMCellsArray(BGMCells*3+2)= jBGM
!        BGMCellsArray(BGMCells*3+3)= kBGM
!        BGMCells=BGMCells+1
!      END IF
!    END DO !iBGM
!  END DO !jBGM
!END DO !kBGM

!! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved
!! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

!! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
!ALLOCATE(CellProcNum(BGMimin-nPaddingCellsX:BGMimax+nPaddingCellsX, &
!                     BGMjmin-nPaddingCellsY:BGMjmax+nPaddingCellsY, &
!                     BGMkmin-nPaddingCellsZ:BGMkmax+nPaddingCellsZ))
!CellProcNum = 0
!Procs = 0
!DO j=1, SUM(ReducedNbrOfBGMCells)*3-2, 3
!   IF((ReducedBGMArray(j).GE.BGMimin-nPaddingCellsX).AND.(ReducedBGMArray(j).LE.BGMimax+nPaddingCellsX))THEN
!     IF((ReducedBGMArray(j+1).GE.BGMjmin-nPaddingCellsY).AND.(ReducedBGMArray(j+1).LE.BGMjmax+nPaddingCellsY))THEN
!       IF((ReducedBGMArray(j+2).GE.BGMkmin-nPaddingCellsZ).AND.(ReducedBGMArray(j+2).LE.BGMkmax+nPaddingCellsZ))THEN
!        CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
!             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
!        Procs = MAX(Procs, CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)))
!       END IF
!     END IF
!   END IF
!END DO
!! allocate the temporary array
!ALLOCATE(CellProcList(BGMimin-nPaddingCellsX:BGMimax+nPaddingCellsX, &
!                      BGMjmin-nPaddingCellsY:BGMjmax+nPaddingCellsY, &
!                      BGMkmin-nPaddingCellsZ:BGMkmax+nPaddingCellsZ, &
!                      1:Procs))
!CellProcList = -1

!! fill array with proc numbers

!CellProcNum = 0
!j_offset = 0
!DO CurrentProc = 0,PartMPI%nProcs-1
!  DO j = 1+j_offset, j_offset+ReducedNbrOfBGMCells(CurrentProc)*3-2,3
!    CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) = &
!             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2)) + 1
!    CellProcList(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2), &
!             CellProcNum(ReducedBGMArray(j),ReducedBGMArray(j+1),ReducedBGMArray(j+2))) = CurrentProc
!  END DO
!  j_offset = j_offset + ReducedNbrOfBGMCells(CurrentProc)*3
!END DO

!! fill real array
!DO Cell=0, BGMCells-1
!  TempProcList=0
!  DO iBGM = BGMCellsArray(Cell*3+1)-nPaddingCellsX, BGMCellsArray(Cell*3+1)+nPaddingCellsX
!    DO jBGM = BGMCellsArray(Cell*3+2)-nPaddingCellsY, BGMCellsArray(Cell*3+2)+nPaddingCellsY
!      DO kBGM = BGMCellsArray(Cell*3+3)-nPaddingCellsZ, BGMCellsArray(Cell*3+3)+nPaddingCellsZ
!        DO m = 1,CellProcNum(iBGM,jBGM,kBGM)
!          TempProcList(CellProcList(iBGM,jBGM,kBGM,m))=1       ! every proc that is within the stencil gets a 1
!        END DO ! m
!        kk = kBGM
!      END DO !l
!      jj = jBGM
!    END DO !k
!    ii = iBGM
!  END DO !i
!  Procs=SUM(TempProcList)
!  IF (Procs.NE.0) THEN
!    ALLOCATE(GEO%FIBGM(ii-nPaddingCellsX,jj-nPaddingCellsY,kk-nPaddingCellsZ)%PaddingProcs(1:Procs+1))
!    GEO%FIBGM(ii-nPaddingCellsX,jj-nPaddingCellsY,kk-nPaddingCellsZ)%PaddingProcs(1) = Procs
!    j=2
!    DO m=0,PartMPI%nProcs-1
!      IF (TempProcList(m) .EQ. 1) THEN
!        GEO%FIBGM(ii-nPaddingCellsX,jj-nPaddingCellsY,kk-nPaddingCellsZ)%PaddingProcs(j)=m
!        j=j+1
!      END IF
!    END DO !m
!  END IF
!END DO !Cell
!DEALLOCATE(ReducedBGMArray, BGMCellsArray, CellProcList, GlobalBGMCellsArray, CellProcNum)
!#endif /*USE_MPI*/

!END SUBROUTINE GetFIBGM



!SUBROUTINE AddHALOCellsToFIBGM(ElemToBGM,HaloElemToBGM)
!!===================================================================================================================================
!! remap all elements including halo-elements into FIBGM
!!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Globals                                                        ! ,            ONLY : UNIT_StdOut
!USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
!USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,nTotalElems
!USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
!USE MOD_Particle_Tracking_Vars ,ONLY: Distance,ListDistance
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!INTEGER,INTENT(IN)    :: mode
!INTEGER,INTENT(IN)               :: ElemToBGM(1:6,1:PP_nElems)
!INTEGER,INTENT(IN),OPTIONAL      :: HaloElemToBGM(1:6,PP_nElems+1:nTotalElems)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                          :: BGMimin,BGMimax,BGMjmin,BGMjmax,BGMkmin,BGMkmax,Allocstat
!REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax
!INTEGER                          :: iBGM,jBGM,kBGM,iElem
!INTEGER                          :: BGMCellXmax,BGMCellXmin
!INTEGER                          :: BGMCellYmax,BGMCellYmin
!INTEGER                          :: BGMCellZmax,BGMCellZmin
!LOGICAL, ALLOCATABLE             :: ElementFound(:)
!INTEGER                          :: maxnBGMElems
!!===================================================================================================================================


!! current min,max
!BGMimax=GEO%FIBGMimax
!BGMimin=GEO%FIBGMimin
!BGMjmax=GEO%FIBGMjmax
!BGMjmin=GEO%FIBGMjmin
!BGMkmax=GEO%FIBGMkmax
!BGMkmin=GEO%FIBGMkmin

!GEO%TFIBGMimax =GEO%FIBGMimax
!GEO%TFIBGMimin =GEO%FIBGMimin
!GEO%TFIBGMjmax =GEO%FIBGMjmax
!GEO%TFIBGMjmin =GEO%FIBGMjmin
!GEO%TFIBGMkmax =GEO%FIBGMkmax
!GEO%TFIBGMkmin =GEO%FIBGMkmin

!BGMCellXmax = BGMimax
!BGMCellXmin = BGMimin
!BGMCellYmax = BGMjmax
!BGMCellYmin = BGMjmin
!BGMCellZmax = BGMkmax
!BGMCellZmin = BGMkmin


!DO iElem=1,nTotalElems
!  IF(iElem.LE.PP_nElems)THEN
!    BGMCellXmin = ElemToBGM(1,iElem)
!    BGMCellXmax = ElemToBGM(2,iElem)
!    BGMCellYmin = ElemToBGM(3,iElem)
!    BGMCellYmax = ElemToBGM(4,iElem)
!    BGMCellZmin = ElemToBGM(5,iElem)
!    BGMCellZmax = ElemToBGM(6,iElem)
!  ELSE
!    IF(.NOT.GEO%directions(1)) BGMCellXmin = HaloElemToBGM(1,iElem)
!    IF(.NOT.GEO%directions(1)) BGMCellXmax = HaloElemToBGM(2,iElem)
!    IF(.NOT.GEO%directions(2)) BGMCellYmin = HaloElemToBGM(3,iElem)
!    IF(.NOT.GEO%directions(2)) BGMCellYmax = HaloElemToBGM(4,iElem)
!    IF(.NOT.GEO%directions(3)) BGMCellZmin = HaloElemToBGM(5,iElem)
!    IF(.NOT.GEO%directions(3)) BGMCellZmax = HaloElemToBGM(6,iElem)
!  END IF

!  BGMimin=MIN(BGMimin,BGMCellXmin)
!  BGMimax=MAX(BGMimax,BGMCellXmax)
!  BGMjmin=MIN(BGMjmin,BGMCellYmin)
!  BGMjmax=MAX(BGMjmax,BGMCellYmax)
!  BGMkmin=MIN(BGMkmin,BGMCellZmin)
!  BGMkmax=MAX(BGMkmax,BGMCellZmax)

!END DO ! iElem = nElems+1,nTotalElems

!GEO%TFIBGMimax =BGMimax
!GEO%TFIBGMimin =BGMimin
!GEO%TFIBGMjmax =BGMjmax
!GEO%TFIBGMjmin =BGMjmin
!GEO%TFIBGMkmax =BGMkmax
!GEO%TFIBGMkmin =BGMkmin

!ALLOCATE(GEO%TFIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax), STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) THEN
!    CALL abort(&
!__STAMP__&
!,' ERROR in AddElemsToTFIBGM: Cannot allocate GEO%TFIBGM!')
!END IF

!ALLOCATE( ElementFound(1:nTotalElems) )
!ElementFound = .FALSE.

!! null number of elements per BGM-Cell
!DO kBGM = BGMkmin,BGMkmax
!  DO jBGM = BGMjmin,BGMjmax
!    DO iBGM = BGMimin,BGMimax
!       GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = 0
!    END DO ! kBGM
!  END DO ! jBGM
!END DO ! iBGM


!!--- compute number of elements in each background cell
!DO iElem=1,PP_nElems
!  !--- find minimum and maximum BGM cell for current element
!  ! here fancy stuff, because element could be wide out of element range
!  BGMCellXmin = MIN(MAX(ElemToBGM(1,iElem),BGMimin),BGMimax)
!  BGMCellXmax = MAX(MIN(ElemToBGM(2,iElem),BGMimax),BGMimin)
!  BGMCellYmin = MIN(MAX(ElemToBGM(3,iElem),BGMjmin),BGMjmax)
!  BGMCellYmax = MAX(MIN(ElemToBGM(4,iElem),BGMjmax),BGMjmin)
!  BGMCellZmin = MIN(MAX(ElemToBGM(5,iElem),BGMkmin),BGMkmax)
!  BGMCellZmax = MAX(MIN(ElemToBGM(6,iElem),BGMkmax),BGMkmin)
!  ! add ecurrent element to number of BGM-elems
!  DO kBGM = BGMCellZmin,BGMCellZmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO iBGM = BGMCellXmin,BGMCellXmax
!         GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1
!         ElementFound(iElem) = .TRUE.
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem

!DO iElem=PP_nElems+1,nTotalElems
!  !--- find minimum and maximum BGM cell for current element
!  ! here fancy stuff, because element could be wide out of element range
!  BGMCellXmin = MIN(MAX(HaloElemToBGM(1,iElem),BGMimin),BGMimax)
!  BGMCellXmax = MAX(MIN(HaloElemToBGM(2,iElem),BGMimax),BGMimin)
!  BGMCellYmin = MIN(MAX(HaloElemToBGM(3,iElem),BGMjmin),BGMjmax)
!  BGMCellYmax = MAX(MIN(HaloElemToBGM(4,iElem),BGMjmax),BGMjmin)
!  BGMCellZmin = MIN(MAX(HaloElemToBGM(5,iElem),BGMkmin),BGMkmax)
!  BGMCellZmax = MAX(MIN(HaloElemToBGM(6,iElem),BGMkmax),BGMkmin)
!  ! add ecurrent element to number of BGM-elems
!  DO kBGM = BGMCellZmin,BGMCellZmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO iBGM = BGMCellXmin,BGMCellXmax
!        GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1
!        ElementFound(iElem) = .TRUE.
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem


!!--- allocate mapping variable and clean number for mapping (below)
!DO kBGM = BGMkmin,BGMkmax
!  DO jBGM = BGMjmin,BGMjmax
!    DO iBGM = BGMimin,BGMimax
!      IF(GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem.EQ.0) CYCLE
!      ALLOCATE(GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem))
!      GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = 0
!    END DO ! kBGM
!  END DO ! jBGM
!END DO ! iBGM

!!--- map elements to background cells
!DO iElem=1,PP_nElems
!  !--- find minimum and maximum BGM cell for current element
!  ! here fancy stuff, because element could be wide out of element range
!  BGMCellXmin = MIN(MAX(ElemToBGM(1,iElem),BGMimin),BGMimax)
!  BGMCellXmax = MAX(MIN(ElemToBGM(2,iElem),BGMimax),BGMimin)
!  BGMCellYmin = MIN(MAX(ElemToBGM(3,iElem),BGMjmin),BGMjmax)
!  BGMCellYmax = MAX(MIN(ElemToBGM(4,iElem),BGMjmax),BGMjmin)
!  BGMCellZmin = MIN(MAX(ElemToBGM(5,iElem),BGMkmin),BGMkmax)
!  BGMCellZmax = MAX(MIN(ElemToBGM(6,iElem),BGMkmax),BGMkmin)

!  ! add current Element to BGM-Elem
!  DO iBGM = BGMCellXmin,BGMCellXmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO kBGM = BGMCellZmin,BGMCellZmax
!        GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1
!        GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem
!DO iElem=PP_nElems+1,nTotalElems
!  !--- find minimum and maximum BGM cell for current element
!  ! here fancy stuff, because element could be wide out of element range
!  BGMCellXmin = MIN(MAX(HaloElemToBGM(1,iElem),BGMimin),BGMimax)
!  BGMCellXmax = MAX(MIN(HaloElemToBGM(2,iElem),BGMimax),BGMimin)
!  BGMCellYmin = MIN(MAX(HaloElemToBGM(3,iElem),BGMjmin),BGMjmax)
!  BGMCellYmax = MAX(MIN(HaloElemToBGM(4,iElem),BGMjmax),BGMjmin)
!  BGMCellZmin = MIN(MAX(HaloElemToBGM(5,iElem),BGMkmin),BGMkmax)
!  BGMCellZmax = MAX(MIN(HaloElemToBGM(6,iElem),BGMkmax),BGMkmin)

!  ! add current Element to BGM-Elem
!  DO iBGM = BGMCellXmin,BGMCellXmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO kBGM = BGMCellZmin,BGMCellZmax
!        GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1
!        GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem


!DO iElem=1,PP_nElems
!  IF(.NOT.ElementFound(iElem))THEN
!    !--- find minimum and maximum BGM cell for current element
!    ! here fancy stuff, because element could be wide out of element range
!    BGMCellXmin = ElemToBGM(1,iElem)
!    BGMCellXmax = ElemToBGM(2,iElem)
!    BGMCellYmin = ElemToBGM(3,iElem)
!    BGMCellYmax = ElemToBGM(4,iElem)
!    BGMCellZmin = ElemToBGM(5,iElem)
!    BGMCellZmax = ElemToBGM(6,iElem)

!    IPWRITE(UNIT_stdOut,*) ' TFIBGM , iElem'
!    IPWRITE(UNIT_stdOut,*) 'xmin',GEO%xmin
!    IPWRITE(UNIT_stdOut,*) 'xmax',GEO%xmax
!    IPWRITE(UNIT_stdOut,*) 'ymin',GEO%ymin
!    IPWRITE(UNIT_stdOut,*) 'ymax',GEO%ymax
!    IPWRITE(UNIT_stdOut,*) 'zmin',GEO%zmin
!    IPWRITE(UNIT_stdOut,*) 'zmax',GEO%zmax
!    IPWRITE(UNIT_stdOut,*) ' BGM , iBGM'
!    IPWRITE(UNIT_stdOut,*) 'xmin', BGMimin,BGMCellXmin
!    IPWRITE(UNIT_stdOut,*) 'xmax', BGMimax,BGMCellXmax
!    IPWRITE(UNIT_stdOut,*) 'ymin', BGMjmin,BGMCellYmin
!    IPWRITE(UNIT_stdOut,*) 'ymax', BGMjmax,BGMCellYmax
!    IPWRITE(UNIT_stdOut,*) 'zmin', BGMkmin,BGMCellYmin
!    IPWRITE(UNIT_stdOut,*) 'zmax', BGMkmax,BGMCellYmax
!    CALL abort(&
!__STAMP__&
!,' Element not located in FIBGM! iElem, myRank',iElem,REAL(PartMPI%MyRank))
!  END IF
!END DO ! iElem

!DO iElem=PP_nElems+1,nTotalElems
!  IF(.NOT.ElementFound(iElem))THEN
!    !--- find minimum and maximum BGM cell for current element
!    ! here fancy stuff, because element could be wide out of element range
!    BGMCellXmin = HaloElemToBGM(1,iElem)
!    BGMCellXmax = HaloElemToBGM(2,iElem)
!    BGMCellYmin = HaloElemToBGM(3,iElem)
!    BGMCellYmax = HaloElemToBGM(4,iElem)
!    BGMCellZmin = HaloElemToBGM(5,iElem)
!    BGMCellZmax = HaloElemToBGM(6,iElem)

!    IPWRITE(UNIT_stdOut,*) ' TFIBGM , iElem'
!    IPWRITE(UNIT_stdOut,*) 'xmin',GEO%xmin,xmin
!    IPWRITE(UNIT_stdOut,*) 'xmax',GEO%xmax,xmax
!    IPWRITE(UNIT_stdOut,*) 'ymin',GEO%ymin,ymin
!    IPWRITE(UNIT_stdOut,*) 'ymax',GEO%ymax,ymax
!    IPWRITE(UNIT_stdOut,*) 'zmin',GEO%zmin,zmin
!    IPWRITE(UNIT_stdOut,*) 'zmax',GEO%zmax,zmax
!    IPWRITE(UNIT_stdOut,*) ' BGM , iBGM'
!    IPWRITE(UNIT_stdOut,*) 'xmin', BGMimin,BGMCellXmin
!    IPWRITE(UNIT_stdOut,*) 'xmax', BGMimax,BGMCellXmax
!    IPWRITE(UNIT_stdOut,*) 'ymin', BGMjmin,BGMCellYmin
!    IPWRITE(UNIT_stdOut,*) 'ymax', BGMjmax,BGMCellYmax
!    IPWRITE(UNIT_stdOut,*) 'zmin', BGMkmin,BGMCellYmin
!    IPWRITE(UNIT_stdOut,*) 'zmax', BGMkmax,BGMCellYmax
!    CALL abort(&
!__STAMP__&
!,' Element not located in FIBGM! iElem, myRank',iElem,REAL(PartMPI%MyRank))
!  END IF
!END DO ! iElem


!DEALLOCATE(Elementfound)

!! and get max number of bgm-elems
!maxnBGMElems=0
!DO kBGM = GEO%TFIBGMkmin,GEO%TFIBGMkmax
!  DO jBGM = GEO%TFIBGMjmin,GEO%TFIBGMjmax
!    DO iBGM = GEO%TFIBGMimin,GEO%TFIBGMimax
!      !maxnBGMElems=MAX(maxnBGMElems,GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem)
!      maxnBGMElems=MAX(maxnBGMElems,GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem)
!    END DO ! kBGM
!  END DO ! jBGM
!END DO ! iBGM
!ALLOCATE(Distance    (1:maxnBGMElems) &
!        ,ListDistance(1:maxnBGMElems) )


!END SUBROUTINE AddHALOCellsToFIBGM





!SUBROUTINE ReShapeBezierSides()
!!===================================================================================================================================
!! Init of Particle mesh
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalBCSides,PartBCSideList,nTotalSides,nPartPeriodicSides
!USE MOD_Mesh_Vars              ,ONLY: nSides,nBCSides,NGeo,BC
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,SideType,SideDistance,SideNormVec
!USE MOD_Particle_Surfaces_Vars ,ONLY: SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER           :: ALLOCSTAT
!INTEGER           :: iSide,nOldBCSides,newBCSideID,BCInc,nPeriodicSidesTmp
!REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummySideSlabNormals                  ! normal vectors of bounding slab box
!REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummySideSlabIntervals               ! intervalls beta1, beta2, beta3
!LOGICAL,ALLOCATABLE,DIMENSION(:)   :: DummyBoundingBoxIsEmpty
!REAL,ALLOCATABLE                   :: DummyBezierControlPoints3D(:,:,:,:)
!INTEGER,ALLOCATABLE                :: DummySideType(:)
!REAL,ALLOCATABLE                   :: DummySideDistance(:)
!REAL,ALLOCATABLE                   :: DummySideNormVec(:,:)
!!===================================================================================================================================

!nPeriodicSidesTmp=0
!DO iSide=nBCSides+1,nSides+nPartPeriodicSides
!  IF(BC(iSide).NE.0)THEN
!    ! different list, contains ALL periodic sides (inner and duplicated)
!    nPeriodicSidesTmp=nPeriodicSidesTmp+1
!  END IF
!END DO

!! now, shrink partbcsidelist
!nOldBCSides  =nTotalBCSides
!nTotalBCSides=nTotalSides-nPartPeriodicSides-nSides+nBCSides+nPeriodicSidesTmp

!IF(nTotalBCSides.EQ.0) RETURN

!! allocate & fill dummy
!! BezierControlPoints3D
!ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides))
!IF (.NOT.ALLOCATED(DummyBezierControlPoints3d)) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate DummyBezierControlPoints3D in ReshapeBezierSides')
!IF (SIZE(DummyBezierControlPoints3D).NE.SIZE(BezierControlPoints3D)) CALL abort(&
!__STAMP__&
!,'size of DummyBezierControlPoionts3D and BezierControlPoints3D not equal!')
!DummyBezierControlPoints3d=BezierControlPoints3D
!DEALLOCATE(BezierControlPoints3D)
!ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nTotalBCSides),STAT=ALLOCSTAT)
!BezierControlPoints3D=0.
!IF (ALLOCSTAT.NE.0) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate BezierControlPoints3D in ReshapeBezierSides')
!! SideSlabNormals
!ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nTotalSides))
!IF (.NOT.ALLOCATED(DummySideSlabNormals)) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate DummySideSlabNormals in ReshapeBezierSides')
!DummySideSlabNormals=SideSlabNormals
!DEALLOCATE(SideSlabNormals)
!ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate SideSlabNormals in ReshapeBezierSides')
!SideSlabNormals=0.
!! SideSlabIntervals
!ALLOCATE(DummySideSlabIntervals(1:6,1:nTotalSides))
!IF (.NOT.ALLOCATED(DummySideSlabIntervals)) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate DummySideSlabIntervals in ReshapeBezierSides')
!DummySideSlabIntervals=SideSlabIntervals
!DEALLOCATE(SideSlabIntervals)
!ALLOCATE(SideSlabIntervals(1:6,1:nTotalBCSides),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate ElemIndex in ReshapeBezierSides')
!SideSlabIntervals=0.
!! BoundingBoxIsEmpty
!ALLOCATE(DummyBoundingBoxIsEmpty(1:nTotalSides))
!IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate DummyBoundingBoxIsEmpty in ReshapeBezierSides')
!DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
!DEALLOCATE(BoundingBoxIsEmpty)
!ALLOCATE(BoundingBoxIsEmpty(1:nTotalBCSides),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!__STAMP__& !wunderschoen!!!
!,'Could not allocate BoundingBoxIsEmpty in ReshapeBezierSides')
!BoundingBoxIsEmpty=.FALSE.
!! side type
!ALLOCATE(DummySideType(1:nOldBCSides))
!IF (.NOT.ALLOCATED(DummySideType)) CALL abort(&
!    __STAMP__&
! ,'Could not allocate DummySideType in ReshapeBezierSides')
!DummySideType=SideType
!DEALLOCATE(SideType)
!ALLOCATE(SideType(1:nTotalBCSides),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!    __STAMP__&
! ,'Could not reallocate SideType in ReshapeBezierSides')
!SideType=-1

!ALLOCATE(DummySideDistance(1:nOldBCSides))
!IF (.NOT.ALLOCATED(DummySideDistance)) CALL abort(&
!    __STAMP__&
! ,'Could not allocate DummySideDistance in ReshapeBezierSides')
!DummySideDistance=SideDistance
!DEALLOCATE(SideDistance)
!ALLOCATE(SideDistance(1:nTotalBCSides),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!    __STAMP__&
! ,'Could not reallocate SideDistance in ReshapeBezierSides')
!SideDistance=0.

!ALLOCATE(DummySideNormVec(1:3,1:nOldBCSides))
!IF (.NOT.ALLOCATED(DummySideNormVec)) CALL abort(&
!    __STAMP__&
! ,'Could not allocate DummySideNormVec in ReshapeBezierSides')
!DummySideNormVec=SideNormVec
!DEALLOCATE(SideNormVec)
!ALLOCATE(SideNormVec(1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!    __STAMP__&
! ,'Could not reallocate SideNormVec in ReshapeBezierSides')
!SideNormVec=0.


!BCInc=0
!!DO iSide=1,nSides
!newBCSideID=0
!DO iSide=1,nBCSides
!  newBCSideID=newBCSideID+1
!  BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
!  SideSlabNormals          (1:3,1:3,          newBCSideID) =DummySideSlabNormals         (1:3,1:3,           iSide)
!  SideSlabIntervals       (1:6,              newBCSideID) =DummySideSlabIntervals      (1:6,               iSide)
!  BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
!  SideType(newBCSideID) = DummySideType(newBCSideID)
!  SideDistance(newBCSideID) = DummySideDistance(newBCSideID)
!  SideNormVec(1:3,newBCSideID) = DummySideNormVec(1:3,newBCSideID)
!END DO ! iSide

!DO iSide=nBCSides+1,nSides+nPartPeriodicSides
!  IF(BC(iSide).EQ.0) CYCLE
!  newBCSideID=newBCSideID+1
!  BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
!  SideSlabNormals          (1:3,1:3,          newBCSideID) =DummySideSlabNormals         (1:3,1:3,           iSide)
!  SideSlabIntervals       (1:6,              newBCSideID) =DummySideSlabIntervals      (1:6,               iSide)
!  BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
!  SideType(newBCSideID) = DummySideType(newBCSideID)
!  SideDistance(newBCSideID) = DummySideDistance(newBCSideID)
!  SideNormVec(1:3,newBCSideID) = DummySideNormVec(1:3,newBCSideID)
!END DO ! iSide

!DO iSide=nSides+nPartPeriodicSides+1,nTotalSides
!  newBCSideID=newBCSideID+1
!  BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
!  SideSlabNormals          (1:3,1:3,          newBCSideID) =DummySideSlabNormals         (1:3,1:3,           iSide)
!  SideSlabIntervals       (1:6,              newBCSideID) =DummySideSlabIntervals      (1:6,               iSide)
!  BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
!  SideType(newBCSideID) = DummySideType(newBCSideID)
!  SideDistance(newBCSideID) = DummySideDistance(newBCSideID)
!  SideNormVec(1:3,newBCSideID) = DummySideNormVec(1:3,newBCSideID)
!END DO ! iSide

!! create new mapping
!SDEALLOCATE(PartBCSideList)
!ALLOCATE(PartBCSideList(1:nTotalSides))
!PartBCSideList=-1

!newBCSideID=0
!DO iSide=1,nBCSides
!  newBCSideID=newBCSideID+1
!  PartBCSideList(iSide)=newBCSideID
!END DO

!DO iSide=nBCSides+1,nSides+nPartPeriodicSides
!  IF(BC(iSide).EQ.0) CYCLE
!  newBCSideID=newBCSideID+1
!  PartBCSideList(iSide)=newBCSideID
!END DO ! iSide

!DO iSide=nSides+nPartPeriodicSides+1,nTotalSides
!  newBCSideID=newBCSideID+1
!  PartBCSideList(iSide)=newBCSideID
!END DO

!! deallocate dummy buffer
!DEALLOCATE(DummyBezierControlPoints3D)
!DEALLOCATE(DummySideSlabNormals)
!DEALLOCATE(DummySideSlabIntervals)
!DEALLOCATE(DummyBoundingBoxIsEmpty)
!DEALLOCATE(DummySideType)
!DEALLOCATE(DummySideDistance)
!DEALLOCATE(DummySideNormVec)


!END SUBROUTINE ReShapeBezierSides



!SUBROUTINE MapElemToFIBGM()
!!----------------------------------------------------------------------------------------------------------------------------------!
!! here, the FIBGM range for each element is stored
!! short list for intersection tracking, longer list for ref mapping tracking
!!----------------------------------------------------------------------------------------------------------------------------------!
!! MODULES                                                                                                                          !
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,nTotalElems
!USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo
!USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!!----------------------------------------------------------------------------------------------------------------------------------!
!! insert modules here
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER           :: ALLOCSTAT,iElem,lastElem
!REAL              :: xmin,ymin,zmin,xmax,ymax,zmax
!INTEGER           :: BGMimax,BGMimin,BGMjmax,BGMjmin,BGMkmax,BGMkmin
!INTEGER           :: BGMCellXmax,BGMCellXmin,BGMCellYmax,BGMCellYmin,BGMCellZmax,BGMCellZmin
!!===================================================================================================================================

!!IF(.NOT.DoRefMapping) RETURN

!IF(DoRefMapping) THEN
!  LastElem=nTotalElems
!ELSE
!  LastElem=PP_nElems
!END IF

!ALLOCATE(GEO%ElemToFIBGM(1:6,1:LastElem),STAT=ALLOCSTAT )
!IF (ALLOCSTAT.NE.0) CALL abort(&
!__STAMP__&
!,'  Cannot allocate GEO%ElemToFIBGM!')

!! because I copy and past
!BGMimax=GEO%FIBGMimax
!BGMimin=GEO%FIBGMimin
!BGMjmax=GEO%FIBGMjmax
!BGMjmin=GEO%FIBGMjmin
!BGMkmax=GEO%FIBGMkmax
!BGMkmin=GEO%FIBGMkmin

!DO iElem=1,LastElem
!  xmin=HUGE(1.)
!  ymin=HUGE(1.)
!  zmin=HUGE(1.)
!  xmax=-HUGE(1.)
!  ymax=-HUGE(1.)
!  zmax=-HUGE(1.)
!  xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
!  xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
!  ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
!  ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
!  zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
!  zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
!  !--- find minimum and maximum BGM cell for current element
!  IF(GEO%nPeriodicVectors.EQ.0)THEN
!    BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
!    BGMCellXmax = MIN(BGMCellXmax,BGMimax)
!    BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
!    BGMCellXmin = MAX(BGMCellXmin,BGMimin)
!    BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
!    BGMCellYmax = MIN(BGMCellYmax,BGMjmax)
!    BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
!    BGMCellYmin = MAX(BGMCellYmin,BGMjmin)
!    BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
!    BGMCellZmax = MIN(BGMCellZmax,BGMkmax)
!    BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
!    BGMCellZmin = MAX(BGMCellZmin,BGMkmin)
!  ELSE
!    ! here fancy stuff, because element could be wide out of element range
!    BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
!    BGMCellXmax = MAX(MIN(BGMCellXmax,BGMimax),BGMimin)
!    BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
!    BGMCellXmin = MIN(MAX(BGMCellXmin,BGMimin),BGMimax)
!    BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
!    BGMCellYmax = MAX(MIN(BGMCellYmax,BGMjmax),BGMjmin)
!    BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
!    BGMCellYmin = MIN(MAX(BGMCellYmin,BGMjmin),BGMjmax)
!    BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
!    BGMCellZmax = MAX(MIN(BGMCellZmax,BGMkmax),BGMkmin)
!    BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
!    BGMCellZmin = MIN(MAX(BGMCellZmin,BGMkmin),BGMkmax)
!  END IF
!  GEO%ElemToFIBGM(1,iElem)=BGMCellXmin
!  GEO%ElemToFIBGM(3,iElem)=BGMCellYmin
!  GEO%ElemToFIBGM(5,iElem)=BGMCellZmin

!  GEO%ElemToFIBGM(2,iElem)=BGMCellXmax
!  GEO%ElemToFIBGM(4,iElem)=BGMCellYmax
!  GEO%ElemToFIBGM(6,iElem)=BGMCellZmax
!END DO ! iElem=1,nTotalElems

!END SUBROUTINE MapElemToFIBGM


!SUBROUTINE CheckIfElemCurved(IsCurved,XCL_NGeo)
!!===================================================================================================================================
!! check if element is curved
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Mesh_Vars   ,ONLY: NGeo,Vdm_CLNGeo1_CLNGeo
!USE MOD_ChangeBasis ,ONLY: changeBasis3D
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!REAL,INTENT(IN)      :: XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo)
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!LOGICAL,INTENT(OUT)  :: IsCurved
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL                 :: XCL_NGeo1(1:3,0:1,0:1,0:1)
!REAL                 :: XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0:NGeo)
!INTEGER              :: NGeo3
!!===================================================================================================================================
!
!IsCurved=.FALSE.
!
!! fill dummy
!XCL_NGeo1(1:3,0,0,0) = XCL_NGeo(1:3, 0  , 0  , 0  )
!XCL_NGeo1(1:3,1,0,0) = XCL_NGeo(1:3,NGeo, 0  , 0  )
!XCL_NGeo1(1:3,0,1,0) = XCL_NGeo(1:3, 0  ,NGeo, 0  )
!XCL_NGeo1(1:3,1,1,0) = XCL_NGeo(1:3,NGeo,NGeo, 0  )
!XCL_NGeo1(1:3,0,0,1) = XCL_NGeo(1:3, 0  , 0  ,NGeo)
!XCL_NGeo1(1:3,1,0,1) = XCL_NGeo(1:3,NGeo, 0  ,NGeo)
!XCL_NGeo1(1:3,0,1,1) = XCL_NGeo(1:3, 0  ,NGeo,NGeo)
!XCL_NGeo1(1:3,1,1,1) = XCL_NGeo(1:3,NGeo,NGeo,NGeo)
!
!CALL ChangeBasis3D(3,1,NGeo,Vdm_CLNGeo1_CLNGeo,XCL_NGeo1,XCL_NGeoNew)
!NGeo3=(NGeo+1)*(NGeo+1)*(NGeo+1)
!
!! check 3D points
!CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeo,IsCurved)
!
!IF(.NOT.IsCurved)THEN
!  ! set all elem sides to blabla
!END IF
!
!END SUBROUTINE CheckIfElemCurved


!SUBROUTINE GetBCElemMap()
!!===================================================================================================================================
!! 1) Add BC and Halo-BC sides in halo_eps distance to a certain element
!! 2) build epsOneCell for each element
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
!USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo,nSides,NGeo,nBCSides,sJ,nElems
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalSides,IsTracingBCElem,nTotalElems,nTotalBCElems
!USE MOD_Particle_Mesh_Vars     ,ONLY: TracingBCInnerSides,TracingBCTotalSides
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide,BCElem,PartSideToElem,PartBCSideList,GEO
!USE MOD_Particle_Mesh_Vars     ,ONLY: RefMappingEps,epsOneCell
!USE MOD_Particle_Surfaces_Vars ,ONLY: sVdm_Bezier
!USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps,halo_eps2
!USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
!USE MOD_Analyze_Vars           ,ONLY: CalcMeshInfo
!#if USE_MPI
!USE MOD_Mesh_Vars              ,ONLY: BC
!#endif /*USE_MPI*/
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                                  :: iElem,firstBezierPoint,lastBezierPoint
!INTEGER                                  :: iSide,p,q,SideID,ilocSide,BCSideID2,BCSideID
!INTEGER                                  :: nSideCount, s,r
!INTEGER,ALLOCATABLE                      :: SideIndex(:)
!REAL,DIMENSION(1:3)                      :: v1,NodeX,Vec1
!REAL,DIMENSION(1:3,0:NGeo,0:NGeo)        :: xNodes
!INTEGER                                  :: nLoop,iTest,nTest
!REAL                                     :: scaleJ, Distance ,maxScaleJ,dx,dy,dz
!LOGICAL                                  :: fullMesh, leave
!!===================================================================================================================================
!ALLOCATE(IsTracingBCElem(nTotalElems))
!IsTracingBCElem=.FALSE.
!IF(DoRefMapping)THEN
!  ALLOCATE(TracingBCInnerSides(nTotalElems))
!  TracingBCInnerSides=0
!  ALLOCATE(TracingBCTotalSides(nTotalElems))
!  TracingBCTotalSides=0
!END IF

!! decide if element:
!! DoRefMapping=T
!! a) HAS own bc faces
!! b) HAS bc-face in halo_eps distance
!! DoRefMapping=F
!! a) HAS own bc faces
!IF(DoRefMapping)THEN
!  ! mark elements as bc element if they have a local-BC side
!  nTotalBCElems=0
!  DO iElem=1,nTotalElems
!    DO ilocSide=1,6
!      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!      IF (SideID.LE.0) CYCLE
!      IF((SideID.LE.nBCSides).OR.(SideID.GT.nSides))THEN
!        IF(.NOT.IsTracingBCElem(iElem))THEN
!          IsTracingBCElem(iElem)=.TRUE.
!          nTotalBCElems=nTotalBCElems+1
!        END IF ! count only single
!      END IF
!    END DO ! ilocSide
!  END DO ! iElem

!  ! for simplifications
!  ! get distance of diagonal of mesh
!  V1(1) = GEO%xmaxglob-GEO%xminglob
!  V1(2) = GEO%ymaxglob-GEO%yminglob
!  V1(3) = GEO%zmaxglob-GEO%zminglob
!  Distance=DOT_PRODUCT(V1,V1)
!  fullMesh=.FALSE.
!  ! build list with elements in halo-eps vicinity around bc-elements
!  IF(Distance.LE.halo_eps2) fullMesh=.TRUE.
!  ! allocate the types for the element to bc-side mapping
!  ALLOCATE( BCElem(1:nTotalElems) )
!  ALLOCATE( SideIndex(1:nTotalSides) )
!  ! for fullMesh, each element requires ALL BC faces
!  IF(fullMesh)THEN
!    DO iElem=1,nTotalElems
!      ! mark my sides
!      BCElem(iElem)%nInnerSides=0
!      DO ilocSide=1,6
!        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!        IF(SideID.LE.0) CYCLE
!        IF(PartBCSideList(SideID).EQ.-1) CYCLE
!        BCElem(iElem)%nInnerSides = BCElem(iElem)%nInnerSides+1
!      END DO ! ilocSide=1,6
!      BCElem(iElem)%lastSide=BCElem(iElem)%nInnerSides
!      ! loop over all sides, exclusive of own sides
!      SideIndex=0
!      DO iSide=1,nTotalSides
!        ! only bc sides
!        BCSideID  =PartBCSideList(iSide)
!        IF(BCSideID.EQ.-1) CYCLE
!        ! ignore sides of the same element
!        IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE
!        IF(SideIndex(iSide).EQ.0)THEN
!          BCElem(iElem)%lastSide=BCElem(iElem)%lastSide+1
!          SideIndex(iSide)=BCElem(iElem)%lastSide
!        END IF
!      END DO ! iSide=1,nTotalSides
!      IF(BCElem(iElem)%lastSide.EQ.0) CYCLE
!      ! set true, only required for elements without an own bc side
!      IsTracingBCElem(iElem)=.TRUE.
!      ! allocate complete side list
!      ALLOCATE( BCElem(iElem)%BCSideID(BCElem(iElem)%lastSide) )
!      ! 1) inner sides
!      nSideCount=0
!      IF(BCElem(iElem)%nInnerSides.GT.0)THEN
!        DO ilocSide=1,6
!          SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!          IF(SideID.LE.0) CYCLE
!          BCSideID=PartBCSideList(SideID)
!          IF(BCSideID.LE.0) CYCLE
!          nSideCount=nSideCount+1
!          BCElem(iElem)%BCSideID(nSideCount)=SideID
!        END DO ! ilocSide
!      END IF ! nInnerSides.GT.0
!      ! 2) outer sides
!      DO iSide=1,nTotalSides
!        IF(SideIndex(iSide).GT.0)THEN
!          nSideCount=nSideCount+1
!          BCElem(iElem)%BCSideID(nSideCount)=iSide !iSide
!        END IF
!      END DO  ! iSide
!    END DO ! iElem=1,nTotalElems
!  ELSE ! .NOT. fullMesh
!    ! each element requires only the sides in its halo region
!    DO iElem=1,nTotalElems
!      ! mark my sides
!      BCElem(iElem)%nInnerSides=0
!      DO ilocSide=1,6
!        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!        IF(SideID.LE.0) CYCLE
!        IF(PartBCSideList(SideID).EQ.-1) CYCLE
!        BCElem(iElem)%nInnerSides = BCElem(iElem)%nInnerSides+1
!      END DO ! ilocSide=1,6
!      BCElem(iElem)%lastSide=BCElem(iElem)%nInnerSides
!      ! loop over all sides, to reduce required storage, if a side is marked once,
!      ! it does not have to be checked for further sides
!      SideIndex=0
!      DO ilocSide=1,6
!        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!        BCSideID2=SideID
!        IF(SideID.GT.0) BCSideID2=PartBCSideList(SideID)
!        IF (BCSideID2.GT.0) THEN
!          xNodes(:,:,:)=BezierControlPoints3D(:,:,:,PartBCSideList(SideID))
!          SELECT CASE(ilocSide)
!          CASE(XI_MINUS,XI_PLUS)
!            firstBezierPoint=0
!            lastBezierPoint=NGeo
!          CASE DEFAULT
!            firstBezierPoint=1
!            lastBezierPoint=NGeo-1
!          END SELECT
!        ELSE
!          SELECT CASE(ilocSide)
!          CASE(XI_MINUS)
!            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,iElem),xNodes(:,:,:))
!            firstBezierPoint=0
!            lastBezierPoint=NGeo
!          CASE(XI_PLUS)
!            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,iElem),xNodes(:,:,:))
!            firstBezierPoint=0
!            lastBezierPoint=NGeo
!          CASE(ETA_MINUS)
!            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,iElem),xNodes(:,:,:))
!            firstBezierPoint=1
!            lastBezierPoint=NGeo-1
!          CASE(ETA_PLUS)
!            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,iElem),xNodes(:,:,:))
!            firstBezierPoint=1
!            lastBezierPoint=NGeo-1
!          CASE(ZETA_MINUS)
!            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,iElem),xNodes(:,:,:))
!            firstBezierPoint=1
!            lastBezierPoint=NGeo-1
!          CASE(ZETA_PLUS)
!            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,iElem),xNodes(:,:,:))
!            firstBezierPoint=1
!            lastBezierPoint=NGeo-1
!          END SELECT
!        END IF
!        DO iSide=1,nTotalSides
!          ! only bc sides
!          BCSideID  =PartBCSideList(iSide)
!          IF(BCSideID.EQ.-1) CYCLE
!          ! ignore sides of the same element
!          IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE
!          IF(SideIndex(iSide).EQ.0)THEN
!            leave=.FALSE.
!            nTest=1
!            DO iTest=1,nTest
!              Vec1=0.
!              ! all points of bc side
!              DO q=firstBezierPoint,lastBezierPoint
!                DO p=firstBezierPoint,lastBezierPoint
!                  NodeX(:) = BezierControlPoints3D(:,p,q,BCSideID)+Vec1
!                  !all nodes of current side
!                  DO s=firstBezierPoint,lastBezierPoint
!                    DO r=firstBezierPoint,lastBezierPoint
!                      dX=ABS(xNodes(1,r,s)-NodeX(1))
!                      IF(dX.GT.halo_eps) CYCLE
!                      dY=ABS(xNodes(2,r,s)-NodeX(2))
!                      IF(dY.GT.halo_eps) CYCLE
!                      dZ=ABS(xNodes(3,r,s)-NodeX(3))
!                      IF(dZ.GT.halo_eps) CYCLE
!                      IF(SQRT(dX*dX+dY*dY+dZ*dZ).LE.halo_eps)THEN
!                        IF(SideIndex(iSide).EQ.0)THEN
!                          BCElem(iElem)%lastSide=BCElem(iElem)%lastSide+1
!                          SideIndex(iSide)=BCElem(iElem)%lastSide
!                          leave=.TRUE.
!                          EXIT
!                        END IF
!                      END IF
!                    END DO ! r
!                    IF(leave) EXIT
!                  END DO ! s
!                  IF(leave) EXIT
!                END DO ! p
!                IF(leave) EXIT
!              END DO ! q
!              IF(leave) EXIT
!            END DO ! iTest=1,nTest
!          END IF ! SideIndex(iSide).EQ.0
!        END DO ! iSide=1,nTotalSides
!      END DO ! ilocSide=1,6
!      IF(BCElem(iElem)%lastSide.EQ.0) CYCLE
!      ! set true, only required for elements without an own bc side
!      IsTracingBCElem(iElem)=.TRUE.
!      ! allocate complete side list
!      ALLOCATE( BCElem(iElem)%BCSideID(BCElem(iElem)%lastSide) )
!      ! 1) inner sides
!      nSideCount=0
!      IF(BCElem(iElem)%nInnerSides.GT.0)THEN
!        DO ilocSide=1,6
!          SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!          IF(SideID.LE.0) CYCLE
!          BCSideID=PartBCSideList(SideID)
!          IF(BCSideID.LE.0) CYCLE
!          nSideCount=nSideCount+1
!          BCElem(iElem)%BCSideID(nSideCount)=SideID
!        END DO ! ilocSide
!      END IF ! nInnerSides.GT.0
!      ! 2) outer sides
!      DO iSide=1,nTotalSides
!        IF(SideIndex(iSide).GT.0)THEN
!          nSideCount=nSideCount+1
!          BCElem(iElem)%BCSideID(nSideCount)=iSide !iSide
!        END IF
!      END DO  ! iSide
!    END DO ! iElem=1,nTotalElems
!  END IF ! fullMesh
!ELSE ! .NOT.DoRefMapping
!  ! tracing
!  ! mark only elements with bc-side
!  nTotalBCElems=0
!  DO iElem=1,nTotalElems
!    DO ilocSide=1,6
!      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!      IF (SideID.LE.0) CYCLE
!      IF(SideID.LE.nBCSides)THEN ! non-halo elements
!        IF(.NOT.IsTracingBCElem(iElem))THEN
!          IsTracingBCElem(iElem)=.TRUE.
!          nTotalBCElems=nTotalBCElems+1
!        END IF ! count only single
!      END IF
!#if USE_MPI
!      IF(SideID.GT.nSides)THEN ! halo elements
!        IF(BC(SideID).NE.0)THEN
!          IF(.NOT.IsTracingBCElem(iElem))THEN
!            IsTracingBCElem(iElem)=.TRUE.
!            nTotalBCElems=nTotalBCElems+1
!          END IF ! count only single
!        END IF
!      END IF ! SideID.GT.nSides
!#endif
!    END DO ! ilocSide
!  END DO ! iElem
!END IF

!IF(DoRefMapping)THEN
!  DO iElem=1,nTotalElems
!    TracingBCInnerSides(iElem) = BCElem(iElem)%nInnerSides
!    TracingBCTotalSides(iElem) = BCElem(iElem)%lastSide
!  END DO ! iElem

!  IF(CalcMeshInfo)THEN
!    CALL AddToElemData(ElementOut,'TracingBCInnerSides',IntArray=TracingBCInnerSides(1:nElems))
!    CALL AddToElemData(ElementOut,'TracingBCTotalSides',IntArray=TracingBCTotalSides(1:nElems))
!  END IF
!END IF

!IF(CalcMeshInfo)THEN
!  CALL AddToElemData(ElementOut,'IsTracingBCElem'    ,LogArray=IsTracingBCElem(    1:nElems))
!END IF

!! finally, build epsonecell per element
!IF(DoRefMapping)THEN
!  ALLOCATE(epsOneCell(1:nTotalElems))
!ELSE
!  ALLOCATE(epsOneCell(1:PP_nElems))
!END IF
!epsOneCell=0.

!nLoop=nTotalElems
!IF(.NOT.DoRefMapping) nLoop=PP_nElems
!maxScaleJ=0.
!DO iElem=1,PP_nElems
!  scaleJ=MAXVAL(sJ(:,:,:,iElem))/MINVAL(sJ(:,:,:,iElem))
!  epsOneCell(iElem)=1.0+SQRT(3.0*scaleJ*RefMappingEps)
!  maxScaleJ=MAX(scaleJ,maxScaleJ)
!END DO ! iElem=1,nLoop
!DO iElem=PP_nElems+1,nLoop
!  epsOneCell(iElem)=1.0+SQRT(maxScaleJ*RefMappingEps)
!END DO ! iElem=1,nLoop

!IF(CalcMeshInfo)THEN
!  CALL AddToElemData(ElementOut,'epsOneCell',RealArray=epsOneCell(1:nElems))
!END IF

!END SUBROUTINE GetBCElemMap


!SUBROUTINE GetShapeFunctionBCElems()
!!===================================================================================================================================
!! Identify all elements that are close to boundaries, where the deposition via shape function would cause the shape function sphere
!! to be truncated by the boundary. In this case, a local deposition is used in that cell for "inner" parts, i.e., shape functions
!! that extend into the element by exterior particle shape functions are still deposited via the shape function.
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
!USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo,NGeo,nElems,BC
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalSides,IsLocalDepositionBCElem,nTotalElems
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide,PartSideToElem,PartBCSideList,SidePeriodicType
!USE MOD_Particle_Surfaces_Vars ,ONLY: sVdm_Bezier
!USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
!USE MOD_PICDepo_Vars           ,ONLY: r_sf,DepositionType,sf1d_dir
!USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps
!USE MOD_Mesh_Vars              ,ONLY: BoundaryType
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                                  :: iElem,firstBezierPoint,lastBezierPoint
!INTEGER                                  :: iSide,p,q,SideID,ilocSide,BCSideID2,BCSideID
!INTEGER                                  :: s,r
!INTEGER,ALLOCATABLE                      :: SideIndex(:)
!REAL,DIMENSION(1:3)                      :: NodeX
!REAL,DIMENSION(1:3,0:NGeo,0:NGeo)        :: xNodes
!REAL                                     :: dx,dy,dz
!LOGICAL                                  :: leave,BCElemSF
!!===================================================================================================================================
!! allocate for local elements + halo elements
!ALLOCATE(IsLocalDepositionBCElem(nTotalElems))
!IsLocalDepositionBCElem=.FALSE.
!! Only add local elements to element list
!CALL AddToElemData(ElementOut,'IsLocalDepositionBCElem',LogArray=IsLocalDepositionBCElem(1:nElems))

!! Auxiliary integer array for marking sides (set equal to 1 if the side is within shape function radius of considered local side)
!ALLOCATE(SideIndex(1:nTotalSides))

!! =============================
!! Workflow:
!!  0.  Sanity Check: halo distance must be equal to or larger than the shape function radius
!!  1.  Loop over all elements (including halo elements)
!!  1.1  Mark rank-local elements with BC sides without calculating the distance between the nodes (more efficient
!!  1.2  Check distance (loop all local sides and measure the distance to all nTotalSides)
!!  1.3  Loop over all sides (including halo sides) and compare distance to the current iLocSide of the element
!!==============================

!! 0.   Check halo distance vs. shape function radius, because the halo region is used for checking the shape function deposition
!IF(halo_eps.LT.r_sf)THEN
!  SWRITE(UNIT_StdOut,'(132("*"))')
!  SWRITE(UNIT_StdOut,'(A)') ' Warning in GetShapeFunctionBCElems: halo_eps is less than r_sh, which may result in wrong '//&
!                            'deposition elements.\n Check IsLocalDepositionBCElem in state file!'
!  SWRITE(UNIT_StdOut,'(A,ES25.14E3)') '  halo_eps : ',halo_eps
!  SWRITE(UNIT_StdOut,'(A,ES25.14E3)') '  r_sf     : ',r_sf
!  SWRITE(UNIT_StdOut,'(A)') ' Consider increasing the halo velocity to remove this warning.'
!  SWRITE(UNIT_StdOut,'(132("*"))')
!END IF



!! 1.  Loop over all elements (including halo elements)
!DO iElem=1,nTotalElems
!  BCElemSF=.FALSE.
!  ! 1.1  Mark rank-local elements with BC sides without calculating the distance between the nodes (more efficient)
!  DO ilocSide=1,6
!    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(SideID.LE.0)                             CYCLE ! Skip sides?
!    IF(PartBCSideList(SideID).EQ.-1)            CYCLE ! Skip non-BC sides
!    IF(SidePeriodicType(SideID).NE.0)           CYCLE ! Skip periodic-BC sides
!    IF(BoundaryType(BC(SideID),BC_TYPE).EQ.100) CYCLE ! Skip inner-BC sides (must be labelled with 100)


!    ! Skip BC sides for shape_function_2d
!    IF(TRIM(DepositionType).EQ.'shape_function_2d')THEN
!      ASSOCIATE ( &
!            x1 => BezierControlPoints3D(sf1d_dir , 0    , 0    , PartBCSideList(SideID))   , &
!            x2 => BezierControlPoints3D(sf1d_dir , 0    , NGeo , PartBCSideList(SideID))   , &
!            x3 => BezierControlPoints3D(sf1d_dir , NGeo , 0    , PartBCSideList(SideID))   , &
!            x4 => BezierControlPoints3D(sf1d_dir , NGeo , NGeo , PartBCSideList(SideID)) )
!        ! Check if all corner points are equal is the "sf1d_dir" direction: Skip this side if true
!        IF((ALMOSTEQUALRELATIVE(x1,x2,1e-6).AND.&
!            ALMOSTEQUALRELATIVE(x1,x3,1e-6).AND.&
!            ALMOSTEQUALRELATIVE(x1,x4,1e-6))) CYCLE
!      END ASSOCIATE
!    END IF
!    IsLocalDepositionBCElem(iElem)=.TRUE.
!    EXIT
!  END DO ! ilocSide=1,6

!  IF(IsLocalDepositionBCElem(iElem)) CYCLE ! finished: next element



!  ! 1.2  Check distance (loop all local sides and measure the distance to all nTotalSides)
!  ! loop over all sides, to reduce required storage, if a side is marked once, it does not have to be checked for further sides
!  SideIndex=0
!  DO ilocSide=1,6
!    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(SideID.GT.0)THEN
!      BCSideID2=PartBCSideList(SideID)
!    ELSE
!      BCSideID2=SideID
!    END IF

!    IF(BCSideID2.GT.0) THEN
!      xNodes(:,:,:)=BezierControlPoints3D(:,:,:,PartBCSideList(SideID))
!      SELECT CASE(ilocSide)
!      CASE(XI_MINUS,XI_PLUS)
!        firstBezierPoint=0
!        lastBezierPoint=NGeo
!      CASE DEFAULT
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      END SELECT
!    ELSE
!      SELECT CASE(ilocSide)
!      CASE(XI_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=0
!        lastBezierPoint=NGeo
!      CASE(XI_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=0
!        lastBezierPoint=NGeo
!      CASE(ETA_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      CASE(ETA_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      CASE(ZETA_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      CASE(ZETA_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      END SELECT
!    END IF


!    ! 1.3  Loop over all sides (including halo sides) and compare distance to the current iLocSide of the element
!    DO iSide=1,nTotalSides
!      BCSideID=PartBCSideList(iSide) ! only bc sides
!      IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE ! Skip sides of the same element
!      IF(BCSideID.EQ.-1)                             CYCLE ! Skip non-BC sides
!      IF(SidePeriodicType(iSide).NE.0)               CYCLE ! Skip periodic sides. Note that side = iSide and not BCSideID
!      IF(BoundaryType(BC(iSide),BC_TYPE).EQ.100)     CYCLE ! Skip inner-BC sides (must be labelled with 100)

!      ! Skip BC sides for shape_function_2d
!      IF(TRIM(DepositionType).EQ.'shape_function_2d')THEN
!        ASSOCIATE ( &
!              x1 => BezierControlPoints3D(sf1d_dir , 0    , 0    , BCSideID)   , &
!              x2 => BezierControlPoints3D(sf1d_dir , 0    , NGeo , BCSideID)   , &
!              x3 => BezierControlPoints3D(sf1d_dir , NGeo , 0    , BCSideID)   , &
!              x4 => BezierControlPoints3D(sf1d_dir , NGeo , NGeo , BCSideID) )
!          ! Check if all corner points are equal is the "sf1d_dir" direction: Skip this side if true
!          IF((ALMOSTEQUALRELATIVE(x1,x2,1e-6).AND.&
!              ALMOSTEQUALRELATIVE(x1,x3,1e-6).AND.&
!              ALMOSTEQUALRELATIVE(x1,x4,1e-6))) CYCLE
!        END ASSOCIATE
!      END IF

!      IF(SideIndex(iSide).EQ.0)THEN
!        leave=.FALSE.
!        ! all points of bc side
!        DO q=firstBezierPoint,lastBezierPoint
!          DO p=firstBezierPoint,lastBezierPoint
!            NodeX(:) = BezierControlPoints3D(:,p,q,BCSideID)
!            !all nodes of current side
!            DO s=firstBezierPoint,lastBezierPoint
!              DO r=firstBezierPoint,lastBezierPoint
!                dX=ABS(xNodes(1,r,s)-NodeX(1))
!                IF(dX.GT.r_sf) CYCLE
!                dY=ABS(xNodes(2,r,s)-NodeX(2))
!                IF(dY.GT.r_sf) CYCLE
!                dZ=ABS(xNodes(3,r,s)-NodeX(3))
!                IF(dZ.GT.r_sf) CYCLE
!                IF(SQRT(dX*dX+dY*dY+dZ*dZ).LE.r_sf)THEN
!                  IF(SideIndex(iSide).EQ.0)THEN
!                    BCElemSF=.TRUE.
!                    SideIndex(iSide)=1 ! mark with number .NE. 0
!                    leave=.TRUE.
!                    EXIT
!                  END IF
!                END IF
!              END DO ! r
!              IF(leave) EXIT
!            END DO ! s
!            IF(leave) EXIT
!          END DO ! p
!          IF(leave) EXIT
!        END DO ! q
!        IF(leave) EXIT
!      END IF ! SideIndex(iSide).EQ.0
!    END DO ! iSide=1,nTotalSides
!  END DO ! ilocSide=1,6

!  ! set true, only required for elements without an own bc side
!  IF(BCElemSF) IsLocalDepositionBCElem(iElem)=.TRUE.
!END DO ! iElem=1,nTotalElems




!END SUBROUTINE GetShapeFunctionBCElems


!SUBROUTINE GetShapeFunctionBCElems_OLD()
!!===================================================================================================================================
!! Identify all elements that are close to boundaries, where the deposition via shape function would cause the shape function sphere
!! to be truncated by the boundary. In this case, a local deposition is used in that cell for "inner" parts, i.e., shape functions
!! that extend into the element by exterior particle shape functions are still deposited via the shape function.
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
!USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo,NGeo,nElems
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalSides,IsLocalDepositionBCElem,nTotalElems
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide,PartSideToElem,PartBCSideList,SidePeriodicType
!USE MOD_Particle_Surfaces_Vars ,ONLY: sVdm_Bezier
!USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
!USE MOD_PICDepo_Vars           ,ONLY: r_sf,DepositionType,sf1d_dir
!USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                                  :: iElem,firstBezierPoint,lastBezierPoint
!INTEGER                                  :: iSide,p,q,SideID,ilocSide,BCSideID2,BCSideID
!INTEGER                                  :: s,r
!INTEGER,ALLOCATABLE                      :: SideIndex(:)
!REAL,DIMENSION(1:3)                      :: NodeX,Vec1
!REAL,DIMENSION(1:3,0:NGeo,0:NGeo)        :: xNodes
!REAL                                     :: dx,dy,dz
!LOGICAL                                  :: leave,BCElemSF
!!===================================================================================================================================
!! allocate for local elements + halo elements
!ALLOCATE(IsLocalDepositionBCElem(nTotalElems))
!IsLocalDepositionBCElem=.FALSE.
!! only add local elements to element list
!CALL AddToElemData(ElementOut,'IsLocalDepositionBCElem'    ,LogArray=IsLocalDepositionBCElem(    1:nElems))

!! =============================
!! Workflow:
!!
!!  0.  Check halo distance vs. shape function radius
!!  1.  Check local BC sides
!!  2.  Check halo BC sides: each element requires only the sides in its halo region
!!==============================

!! 0.   Check halo distance vs. shape function radius, because the halo region is used for checking the shape function deposition
!IF(halo_eps.LT.r_sf)THEN
!  SWRITE(UNIT_StdOut,'(132("*"))')
!  SWRITE(UNIT_StdOut,'(A)') ' Warning in GetShapeFunctionBCElems: halo_eps is less than r_sh, which may result in wrong '//&
!                            'deposition elements.\n Check IsLocalDepositionBCElem in state file!'
!  SWRITE(UNIT_StdOut,'(A,ES25.14E3)') '  halo_eps : ',halo_eps
!  SWRITE(UNIT_StdOut,'(A,ES25.14E3)') '  r_sf     : ',r_sf
!  SWRITE(UNIT_StdOut,'(A)') ' Consider increasing the halo velocity to remove this warning.'
!  SWRITE(UNIT_StdOut,'(132("*"))')
!END IF

!!    ! 1.   Check local BC sides:  mark elements as bc element if they have a local-BC side (skip periodic sides)
!!    DO iElem=1,nTotalElems
!!      DO ilocSide=1,6
!!        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!!        IF (SideID.LE.0)                               CYCLE ! Skip ?
!!        IF((SideID.LE.nBCSides).OR.(SideID.GT.nSides)) THEN  ! Don't skip ? and ?
!!          IF(SidePeriodicType(SideID).NE.0)            CYCLE ! skip periodic sides
!!          ! Skip BC sides for shape_function_2d
!!          IF(TRIM(DepositionType).EQ.'shape_function_2d')THEN
!!            ASSOCIATE ( &
!!                  x1 => BezierControlPoints3D(sf1d_dir , 0    , 0    , PartBCSideList(SideID))   , &
!!                  x2 => BezierControlPoints3D(sf1d_dir , 0    , NGeo , PartBCSideList(SideID))   , &
!!                  x3 => BezierControlPoints3D(sf1d_dir , NGeo , 0    , PartBCSideList(SideID))   , &
!!                  x4 => BezierControlPoints3D(sf1d_dir , NGeo , NGeo , PartBCSideList(SideID)) )
!!              ! Check if all corner points are equal is the "sf1d_dir" direction: Skip this side if true
!!              IF((ALMOSTEQUALRELATIVE(x1,x2,1e-6).AND.&
!!                  ALMOSTEQUALRELATIVE(x1,x3,1e-6).AND.&
!!                  ALMOSTEQUALRELATIVE(x1,x4,1e-6))) CYCLE
!!            END ASSOCIATE
!!          END IF
!!          IsLocalDepositionBCElem(iElem)=.TRUE.
!!        END IF
!!      END DO ! ilocSide
!!    END DO ! iElem

!! 2.   Check halo BC sides: each element requires only the sides in its halo region
!ALLOCATE( SideIndex(1:nTotalSides) )
!DO iElem=1,nTotalElems

!  !   IF(IsLocalDepositionBCElem(iElem)) CYCLE ! identified in previous step

!  ! mark my sides
!  BCElemSF=.FALSE.
!  DO ilocSide=1,6
!    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(SideID.LE.0)                   CYCLE ! Skip ?
!    IF(PartBCSideList(SideID).EQ.-1)  CYCLE ! Skip non-BC sides
!    IF(SidePeriodicType(SideID).NE.0) CYCLE ! Skip periodic-BC sides
!    !IF(SideID.GT.nBCSides)            CYCLE ! Skip non-BC sides -> already done in PartBCSideList(SideID).EQ.-1 ?

!    ! Skip BC sides for shape_function_2d
!    IF(TRIM(DepositionType).EQ.'shape_function_2d')THEN
!      ASSOCIATE ( &
!            x1 => BezierControlPoints3D(sf1d_dir , 0    , 0    , PartBCSideList(SideID))   , &
!            x2 => BezierControlPoints3D(sf1d_dir , 0    , NGeo , PartBCSideList(SideID))   , &
!            x3 => BezierControlPoints3D(sf1d_dir , NGeo , 0    , PartBCSideList(SideID))   , &
!            x4 => BezierControlPoints3D(sf1d_dir , NGeo , NGeo , PartBCSideList(SideID)) )
!        ! Check if all corner points are equal is the "sf1d_dir" direction: Skip this side if true
!        IF((ALMOSTEQUALRELATIVE(x1,x2,1e-6).AND.&
!            ALMOSTEQUALRELATIVE(x1,x3,1e-6).AND.&
!            ALMOSTEQUALRELATIVE(x1,x4,1e-6))) CYCLE
!      END ASSOCIATE
!    END IF
!    IsLocalDepositionBCElem(iElem)=.TRUE.
!    EXIT
!  END DO ! ilocSide=1,6

!  IF(IsLocalDepositionBCElem(iElem)) CYCLE ! finished: next element



!  ! 3.  Check distance
!  ! loop over all sides, to reduce required storage, if a side is marked once,
!  ! it does not have to be checked for further sides
!  SideIndex=0
!  DO ilocSide=1,6
!    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    !BCSideID2=SideID
!    !IF(SideID.GT.0) BCSideID2=PartBCSideList(SideID)
!    IF(SideID.GT.0)THEN
!      BCSideID2=PartBCSideList(SideID)
!    ELSE
!      BCSideID2=SideID
!    END IF

!    IF(BCSideID2.GT.0) THEN
!      xNodes(:,:,:)=BezierControlPoints3D(:,:,:,PartBCSideList(SideID))
!      SELECT CASE(ilocSide)
!      CASE(XI_MINUS,XI_PLUS)
!        firstBezierPoint=0
!        lastBezierPoint=NGeo
!      CASE DEFAULT
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      END SELECT
!    ELSE
!      SELECT CASE(ilocSide)
!      CASE(XI_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=0
!        lastBezierPoint=NGeo
!      CASE(XI_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=0
!        lastBezierPoint=NGeo
!      CASE(ETA_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      CASE(ETA_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      CASE(ZETA_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      CASE(ZETA_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,iElem),xNodes(:,:,:))
!        firstBezierPoint=1
!        lastBezierPoint=NGeo-1
!      END SELECT
!    END IF
!    DO iSide=1,nTotalSides
!      BCSideID=PartBCSideList(iSide) ! only bc sides
!      IF(BCSideID.EQ.-1)                             CYCLE ! Skip non-BC sides
!      IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE ! Ignore sides of the same element
!      IF(SidePeriodicType(iSide).NE.0)               CYCLE ! Skip periodic sides. Note that side = iSide and not BCSideID

!      ! Skip BC sides for shape_function_2d
!      IF(TRIM(DepositionType).EQ.'shape_function_2d')THEN
!        ASSOCIATE ( &
!              x1 => BezierControlPoints3D(sf1d_dir , 0    , 0    , BCSideID)   , &
!              x2 => BezierControlPoints3D(sf1d_dir , 0    , NGeo , BCSideID)   , &
!              x3 => BezierControlPoints3D(sf1d_dir , NGeo , 0    , BCSideID)   , &
!              x4 => BezierControlPoints3D(sf1d_dir , NGeo , NGeo , BCSideID) )
!          ! Check if all corner points are equal is the "sf1d_dir" direction: Skip this side if true
!          IF((ALMOSTEQUALRELATIVE(x1,x2,1e-6).AND.&
!              ALMOSTEQUALRELATIVE(x1,x3,1e-6).AND.&
!              ALMOSTEQUALRELATIVE(x1,x4,1e-6))) CYCLE
!        END ASSOCIATE
!      END IF

!      IF(SideIndex(iSide).EQ.0)THEN
!        leave=.FALSE.
!        Vec1=0.

!        ! all points of bc side
!        DO q=firstBezierPoint,lastBezierPoint
!          DO p=firstBezierPoint,lastBezierPoint
!            NodeX(:) = BezierControlPoints3D(:,p,q,BCSideID)+Vec1
!            !all nodes of current side
!            DO s=firstBezierPoint,lastBezierPoint
!              DO r=firstBezierPoint,lastBezierPoint
!                dX=ABS(xNodes(1,r,s)-NodeX(1))
!                IF(dX.GT.r_sf) CYCLE
!                dY=ABS(xNodes(2,r,s)-NodeX(2))
!                IF(dY.GT.r_sf) CYCLE
!                dZ=ABS(xNodes(3,r,s)-NodeX(3))
!                IF(dZ.GT.r_sf) CYCLE
!                IF(SQRT(dX*dX+dY*dY+dZ*dZ).LE.r_sf)THEN
!                  IF(SideIndex(iSide).EQ.0)THEN
!                    BCElemSF=.TRUE.
!                    SideIndex(iSide)=1 ! mark with number .NE. 0
!                    leave=.TRUE.
!                    EXIT
!                  END IF
!                END IF
!              END DO ! r
!              IF(leave) EXIT
!            END DO ! s
!            IF(leave) EXIT
!          END DO ! p
!          IF(leave) EXIT
!        END DO ! q
!        IF(leave) EXIT
!      END IF ! SideIndex(iSide).EQ.0
!    END DO ! iSide=1,nTotalSides
!  END DO ! ilocSide=1,6

!  ! set true, only required for elements without an own bc side
!  IF(BCElemSF) IsLocalDepositionBCElem(iElem)=.TRUE.
!END DO ! iElem=1,nTotalElems




!END SUBROUTINE GetShapeFunctionBCElems_OLD


!USBROUTINE CalcElemAndSideNum()
!!===================================================================================================================================
!! calculate number of different elements types and side types for local and halo region
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
!USE MOD_Particle_Surfaces_Vars ,ONLY: SideType
!USE MOD_Particle_Mesh_Vars     ,ONLY: ElemCurved
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalSides,IsTracingBCElem,nTotalElems
!USE MOD_Particle_Mesh_Vars     ,ONLY: nPartSides
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalBCSides
!#if USE_MPI
!USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
!USE MOD_Particle_MPI_HALO      ,ONLY: WriteParticlePartitionInformation
!#endif /*USE_MPI*/
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                                  :: iElem
!INTEGER                                  :: iSide
!INTEGER                                  :: nBCElems,nBCelemsTot
!INTEGER                                  :: nPlanarRectangular, nPlanarNonRectangular,nPlanarCurved,nBilinear,nCurved
!INTEGER                                  :: nPlanarRectangularTot, nPlanarNonRectangularTot,nPlanarCurvedTot,nBilinearTot,nCurvedTot
!INTEGER                                  :: nLinearElems, nCurvedElems, nCurvedElemsTot
!#if USE_MPI
!INTEGER                                  :: nPlanarRectangularHalo, nPlanarNonRectangularHalo,nPlanarCurvedHalo, &
!                                            nBilinearHalo,nCurvedHalo,nCurvedElemsHalo,nLinearElemsHalo,nBCElemsHalo,nDummy
!#endif /*USE_MPI*/
!INTEGER                                  :: nLoop
!!===================================================================================================================================
!
!! zero counter for side and elem types
!nPlanarRectangular         = 0
!nPlanarNonRectangular      = 0
!nPlanarCurved              = 0
!nBilinear                  = 0
!nCurved                    = 0
!nBCElems                   = 0
!nCurvedElems               = 0
!nLinearElems               = 0
!#if USE_MPI
!nPlanarRectangularHalo     = 0
!nPlanarNonRectangularHalo  = 0
!nPlanarCurvedHalo          = 0
!nBilinearHalo              = 0
!nCurvedHalo                = 0
!nCurvedElemsHalo           = 0
!nLinearElemsHalo           = 0
!nBCElemsHalo               = 0
!#endif /*USE_MPI*/
!
!DO iElem=1,nTotalElems
!  ! count elements by type and in own and halo region
!  IF(iElem.LE.PP_nElems)THEN
!    IF(ElemCurved(iElem))THEN
!      nCurvedElems=nCurvedElems+1
!    ELSE
!      nLinearElems=nLinearElems+1
!    END IF
!    IF(IsTracingBCElem(iElem))THEN
!      nBCElems=nBCElems+1
!    END IF ! count only single
!#if USE_MPI
!  ELSE
!    IF(ElemCurved(iElem)) THEN
!      nCurvedElemsHalo=nCurvedElemsHalo+1
!    ELSE
!      nLinearElemsHalo=nLinearElemsHalo+1
!    END IF
!    IF(IsTracingBCElem(iElem))THEN
!      nBCElemsHalo=nBCElemsHalo+1
!    END IF ! count only single
!#endif /*USE_MPI*/
!  END IF
!END DO
!nLoop = nTotalSides
!IF (DoRefMapping) nLoop = nTotalBCSides
!DO iSide=1,nLoop
!  IF (iSide.LE.nPartSides) THEN
!    SELECT CASE(SideType(iSide))
!    CASE (PLANAR_RECT)
!      nPlanarRectangular=nPlanarRectangular+1
!    CASE (PLANAR_NONRECT)
!      nPlanarNonRectangular=nPlanarNonRectangular+1
!    CASE (BILINEAR)
!      nBilinear = nBilinear+1
!    CASE (PLANAR_CURVED)
!      nPlanarCurved = nPlanarCurved+1
!    CASE (CURVED)
!      nCurved = nCurved+1
!    END SELECT
!#if USE_MPI
!  ELSE IF (iSide.GT.nPartSides) THEN
!    SELECT CASE(SideType(iSide))
!    CASE (PLANAR_RECT)
!      nPlanarRectangularHalo=nPlanarRectangularHalo+1
!    CASE (PLANAR_NONRECT)
!      nPlanarNonRectangularHalo=nPlanarNonRectangularHalo+1
!    CASE (BILINEAR)
!      nBilinearHalo = nBilinearHalo+1
!    CASE (PLANAR_CURVED)
!      nPlanarCurvedHalo = nPlanarCurvedHalo+1
!    CASE (CURVED)
!      nCurvedHalo = nCurvedHalo+1
!    END SELECT
!#endif /*USE_MPI*/
!  END IF
!END DO
!
!#if USE_MPI
!IF(MPIRoot) THEN
!  CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nCurvedElems,nCurvedElemsTot,1,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
!  IF(DoRefMapping) CALL MPI_REDUCE(nBCElems,nBCElemsTot ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!ELSE ! no Root
!  CALL MPI_REDUCE(nPlanarRectangular     ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nPlanarNonRectangular  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nBilinear              ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nPlanarCurved          ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nCurved                ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!  CALL MPI_REDUCE(nCurvedElems           ,nDummy,1,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
!  IF(DoRefMapping) CALL MPI_REDUCE(nBCElems  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
!END IF
!#else
!nPlanarRectangularTot   =nPlanarRectangular
!nPlanarNonRectangularTot=nPlanarNonRectangular
!nBilinearTot            =nBilinear
!nPlanarCurvedTot        =nPlanarCurved
!nCurvedTot              =nCurved
!nCurvedElemsTot         =nCurvedElems
!IF(DorefMapping) nBCElemstot=nBCElems
!#endif /*USE_MPI*/
!
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar-rectangular     faces: ', nPlanarRectangulartot
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar-non-rectangular faces: ', nPlanarNonRectangulartot
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of bi-linear              faces: ', nBilineartot
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar-curved          faces: ', nPlanarCurvedtot
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of curved                 faces: ', nCurvedtot
!! and add number of curved elems
!IF(DoRefMapping)THEN
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of BC-adjoined            elems: ', nBCElemstot
!END IF
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of (bi-)linear            elems: ', nGlobalElems-nCurvedElemsTot
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of curved                 elems: ', nCurvedElemsTot
!SWRITE(UNIT_StdOut,'(132("-"))')
!#if USE_MPI
!CALL WriteParticlePartitionInformation(nPlanarRectangular+nPlanarNonRectangular,nBilinear,nCurved+nPlanarCurved,                    &
!                                       nPlanarRectangularHalo+nPlanarNonRectangularHalo,nBilinearHalo,nCurvedHalo+nPlanarCurvedHalo &
!                                      ,nBCElems,nLinearElems,nCurvedElems,nBCElemsHalo,nLinearElemsHalo,nCurvedElemsHalo)
!#endif
!
!END SUBROUTINE CalcElemAndSideNum


!SUBROUTINE ElemConnectivity()
!!===================================================================================================================================
!! computes the element connectivity between different elements, inclusive the halo region
!! and mortar interfaces
!! CAUTION: the assumption is, that one element is only linked once or twice with another element
!!          one link: normal inner connection or periodic connection
!!          two links: one normal connection PLUS one periodic connection
!!          more than 2 links: funny.
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToElemGlob, PartElemToElemAndSide,nTotalElems,PartElemToSide,PartBCSideList &
!                                 ,SidePeriodicType,ElemToGlobalElemID
!USE MOD_Mesh_Vars              ,ONLY: OffSetElem,BC,BoundaryType,MortarType
!USE MOD_Particle_Surfaces_Vars ,ONLY: SideNormVec
!USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!#if USE_MPI
!USE MOD_MPI_Vars               ,ONLY: OffSetElemMPI
!USE MOD_Particle_MPI_Vars      ,ONLY: PartHaloElemToProc
!#endif /*USE_MPI*/
!USE MOD_Mesh_vars
!!----------------------------------------------------------------------------------------------------------------------------------!
!! insert modules here
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                       :: iElem,ilocSide,iMortar,ilocSide2,iMortar2,NbElemID,ElemID,BCID,SideID,BCSideID
!INTEGER(KIND=8)               :: GlobalElemID
!LOGICAL                       :: found
!REAL                          :: Vec1(1:3)
!#if USE_MPI
!INTEGER                       :: iHaloElem,ProcID
!INTEGER(KIND=8)               :: HaloGlobalElemID
!#endif /*USE_MPI*/
!!===================================================================================================================================
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)')' BUILD MESH-CONNECTIVITY ... '

!SDEALLOCATE(PartElemToElemAndSide)
!ALLOCATE(PartElemToElemAndSide(1:8,1:6,1:nTotalElems))
!                      ! [1]1:4 - MortarNeighborElemID
!                      ! [1]5:8 -       Neighbor locSideID
!                      ! [2]1:6 - locSideID
!                      ! [3]    - nTotalElems
!                      ! if the connections points to an element which is not in MY region (MY elems + halo elems)
!                      ! then this connection points to -1
!ALLOCATE(ElemToGlobalElemID(1:nTotalElems))
!! nullify
!PartElemToElemAndSide=-1
!ElemToGlobalElemID=-1

!! now, map the PartElemToElemGlob to local elemids
!! loop over all Elems and map the neighbor element to local coordinates
!DO iElem=1,nTotalElems
!  IF(iElem.LE.nElems)THEN
!    ElemToGlobalElemID(iElem)=offSetElem+iElem
!#if USE_MPI
!  ELSE
!    ProcID=PartHaloElemToProc(NATIVE_PROC_ID,iElem)
!    ElemToGlobalElemID(iElem)=offSetElemMPI(ProcID) + PartHaloElemToProc(NATIVE_ELEM_ID,iElem)
!#endif /*USE_MPI*/
!  END IF
!  DO ilocSide=1,6
!    DO iMortar=1,4
!      GlobalElemID=PartElemToElemGlob(iMortar,ilocSide,iElem)
!      IF(GlobalElemID.LE.0) CYCLE
!      ! check if the element is in MY range of elements
!      IF((GlobalElemID.GE.OffSetElem+1).AND.(GlobalElemID.LE.(OffSetElem+PP_nElems)))THEN
!        PartElemToElemAndSide(iMortar,ilocSide,iElem)=INT(GlobalElemID-OffSetElem,4)
!        CYCLE
!      END IF
!#if USE_MPI
!      ! neighbor element not found, hence, it can be a halo element
!      DO iHaloElem=PP_nElems+1,nTotalElems
!        ProcID=PartHaloElemToProc(NATIVE_PROC_ID,iHaloElem)
!        HaloGlobalElemID=offSetElemMPI(ProcID) + PartHaloElemToProc(NATIVE_ELEM_ID,iHaloElem)
!        CHECKSAFEINT(HaloGlobalElemID,4)
!        IF(HaloGlobalElemID.EQ.GlobalElemID)THEN
!          PartElemToElemAndSide(iMortar,ilocSide,iElem)=iHaloElem
!          EXIT
!        END IF
!      END DO ! iHaloElem=1,nTotalElems
!#endif /*USE_MPI*/
!    END DO ! iMortar=1,4
!  END DO ! ilocSide=1,6
!END DO ! iElem=1,PP_nElems

!! which local side of neighbor element is connected to MY element
!DO iElem=1,nTotalElems
!  DO ilocSide=1,6
!    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    ! check for ref-mapping or tracing
!    IF(DoRefMapping)THEN
!      IF(SideID.GT.0)THEN
!        BCSideID=PartBCSideList(SideID)
!      ELSE
!        BCSideID=-1
!      END IF
!    ELSE
!      BCSideID=SideID
!    END IF
!    ! disable BCSideID, if it is NOT a periodic side
!    IF(BCSideID.GT.0)THEN ! only BC faces
!      IF(SidePeriodicType(SideID).NE.0)THEN ! only periodic sides
!        Vec1=SideNormVec(1:3,BCSideID)
!        IF(ALMOSTZERO(DOT_PRODUCT(Vec1,Vec1))) CALL abort(&
!__STAMP__&
!        , ' Error in ElemConnectivity. No SideNormVec!',iElem,REAL(ilocSide))
!      ELSE ! disable non-periodic  sides
!        Vec1=0.
!        BCSideID=-1
!      END IF
!    END IF
!    IF(BCSideID.GT.0)THEN ! periodic sides
!      DO iMortar=1,4
!        NBElemID=PartElemToElemAndSide(iMortar,ilocSide,iElem)
!        IF(NBElemID.EQ.-1) CYCLE
!        found=.FALSE.
!        ! loop  over all local sides of neighbor element to find the right face
!        DO ilocSide2=1,6
!          DO iMortar2=1,4
!            ElemID=PartElemToElemAndSide(iMortar2,ilocSide2,NBElemID)
!            IF(ElemID.LE.0) CYCLE
!            IF(ElemID.EQ.iElem) THEN
!              ! check if periodic side
!              SideID=PartElemToSide(E2S_SIDE_ID,ilocSide2,NBElemID)
!              ! check for ref-mapping or tracing
!              IF(DoRefMapping)THEN
!                IF(SideID.GT.0)THEN
!                  BCSideID=PartBCSideList(SideID)
!                ELSE
!                  BCSideID=-1
!                END IF
!              ELSE
!                BCSideID=SideID
!              END IF
!              IF(BCSideID.GT.0)THEN ! only BC faces
!                IF(SidePeriodicType(SideID).NE.0)THEN ! only periodic sides
!                  IF(ALMOSTEQUAL(ABS(DOT_PRODUCT(Vec1,SideNormVec(1:3,BCSideID))),1.0))THEN
!                    ! finally, found matching local sides
!                    PartElemToElemAndSide(iMortar+4,ilocSide,iElem)=ilocSide2
!                    Found=.TRUE.
!                    EXIT
!                  ELSE
!                    CYCLE
!                  END IF
!                ELSE ! disable non-periodic  sides
!                  CYCLE
!                END IF
!              ELSE
!                CYCLE
!              END IF
!            END IF
!          END DO ! iMortar=1,4
!          IF(Found) EXIT
!        END DO ! ilocSide=1,6
!      END DO ! iMortar=1,4
!    ELSE ! non-periodic sides
!      DO iMortar=1,4
!        NBElemID=PartElemToElemAndSide(iMortar,ilocSide,iElem)
!        IF(NBElemID.EQ.-1) CYCLE
!        found=.FALSE.
!        ! loop  over all local sides of neighbor element to find the right face
!        DO ilocSide2=1,6
!          DO iMortar2=1,4
!            ElemID=PartElemToElemAndSide(iMortar2,ilocSide2,NBElemID)
!            IF(ElemID.LE.0) CYCLE
!            IF(ElemID.EQ.iElem) THEN
!              ! check if periodic side
!              SideID=PartElemToSide(E2S_SIDE_ID,ilocSide2,NBElemID)
!              ! check for ref-mapping or tracing
!              IF(DoRefMapping)THEN
!                IF(SideID.GT.0)THEN
!                  BCSideID=PartBCSideList(SideID)
!                ELSE
!                  BCSideID=-1
!                END IF
!              ELSE
!                BCSideID=SideID
!              END IF
!              IF(BCSideID.GT.0)THEN ! BC face?
!                IF(SidePeriodicType(SideID).NE.0)THEN ! only non-periodic sides
!                  CYCLE
!                ELSE ! enable non-periodic  sides
!                  ! finally, found matching local sides
!                  PartElemToElemAndSide(iMortar+4,ilocSide,iElem)=ilocSide2
!                  Found=.TRUE.
!                  EXIT
!                END IF
!              ELSE
!                ! finally, found matching local sides
!                PartElemToElemAndSide(iMortar+4,ilocSide,iElem)=ilocSide2
!                Found=.TRUE.
!                EXIT
!              END IF
!            END IF
!          END DO ! iMortar=1,4
!          IF(Found) EXIT
!        END DO ! ilocSide=1,6
!      END DO ! iMortar=1,4
!    END IF ! periodic sides
!  END DO ! ilocSide=1,6
!END DO ! iElem=1,PP_nElems

!! sanity check
!DO iElem=1,nTotalElems
!  DO ilocSide=1,6
!    DO iMortar=1,4
!      IF((PartElemToElemAndSide(iMortar,ilocSide,iElem).GT.0).AND.(PartElemToElemAndSide(iMortar+4,ilocSide,iElem).EQ.-1))THEN
!        IPWRITE(UNIT_StdOut,*) ' iElem:     ', iElem
!        IPWRITE(UNIT_StdOut,*) ' ilocSide:  ', ilocSide
!        IPWRITE(UNIT_StdOut,*) ' NBElem-ID: ', PartElemToElemAndSide(iMortar,ilocSide,iElem)
!        IPWRITE(UNIT_StdOut,*) ' NBElem-ID: ', PartElemToElemAndSide(iMortar+4,ilocSide,iElem)
!        CALL abort(&
!__STAMP__&
!        , ' Error in ElemConnectivity. Found no neighbor locSideID. iElem,ilocSide',iElem,REAL(ilocSide))
!      END IF
!      IF((PartElemToElemAndSide(iMortar,ilocSide,iElem).EQ.-1).AND.(PartElemToElemAndSide(iMortar+4,ilocSide,iElem).GT.-1))THEN
!        IPWRITE(UNIT_StdOut,*) ' iElem:     ', iElem
!        IPWRITE(UNIT_StdOut,*) ' ilocSide:  ', ilocSide
!        IPWRITE(UNIT_StdOut,*) ' NBElem-ID: ', PartElemToElemAndSide(iMortar,ilocSide,iElem)
!        CALL abort(&
!__STAMP__&
!        , ' Error in ElemConnectivity. Found no neighbor ElemID. iElem,ilocSide',iElem,REAL(ilocSide))
!      END IF
!    END DO ! iMortar=1,4
!  END DO ! ilocSide=1,6
!END DO

!IF(nGlobalMortarSides.EQ.0) THEN
!  ! check is working on CONFORM mesh!
!  DO iElem=1,nTotalElems
!    DO ilocSide=1,6
!      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!      IF(DoRefMapping)THEN
!        IF(SideID.LT.1) CYCLE
!      ELSE
!        IF(SideID.LE.0) CALL abort(&
!__STAMP__&
!         , ' Error in PartElemToSide! No SideID for side!. iElem,ilocSide',iElem,REAL(ilocSide))
!      END IF
!      IF(MortarType(1,SideID).NE.0) CYCLE
!      BCID=BC(SideID)
!      IF(BCID.NE.0)THEN
!        IF(BoundaryType(BCID,BC_TYPE).GT.1) CYCLE
!      END IF
!      IF(PartElemToElemAndSide(1,ilocSide,iElem).LT.1)THEN
!         CALL abort(&
!__STAMP__&
!        , ' Error in ElemConnectivity. Found no neighbor ElemID. iElem,ilocSide',iElem,REAL(ilocSide))
!        END IF
!    END DO ! ilocSide=1,6
!  END DO
!END IF

!#if USE_MPI
!CALL MPI_BARRIER(MPI_COMM_WORLD,iERROR)
!#endif
!SWRITE(UNIT_stdOut,'(A)')' BUILD MESH-CONNECTIVITY SUCCESSFUL '
!SWRITE(UNIT_StdOut,'(132("-"))')

!END SUBROUTINE ElemConnectivity


!SUBROUTINE NodeNeighbourhood()
!!===================================================================================================================================
!!> Subroutine for initialization of neighbourhood with nodes using GEO container
!!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Globals
!USE MOD_Mesh_Vars          ,ONLY: nNodes
!USE MOD_Particle_Mesh_Vars ,ONLY: GEO, PartElemToElemAndSide
!#ifdef CODE_ANALYZE
!#if USE_MPI
!USE MOD_Mesh_Vars          ,ONLY: offsetElem
!USE MOD_MPI_Vars           ,ONLY: offsetElemMPI
!USE MOD_MPI_Shared_Vars
!USE MOD_Particle_MPI_Vars  ,ONLY: PartHaloElemToProc
!#endif /*USE_MPI*/
!#endif /*CODE_ANALYZE*/
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!TYPE tNodeToElem
!  INTEGER :: ElemID(50)
!END TYPE
!TYPE(tNodeToElem)      :: TempNodeToElem(1:nNodes)
!INTEGER                :: TempElemsOnNode(1:nNodes)
!INTEGER                :: Element, iLocSide, k, l,iMortar
!LOGICAL                :: ElemExists
!INTEGER                :: iElem, jNode
!INTEGER                :: iNode
!INTEGER                :: TempHaloElems(1:500)
!INTEGER                :: TempHaloNumElems
!!#if USE_MPI
!!LOGICAL                :: HaloNeighNode(1:nNodes)
!!#endif /*USE_MPI*/
!LOGICAL                :: ElemDone
!REAL                   :: MPINodeCoord(3), ElemCoord(3)
!!===================================================================================================================================
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)')' BUILD NODE-NEIGHBOURHOOD ... '
!
!!#if USE_MPI
!!! set nodes of sides with halo element connected to it as HaloNeighNodes
!!GEO%HaloNeighNode(:) = .FALSE.
!!DO iElem=1,nElems
!!  DO iLocSide = 1,6
!!    IF (PartElemToElemAndSide(1,iLocSide,iElem).GT.PP_nElems) THEN
!!      DO iNode = 1,4
!!        GEO%HaloNeighNode(GEO%ElemSideNodeID(iNode,iLocSide,iElem)) = .TRUE.
!!      END DO
!!    END IF
!!  END DO
!!END DO
!!#ifdef CODE_ANALYZE
!!DO iNode=1,nNodes
!!  IF (GEO%HaloNeighNode(iNode)) THEN
!!    print*,'Rank: ',MyRank,'---- local Node: ',iNode,' is halo node'
!!  END IF
!!END DO
!!#endif /*CODE_ANALYZE*/
!!#endif /*USE_MPI*/
!
!ALLOCATE(GEO%NumNeighborElems(1:PP_nElems))
!ALLOCATE(GEO%ElemToNeighElems(1:PP_nElems))
!GEO%NumNeighborElems(:)=0
!TempElemsOnNode(:)=0
!DO iNode=1,nNodes
!  TempNodeToElem(iNode)%ElemID=-1
!END DO
!
!! find all real neighbour elements for elements with halo neighbours
!! recursively checking connected halo area
!DO iElem =1, PP_nElems
!  TempHaloElems(1:500) = 0
!  TempHaloNumElems = 0
!  ! now check every side for neighbours, add valid neighbour to corresponding array and proceed recursively until neighbourhood
!  ! is finished
!  DO iLocSide = 1, 6
!    DO iMortar=1,4
!      ElemExists = .FALSE.
!      Element = PartElemToElemAndSide(iMortar,iLocSide,iElem)
!      IF (Element.GT.0) THEN !side has neighbour element
!        DO l=1, TempHaloNumElems
!          IF(Element.EQ.TempHaloElems(l)) THEN
!            ElemExists=.TRUE.
!            EXIT
!          END IF
!        END DO
!        IF (.NOT.ElemExists) THEN
!         TempHaloNumElems = TempHaloNumElems + 1
!          TempHaloElems(TempHaloNumElems) = Element
!        END IF
!        CALL RecurseCheckNeighElems(iElem,Element,TempHaloNumElems,TempHaloElems)
!      END IF
!    END DO
!  END DO
!  IF (TempHaloNumElems.LE.0) CALL abort(&
!__STAMP__&
!,'ERROR in FindNeighbourElems! no neighbour elements found for Element',iElem)
!  ! write local variables into global array
!  GEO%NumNeighborElems(iElem) = TempHaloNumElems
!  ALLOCATE(GEO%ElemToNeighElems(iElem)%ElemID(1:GEO%NumNeighborElems(iElem)))
!  GEO%ElemToNeighElems(iElem)%ElemID(1:GEO%NumNeighborElems(iElem)) = TempHaloElems(1:TempHaloNumElems)
!  ! add neighbour elements to respective nodes of iElem if nodes of neighbour are the same to those of iElem
!  DO l=1,GEO%NumNeighborElems(iElem)
!    Element = GEO%ElemToNeighElems(iElem)%ElemID(l)
!    DO iNode = 1, 8
!      ElemExists = .FALSE.
!!      DO k = 1,TempElemsOnNode(ElemNodeID_Shared(iNode,iElem))
!      DO k = 1,TempElemsOnNode(ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+iNode)
!        IF (TempNodeToElem(ElemNodeID_Shared(iNode,iElem))%ElemID(k).EQ.Element) THEN
!          ElemExists = .TRUE.
!          EXIT
!        END IF
!      END DO
!      IF (ElemExists) CYCLE
!      ElemCoord(1:3) = NodeCoords_Shared(1:3,ElemNodeID_Shared(iNode,iElem))
!      DO jNode = 1, 8
!        ElemDone = .FALSE.
!        MPINodeCoord(1:3) = NodeCoords_Shared(1:3,ElemNodeID_Shared(jNode,Element))
!        IF(ALMOSTEQUAL(MPINodeCoord(1),ElemCoord(1)).AND.ALMOSTEQUAL(MPINodeCoord(2),ElemCoord(2)) &
!            .AND.ALMOSTEQUAL(MPINodeCoord(3),ElemCoord(3))) THEN
!          TempElemsOnNode(ElemNodeID_Shared(iNode,iElem)) = TempElemsOnNode(ElemNodeID_Shared(iNode,iElem)) + 1
!          TempNodeToElem(ElemNodeID_Shared(iNode,iElem))%ElemID(TempElemsOnNode(ElemNodeID_Shared(iNode,iElem))) = Element
!          ElemDone = .TRUE.
!        END IF
!        IF (ElemDone) EXIT
!      END DO ! jNode=1,8
!      IF (ElemDone) CYCLE
!    END DO ! iNode=1,8
!  END DO ! l=1,NumNeihborElems
!END DO ! iElem=1,PP_nElems
!
!! check if current element already added to every node of the element and add to missing (elements for corner nodes not added yet)
!! this part needs its own loop over all elements
!DO iElem=1,PP_nElems
!  DO iNode = 1,8
!    ElemExists = .FALSE.
!    DO l = 1,TempElemsOnNode(ElemNodeID_Shared(iNode,iElem))
!      IF (TempNodeToElem(ElemNodeID_Shared(iNode,iElem))%ElemID(l).EQ.iElem) THEN
!        ElemExists = .TRUE.
!        EXIT
!      END IF
!    END DO
!    IF (.NOT.ElemExists) THEN
!      TempElemsOnNode(ElemNodeID_Shared(iNode,iElem)) = TempElemsOnNode(ElemNodeID_Shared(iNode,iElem)) + 1
!      TempNodeToElem(ElemNodeID_Shared(iNode,iElem))%ElemID(TempElemsOnNode(ElemNodeID_Shared(iNode,iElem))) = iElem
!    END IF
!  END DO
!END DO
!
!! write number of elements for corresponding proc global nodes into GEO structure
!ALLOCATE(GEO%ElemsOnNode(1:nNodes))
!ALLOCATE(GEO%NodeToElem(1:nNodes))
!DO iNode=1,nNodes
!  GEO%ElemsOnNode(iNode) = TempElemsOnNode(iNode)
!  ALLOCATE(GEO%NodeToElem(iNode)%ElemID(1:GEO%ElemsOnNode(iNode)))
!  GEO%NodeToElem(iNode)%ElemID(1:GEO%ElemsOnNode(iNode)) = TempNodeToElem(iNode)%ElemID(1:TempElemsOnNode(iNode))
!END DO
!
!! fill array of neighbour nodes to proc global nodes
!ALLOCATE(GEO%NeighNodesOnNode(1:nNodes))
!ALLOCATE(GEO%NodeToNeighNode(1:nNodes))
!DO iNode=1,nNodes
!  TempHaloElems(:) = 0
!  TempHaloNumElems = 0
!  DO iElem=1,GEO%ElemsOnNode(iNode)
!    DO jNode=1,8
!      IF (ElemNodeID_Shared(jNode,GEO%NodeToElem(iNode)%ElemID(iElem)).EQ.iNode) CYCLE
!      ElemExists=.false.
!      DO l=1, TempHaloNumElems
!        IF(ElemNodeID_Shared(jNode,GEO%NodeToElem(iNode)%ElemID(iElem)).EQ.TempHaloElems(l)) THEN
!          ElemExists=.true.
!          CYCLE
!        END IF
!      END DO
!      IF(.NOT.ElemExists) THEN
!        TempHaloNumElems = TempHaloNumElems + 1
!        TempHaloElems(TempHaloNumElems) = ElemNodeID_Shared(jNode,GEO%NodeToElem(iNode)%ElemID(iElem))
!      END IF
!    END DO
!  END DO
!  ! write local variables into global array
!  ALLOCATE(GEO%NodeToNeighNode(iNode)%ElemID(1:TempHaloNumElems))
!  GEO%NeighNodesOnNode(iNode)=TempHaloNumElems
!  GEO%NodeToNeighNode(iNode)%ElemID(1:TempHaloNumElems) = TempHaloElems(1:TempHaloNumElems)
!END DO
!
!#ifdef CODE_ANALYZE
!! write some code analyze output of connectivity
!DO iElem=1,PP_nElems
!#if USE_MPI
!  IPWRITE(UNIT_StdOut,*) '------ Element: ',iElem+offsetElem,' has ',GEO%NumNeighborElems(iElem),' Neighbours'
!  IPWRITE(UNIT_StdOut,*) '------ Neighbours are:'
!  DO l=1,GEO%NumNeighborElems(iElem)
!    IF (GEO%ElemToNeighElems(iElem)%ElemID(l).GT.PP_nElems) THEN
!      IPWRITE(UNIT_StdOut,*) offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,GEO%ElemToNeighElems(iElem)%ElemID(l))) &
!          + PartHaloElemToProc(NATIVE_ELEM_ID,GEO%ElemToNeighElems(iElem)%ElemID(l))
!    ELSE
!      IPWRITE(UNIT_StdOut,*) GEO%ElemToNeighElems(iElem)%ElemID(l) + offsetElem
!    END IF
!  END DO
!#else
!  IPWRITE(UNIT_StdOut,*) '------ Element: ',iElem,' has ',GEO%NumNeighborElems(iElem),' Neighbours'
!  IPWRITE(UNIT_StdOut,*) '------ Neighbours are:',GEO%ElemToNeighElems(iElem)%ElemID(:)
!#endif /*USE_MPI*/
!END DO
!
!DO iNode=1,nNodes
!  IPWRITE(UNIT_StdOut,*) '------ Node: ',iNode,' has: ',GEO%ElemsOnNode(iNode),' Elements'
!  IPWRITE(UNIT_StdOut,*) '------ Node: ',iNode,' has: ',GEO%ElemsOnNode(iNode),' Elements'
!END DO
!#endif /*CODE_ANALYZE*/
!
!SWRITE(UNIT_stdOut,'(A)')' BUILD NODE-NEIGHBOURHOOD SUCCESSFUL '
!SWRITE(UNIT_StdOut,'(132("-"))')
!END SUBROUTINE NodeNeighbourhood
!
!
!RECURSIVE SUBROUTINE RecurseCheckNeighElems(StartElem,HaloElem,TempHaloNumElems,TempHaloElems)
!!===================================================================================================================================
!!> Subroutine for recursively checking halo neighbourhood for connectivity to current elem
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_PreProc
!USE MOD_Particle_Mesh_Vars ,ONLY: GEO, PartElemToElemAndSide
!#if USE_MPI
!USE MOD_MPI_Shared_Vars
!#endif
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN)     :: StartElem
!INTEGER,INTENT(INOUT)  :: HaloElem,TempHaloElems(1:500), TempHaloNumElems
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                :: iNode, jNode
!INTEGER                :: iLocSide, l, iMortar
!INTEGER                :: currentElem
!LOGICAL                :: ElemExists, ElemDone
!REAL                   :: MPINodeCoord(3), ElemCoord(3)
!!===================================================================================================================================
!DO iLocSide = 1,6
!  DO iMortar=1,4
!    ElemExists = .FALSE.
!    currentElem = PartElemToElemAndSide(iMortar,iLocSide,HaloElem)
!    IF (currentElem.GT.0 .AND. currentElem.NE.StartElem) THEN
!    !IF (currentElem.GT.PP_nElems) THEN
!      DO l=1, TempHaloNumElems
!        IF(currentElem.EQ.TempHaloElems(l)) THEN
!          ElemExists=.TRUE.
!          EXIT
!        END IF
!      END DO
!      IF (.NOT.ElemExists) THEN
!        ElemDone = .FALSE.
!        DO iNode = 1, 8
!          DO jNode = 1, 8
!            MPINodeCoord(1:3) = NodeCoords_Shared(1:3,ElemNodeID_Shared(jNode,currentElem))
!            ElemCoord(1:3) = NodeCoords_Shared(1:3,ElemNodeID_Shared(iNode,StartElem))
!            IF(ALMOSTEQUAL(MPINodeCoord(1),ElemCoord(1)).AND.ALMOSTEQUAL(MPINodeCoord(2),ElemCoord(2)) &
!                .AND.ALMOSTEQUAL(MPINodeCoord(3),ElemCoord(3))) THEN
!              TempHaloNumElems = TempHaloNumElems + 1
!              TempHaloElems(TempHaloNumElems) = currentElem
!              ElemDone = .TRUE.
!              CALL RecurseCheckNeighElems(StartElem,currentElem,TempHaloNumElems,TempHaloElems)
!            END IF
!            IF (ElemDone) EXIT
!          END DO
!          IF (ElemDone) EXIT
!        END DO
!      END IF
!    END IF
!  END DO
!END DO
!
!END SUBROUTINE RecurseCheckNeighElems
!
!#if USE_MPI
!SUBROUTINE BuildLocNodeToHaloNodeComm()
!!===================================================================================================================================
!!> build all missing stuff for node communication, like
!!> MPI-neighbor list
!!> PartHaloNodeToProc
!!> send and recv list mapping
!!> Only receiving process knows to which local node the information is send
!!> The sending process does not know the final nodeID
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Mesh_Vars          ,ONLY: nNodes
!USE MOD_Particle_MPI_Vars  ,ONLY: PartMPI
!USE MOD_Particle_MPI_Vars  ,ONLY: PartHaloNodeToProc
!USE MOD_Particle_MPI_Vars  ,ONLY: NodeSendBuf, NodeRecvBuf, NodeExchange
!USE MOD_Particle_Mesh_Vars ,ONLY: nTotalNodes, GEO
!#if USE_MPI
!USE MOD_MPI_Shared_Vars
!#endif
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER             :: NodeIndexToSend(1:nNodes,0:PartMPI%nProcs-1)
!INTEGER             :: nDOF,ALLOCSTAT
!INTEGER             :: iProc, iNode,NodeID,iElem,iSendNode,iRecvNode,iPos,jNode
!INTEGER,ALLOCATABLE :: recv_status_list(:,:)
!INTEGER             :: NativeNodeID, iMPINeighbor
!REAL                :: MPINodeCoord(3), ElemCoord(3)
!!===================================================================================================================================
!
!! get list of mpi Node neighbors in halo-region (first mappings)
!  ! caution:
!  !  mapping1 is only done for halo-nodes with direct connection to local nodes (equal nodes to local)
!ALLOCATE(PartMPI%IsMPINodeNeighbor(0:PartMPI%nProcs-1))
!PartMPI%IsMPINodeNeighbor(:) = .FALSE.
!PartMPI%nMPINodeNeighbors = 0
!
!NodeIndexToSend(:,:) = -1
!IF(nTotalNodes.GT.nNodes)THEN
!  ! get all MPI-neighbors to communicate with (only direct halo border)
!  ! list needed because not all halo nodes are considered for first mapping
!  DO iProc=0,PartMPI%nProcs-1
!    IF(iProc.EQ.PartMPI%MyRank) CYCLE
!    DO iNode=1,nNodes
!      DO iElem=1,GEO%ElemsOnNode(iNode)
!        IF (GEO%NodeToElem(iNode)%ElemID(iElem).LE.PP_nElems) CYCLE
!        DO jNode=1,8
!          IF ( iProc.NE.PartHaloNodeToProc(NATIVE_PROC_ID,ElemNodeID_Shared(jNode,GEO%NodeToElem(iNode)%ElemID(iElem))) ) CYCLE
!          MPINodeCoord(1:3) = NodeCoords_Shared(1:3,ElemNodeID_Shared(jNode,GEO%NodeToElem(iNode)%ElemID(iElem)))
!          ElemCoord(1:3) = NodeCoords_Shared(1:3,iNode)
!          IF(ALMOSTEQUAL(MPINodeCoord(1),ElemCoord(1)).AND.ALMOSTEQUAL(MPINodeCoord(2),ElemCoord(2)) &
!              .AND.ALMOSTEQUAL(MPINodeCoord(3),ElemCoord(3))) THEN
!            NodeIndexToSend(iNode,iProc) = ElemNodeID_Shared(jNode,GEO%NodeToElem(iNode)%ElemID(iElem)) ! index of halo node on local proc
!            IF (.NOT.PartMPI%IsMPINodeNeighbor(iProc)) THEN
!              PartMPI%IsMPINodeNeighbor(iProc) = .TRUE.
!              PartMPI%nMPINodeNeighbors = PartMPI%nMPINodeNeighbors + 1
!            END IF
!            EXIT
!          END IF
!        END DO
!        IF (NodeIndexToSend(iNode,iProc).GT.0) EXIT
!      END DO ! iElem=1,GEO%ElemsOnNode(iNode)
!    END DO ! iNode=1,nNodes
!  END DO ! iProc = 0, PartMPI%nProcs-1
!END IF
!
!! fill list with neighbor proc id and add local neighbor id to PartHaloNodeToProc
!ALLOCATE( PartMPI%MPINodeNeighbor(PartMPI%nMPINodeNeighbors))
!iMPINeighbor=0
!DO iProc=0,PartMPI%nProcs-1
!  ! Check if iProc is my node neighbour
!  IF(PartMPI%IsMPINodeNeighbor(iProc))THEN
!    iMPINeighbor=iMPINeighbor+1
!    ! Mapping of node neighbour proc to global proc id (PartMPI%COMM)
!    PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID=iProc
!    ! Loop all halo nodes
!    DO iNode=nNodes+1,nTotalNodes
!      IF(iProc.EQ.PartHaloNodeToProc(NATIVE_PROC_ID,iNode)) PartHaloNodeToProc(LOCAL_PROC_ID,iNode)=iMPINeighbor
!    END DO ! iNode
!  END IF
!END DO
!
!! array how many nodes have to be communicated
!ALLOCATE(NodeExchange%nNodesSend(1:PartMPI%nMPINodeNeighbors) &
!        ,NodeExchange%nNodesRecv(1:PartMPI%nMPINodeNeighbors) &
!        ,NodeExchange%SendRequest(PartMPI%nMPINodeNeighbors)  &
!        ,NodeExchange%RecvRequest(PartMPI%nMPINodeNeighbors)  )
!NodeExchange%nNodesSend(:) = 0
!NodeExchange%nNodesRecv(:) = 0
!
!! count number of nodes to send to each proc
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  DO iNode=1,nNodes
!    IF (NodeIndexToSend(iNode,PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID).GT.nNodes) THEN
!      NodeExchange%nNodesSend(iMPINeighbor) = NodeExchange%nNodesSend(iMPINeighbor) + 1
!    END IF
!  END DO
!END DO
!
!! open envelope receiving number of send nodes
!ALLOCATE(RECV_STATUS_LIST(1:MPI_STATUS_SIZE,1:PartMPI%nMPINodeNeighbors))
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  CALL MPI_IRECV( NodeExchange%nNodesRecv(iMPINeighbor)            &
!                , 1                                                &
!                , MPI_INTEGER                                      &
!                , PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID &
!                , 1313                                             &
!                , PartMPI%COMM                                     &
!                , NodeExchange%RecvRequest(iMPINeighbor)           &
!                , IERROR )
!END DO ! iMPINeighbor
!
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  CALL MPI_ISEND( NodeExchange%nNodesSend(iMPINeighbor)            &
!                , 1                                                &
!                , MPI_INTEGER                                      &
!                , PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID &
!                , 1313                                             &
!                , PartMPI%COMM                                     &
!                , NodeExchange%SendRequest(iMPINeighbor)           &
!                , IERROR )
!END DO ! iMPINeighbor
!
!
!! 4) Finish Received number of nodes
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  CALL MPI_WAIT(NodeExchange%SendRequest(iMPINeighbor),MPIStatus,IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!,' MPI Communication error', IERROR)
!  CALL MPI_WAIT(NodeExchange%RecvRequest(iMPINeighbor),recv_status_list(:,iMPINeighbor),IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!END DO ! iMPINeighbor
!
!! allocate send and receive buffer for communicating send node mapping
!ALLOCATE(NodeSendBuf(PartMPI%nMPINodeNeighbors))
!ALLOCATE(NodeRecvBuf(PartMPI%nMPINodeNeighbors))
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  IF(NodeExchange%nNodesSend(iMPINeighbor).GT.0)THEN
!    ALLOCATE(NodeSendBuf(iMPINeighbor)%content(NodeExchange%nNodesSend(iMPINeighbor)),STAT=ALLOCSTAT)
!    NodeSendBuf(iMPINeighbor)%content(:)=0.
!  END IF
!  IF(NodeExchange%nNodesRecv(iMPINeighbor).GT.0)THEN
!    ALLOCATE(NodeRecvBuf(iMPINeighbor)%content(NodeExchange%nNodesRecv(iMPINeighbor)),STAT=ALLOCSTAT)
!    NodeRecvBuf(iMPINeighbor)%content(:)=0.
!  END IF
!END DO ! iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!
!! open receive buffer
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  IF(NodeExchange%nNodesRecv(iMPINeighbor).EQ.0) CYCLE
!  CALL MPI_IRECV( NodeRecvBuf(iMPINeighbor)%content                &
!                , NodeExchange%nNodesRecv(iMPINeighbor)            &
!                , MPI_DOUBLE_PRECISION                      &
!                , PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID &
!                , 1414                                      &
!                , PartMPI%COMM                              &
!                , NodeExchange%RecvRequest(iMPINeighbor)           &
!                , IERROR )
!END DO ! iMPINeighbor
!
!! build message
!! after this message, the receiving process knows to which of his nodes it receives and the sending process will know which nodes to
!! send
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  IF(NodeExchange%nNodesSend(iMPINeighbor).EQ.0) CYCLE
!  ALLOCATE(PartMPI%MPINodeNeighbor(iMPINeighbor)%SendList(NodeExchange%nNodesSend(iMPINeighbor)))
!  PartMPI%MPINodeNeighbor(iMPINeighbor)%SendList(:) = 0
!  iSendNode=0
!  iPos=1
!  DO iNode=1,nNodes
!    IF (NodeIndexToSend(iNode,PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID).GT.nNodes) THEN
!      iSendNode=iSendNode+1
!      PartMPI%MPINodeNeighbor(iMPINeighbor)%SendList(iSendNode)=iNode
!      NodeID=PartHaloNodeToProc(NATIVE_ELEM_ID,NodeIndexToSend(iNode,PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID))
!      NodeSendBuf(iMPINeighbor)%content(iPos)=REAL(NodeID)
!      iPos=iPos+1
!    END IF
!  END DO ! iNode=1,nNodes
!  IF(iSendNode.NE.NodeExchange%nNodesSend(iMPINeighbor)) CALL abort(&
!__STAMP__&
!          ,' Message for node-exchange in init too short!',iMPINeighbor)
!  IF(ANY(NodeSendBuf(iMPINeighbor)%content.LE.0))THEN
!    IPWRITE(UNIT_stdOut,*) ' nSendNodes', NodeExchange%nNodesSend(iMPINeighbor), ' to Proc ', iMPINeighbor
!    CALL abort(&
!__STAMP__&
!          ,' Sent Native-NodeID is < zero!')
!  END IF
!END DO
!
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  IF(NodeExchange%nNodesSend(iMPINeighbor).EQ.0) CYCLE
!  CALL MPI_ISEND( NodeSendBuf(iMPINeighbor)%content                &
!                , NodeExchange%nNodesSend(iMPINeighbor)            &
!                , MPI_DOUBLE_PRECISION                      &
!                , PartMPI%MPINodeNeighbor(iMPINeighbor)%COMMProcID &
!                , 1414                                      &
!                , PartMPI%COMM                              &
!                , NodeExchange%SendRequest(iMPINeighbor)           &
!                , IERROR )
!END DO ! iMPINeighbor
!
!! 4) Finish Received indexing of received nodes
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  IF(NodeExchange%nNodesSend(iMPINeighbor).NE.0) THEN
!    CALL MPI_WAIT(NodeExchange%SendRequest(iMPINeighbor),MPIStatus,IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!  IF(NodeExchange%nNodesRecv(iMPINeighbor).NE.0) THEN
!    CALL MPI_WAIT(NodeExchange%RecvRequest(iMPINeighbor),recv_status_list(:,iMPINeighbor),IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!END DO ! iMPINeighbor
!
!! fill list with received Node-IDs
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  IF(NodeExchange%nNodesRecv(iMPINeighbor).EQ.0) CYCLE
!  ALLOCATE(PartMPI%MPINodeNeighbor(iMPINeighbor)%RecvList(NodeExchange%nNodesRecv(iMPINeighbor)))
!  iPos=1
!  DO iRecvNode=1,NodeExchange%nNodesRecv(iMPINeighbor)
!    NativeNodeID   = INT(NodeRecvBuf(iMPINeighbor)%content(iPos))
!    IF(NativeNodeID.GT.nNodes)THEN
!     CALL abort(&
!__STAMP__&
!          ,' Cannot send halo-data to other procs. big error! ', NativeNodeID, REAL(nNodes))
!    END IF
!    PartMPI%MPINodeNeighbor(iMPINeighbor)%RecvList(iRecvNode)=NativeNodeID
!    iPos=iPos+1
!  END DO ! RecvNode=1,NodeExchange%nNodesRecv(iMPINeighbor)
!END DO ! iMPINeighbor
!
!nDOF = 1
!DO iMPINeighbor=1,PartMPI%nMPINodeNeighbors
!  SDEALLOCATE(NodeSendBuf(iMPINeighbor)%content)
!  SDEALLOCATE(NodeRecvBuf(iMPINeighbor)%content)
!  IF(NodeExchange%nNodesSend(iMPINeighbor).GT.0) THEN
!    ALLOCATE(NodeSendBuf(iMPINeighbor)%content(nDOF*NodeExchange%nNodesSend(iMPINeighbor)))
!    NodeSendBuf(iMPINeighbor)%content(:)=0.
!  END IF
!  IF(NodeExchange%nNodesRecv(iMPINeighbor).GT.0) THEN
!    ALLOCATE(NodeRecvBuf(iMPINeighbor)%content(nDOF*NodeExchange%nNodesRecv(iMPINeighbor)) )
!    NodeRecvBuf(iMPINeighbor)%content(:)=0.
!  END IF
!END DO ! iMPINeighbor
!DEALLOCATE(recv_status_list)
!
!CALL MPI_BARRIER(PartMPI%Comm,iError)
!
!
!END SUBROUTINE BuildLocNodeToHaloNodeComm
!#endif /*USE_MPI*/


!SUBROUTINE DuplicateSlavePeriodicSides()
!!===================================================================================================================================
!! duplicate periodic sides from old nPartSides=nSides to nPartSidesNew=nSides+nDuplicatePeriodicSides
!! without MPI
!! 1) loop over all sides and detect periodic sides
!! 2) increase BezierControlPoints and SideXXX from old nSides to nSides+nDuplicatePeriodicSides
!! 3) loop over the OLD sides and copy the corresponding SideXXX. The missing BezierControlPoints (periodic shifted values)
!!    are build from the other element. Now two BezierControlPoints existes which are shifted by the sideperiodicvector
!! 4) shift and map sideperiodicvector and displacement to match new sides
!! with MPI
!! 1) loop over all sides and detect periodic sides
!! 2) increase BezierControlPoints and SideXXX from old nSides to nSides+nDuplicatePeriodicSides
!! 3) loop over the OLD sides and copy the corresponding SideXXX. The missing BezierControlPoints (periodic shifted values)
!!    are build  the other element. Now two BezierControlPoints existes which are shifted by the sideperiodicvector
!! 3b) newSideId depends on localSideID and yourMPISide
!!     a) both periodic sides are on proc:
!!        * duplicate side and two separate sideids with changes in partsidetoelem
!!        * build missing side with own data
!!     b) periodic side is MPI Side
!!        I) mySide (Master)-Side
!!           *  nothing to due, old side can be reused
!!        II) yourSide (Slave)-Side
!!           *  build new Side with own data
!! 4) shift and map sideperiodicvector and displacement to match new sides
!! Note:
!! periodic sides are unique for the DG operator and duplicated for the particle tracking
!! CAUTION:
!! Routine has to be called before MarkAllBCSides
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!USE MOD_GLobals
!USE MOD_Mesh_Vars              ,ONLY: MortarType,BC,NGeo,nBCs,nSides,BoundaryType,MortarSlave2MasterInfo,nElems,XCL_NGeo
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide,PartSideToElem,nTotalSides,SidePeriodicType,nPartPeriodicSides,GEO &
!                                     ,nTotalBCSides,nPartSides
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
!USE MOD_Mesh_Vars              ,ONLY: NGeoElevated
!USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals,RotateMasterToSlave,GetBezierControlPoints3D
!USE MOD_Particle_Surfaces_vars ,ONLY: BezierControlPoints3D,SideSlabIntervals,BezierControlPoints3DElevated &
!                                        ,SideSlabIntervals,SideSlabNormals,BoundingBoxIsEmpty
!USE MOD_Particle_Tracking_Vars ,ONLY: CartesianPeriodic
!USE MOD_Particle_MPI_Vars      ,ONLY: printBezierControlPointsWarnings
!!----------------------------------------------------------------------------------------------------------------------------------!
!! insert modules here
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                              :: iSide,NBElemID,tmpnSides,NBlocSideID,ElemID,newSideID,locSideID,PVID
!INTEGER                              :: BCID,iBC,flip,ilocSide,iElem,SideID,idir
!REAL,ALLOCATABLE                     :: DummyBezierControlPoints3D(:,:,:,:)
!REAL,ALLOCATABLE                     :: DummyBezierControlPoints3DElevated(:,:,:,:)
!REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: DummySideSlabNormals                  ! normal vectors of bounding slab box
!REAL,ALLOCATABLE,DIMENSION(:,:)      :: DummySideSlabIntervals               ! intervalls beta1, beta2, beta3
!!REAL,ALLOCATABLE,DIMENSION(:,:)      :: DummySidePeriodicDisplacement        ! intervalls beta1, beta2, beta3
!LOGICAL,ALLOCATABLE,DIMENSION(:)     :: DummyBoundingBoxIsEmpty
!INTEGER,ALLOCATABLE,DIMENSION(:)     :: DummyBC
!INTEGER,ALLOCATABLE,DIMENSION(:)     :: DummyMortarSlave2MasterInfo
!INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: DummyMortarType
!INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: DummyPartSideToElem
!INTEGER,ALLOCATABLE,DIMENSION(:)     :: DummySidePeriodicType
!LOGICAL                              :: MapPeriodicSides
!REAL                                 :: MinMax(1:2),MinMaxGlob(1:6),xTest(1:3)
!!===================================================================================================================================

!! 1) loop over all sides and detect periodic sides
!nPartPeriodicSides=0
!MapPeriodicSides=.FALSE.
!IF(.NOT.CartesianPeriodic .AND. GEO%nPeriodicVectors.GT.0)THEN
!  DO iSide=1,nSides
!    IF(SidePeriodicType(iSide).NE.0)THEN
!      ! abort if particles are traced over mortar sides
!      IF(MortarSlave2MasterInfo(iSide).NE.-1.OR.MortarType(1,iSide).GE.0) THEN
!        WRITE (*,*) "MortarSlave2MasterInfo(iSide) =", MortarSlave2MasterInfo(iSide)
!        WRITE (*,*) "MortarType(1,iSide)           =", MortarType(1,iSide)
!        CALL abort(&
!          __STAMP__&
!          , ' Periodic tracing over mortar sides is not implemented!')
!      END IF
!      ! ignore MPI sides, these have NOT to be mirrored
!      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
!      IF(ElemID.EQ.-1) THEN
!        ! master side is NOT on proc, hence, the side must NOT BE DUPLICATED
!        MapPeriodicSides=.TRUE.
!        CYCLE
!      END IF
!      NBElemID=PartSideToElem(S2E_NB_ELEM_ID,iSide)
!      ! only master side is on proc, nothing to do
!      IF(NBElemID.LT.1) CYCLE
!      IF(NBElemID.GT.nElems) CYCLE
!      ! if master and slave side are on proc, duplicate
!      nPartPeriodicSides=nPartPeriodicSides+1
!      MapPeriodicSides=.TRUE.
!    END IF
!  END DO
!END IF

!!IF(nPartPeriodicSides.GT.0)THEN
!IF(MapPeriodicSides)THEN
!  ! map min-max glob to local array
!  MinMaxGlob(1)=GEO%xminglob
!  MinMaxGlob(2)=GEO%yminglob
!  MinMaxGlob(3)=GEO%zminglob
!  MinMaxGlob(4)=GEO%xmaxglob
!  MinMaxGlob(5)=GEO%ymaxglob
!  MinMaxGlob(6)=GEO%zmaxglob

!  ! 2) increase BezierControlPoints and SideXXX from old nSides to nSides+nDuplicatePeriodicSides
!  ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides))
!  ALLOCATE(DummyBezierControlPoints3dElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides))
!  ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nTotalSides))
!  ALLOCATE(DummySideSlabIntervals(1:6,1:nTotalSides))
!  !ALLOCATE(DummySidePeriodicDisplacement(1:3,1:nTotalSides))
!  ALLOCATE(DummyBoundingBoxIsEmpty(1:nTotalSides))
!  ALLOCATE(DummyBC(1:nTotalSides))
!  ALLOCATE(DummyMortarType(1:2,1:nTotalSides))
!  ALLOCATE(DummyPartSideToElem(1:5,1:nTotalSides))
!  ALLOCATE(DummySidePeriodicType(1:nTotalSides))
!  ALLOCATE(DummyMortarSlave2MasterInfo(1:nTotalSides))

!  ! copy data to backup
!  DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nTotalSides)
!  DummyBezierControlPoints3dElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides) &
!     = BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides)
!  DummySideSlabNormals(1:3,1:3,1:nTotalSides)                 = SideSlabNormals(1:3,1:3,1:nTotalSides)
!  DummySideSlabIntervals(1:6,1:nTotalSides)                   = SideSlabIntervals(1:6,1:nTotalSides)
!  !DummySidePeriodicDisplacement(1:3,1:nTotalSides)            = SidePeriodicDisplacement(1:3,1:nTotalSides)
!  DummyBoundingBoxIsEmpty(1:nTotalSides)                      = BoundingBoxIsEmpty(1:nTotalSides)
!  DummyBC(1:nTotalSides)                                      = BC(1:nTotalSides)
!  DummyMortarSlave2MasterInfo(1:nTotalSides)                  = MortarSlave2MasterInfo(1:nTotalSides)
!  DummyMortarType(1:2,1:nTotalSides)                          = MortarType(1:2,1:nTotalSides)
!  DummyPartSideTOElem(1:5,1:nTotalSides)                      = PartSideTOElem(1:5,1:nTotalSides)
!  DummySidePeriodicType(1:nTotalSides)                        = SidePeriodicType(1:nTotalSides)

!  ! deallocate old values and reallocate
!  DEALLOCATE(BezierControlPoints3D)
!  DEALLOCATE(BezierControlPoints3DElevated)
!  DEALLOCATE(SideSlabNormals)
!  DEALLOCATE(SideSlabIntervals)
!  DEALLOCATE(BoundingBoxIsEmpty)
!  DEALLOCATE(MortarSlave2MasterInfo)
!  DEALLOCATE(BC)
!  DEALLOCATE(MortarType)
!  DEALLOCATE(PartSideToElem)
!  DEALLOCATE(SidePeriodicType)

!  tmpnSides  =nTotalSides
!  nTotalSides=nTotalSides+nPartPeriodicSides
!  ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nTotalSides))
!  BezierControlPoints3D=-1.
!  ALLOCATE(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides))
!  BezierControlPoints3DElevated=-1.
!  ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalSides))
!  SideSlabNormals=-1.
!  ALLOCATE(SideSlabIntervals(1:6,1:nTotalSides))
!  SideSlabIntervals=-1.
!  ALLOCATE(BoundingBoxIsEmpty(1:nTotalSides))
!  BoundingBoxIsEmpty=.FALSE.
!  ALLOCATE(BC(1:nTotalSides))
!  BC=-3
!  ALLOCATE(MortarSlave2MasterInfo(1:nTotalSides))
!  MortarSlave2MasterInfo=-1
!  ALLOCATE(MortarType(1:2,1:nTotalSides))
!  MortarType=-1
!  ALLOCATE(PartSideToElem(1:5,1:nTotalSides))
!  PartSideToElem=-1
!  ALLOCATE(SidePeriodicType(1:nTotalSides))
!  SidePeriodicType=0
!  !ALLOCATE(SidePeriodicDisplacement(1:3,1:nTotalSides))

!  BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:tmpnSides) = DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:tmpnSides)
!  BezierControlPoints3dElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:tmpnSides) &
!     = DummyBezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:tmpnSides)
!  SideSlabNormals(1:3,1:3,1:tmpnSides)                 = DummySideSlabNormals(1:3,1:3,1:tmpnSides)
!  SideSlabIntervals(1:6,1:tmpnSides)                   = DummySideSlabIntervals(1:6,1:tmpnSides)
!  BoundingBoxIsEmpty(1:tmpnSides)                      = DummyBoundingBoxIsEmpty(1:tmpnSides)
!  !SidePeriodicDisplacement(1:3,1:nTotalSides)          = DummySidePeriodicDisplacement(1:3,1:nTotalSides)
!  BC(1:tmpnSides)                                      = DummyBC(1:tmpnSides)
!  MortarSlave2MasterInfo(1:tmpnSides)                  = DummyMortarSlave2MasterInfo(1:tmpnSides)
!  MortarType(1:2,1:tmpnSides)                          = DummyMortarType(1:2,1:tmpnSides)
!  PartSideToElem(1:5,1:tmpnSides)                      = DummyPartSideTOElem(1:5,1:tmpnSides)
!  SidePeriodicType(1:tmpnSides)                        = DummySidePeriodicType(1:tmpnSides)

!  ! 3) loop over the OLD sides and copy the corresponding SideXXX. The missing BezierControlPoints (periodic shifted values)
!  !    are build from the other element. Now two BezierControlPoints exists which are shifted by the SidePeriodicVector
!  nPartPeriodicSides=0
!  DO iSide=1,tmpnSides
!    IF(SidePeriodicType(iSide).NE.0)THEN
!      NBElemID=PartSideToElem(S2E_NB_ELEM_ID,iSide)
!      IF(NBElemID.LT.1) CYCLE
!      IF(NBElemID.GT.nElems) CYCLE
!      NBlocSideID=PartSideToElem(S2E_NB_LOC_SIDE_ID,iSide)
!      flip=PartSideToElem(S2E_FLIP,iSide)
!      locSideID=PartSideToElem(S2E_LOC_SIDE_ID,iSide)
!      ElemID   =PartSideToElem(S2E_ELEM_ID,iSide)
!      ! 3b) set newSideID and sidedata
!      IF(ElemID.EQ.-1) THEN
!        ! MPI side
!        newSideID=iSide
!        PVID=SidePeriodicType(iSide)
!        SidePeriodicType(newSideID)=-SidePeriodicType(iSide) ! stored the inital alpha value
!      ELSE
!        nPartPeriodicSides=nPartPeriodicSides+1
!        newSideID=tmpnSides+nPartPeriodicSides
!        ! bc
!        BCID = BC(iSide)
!        PVID = BoundaryType(BCID,BC_ALPHA)
!        ! loop over bc to get the NEW BC type
!        DO iBC = 1,nBCs
!          IF(BoundaryType(iBC,BC_ALPHA).EQ.-PVID) THEN
!            BC(newSideID)=iBC
!            EXIT
!          END IF
!        END DO
!        MortarSlave2MasterInfo(newSideID) = DummyMortarSlave2MasterInfo(iSide)
!        PVID=SidePeriodicType(iSide)
!        SidePeriodicType(newSideID)=-SidePeriodicType(iSide) ! stored the inital alpha value
!        ! rebuild sides for sanity
!        CALL GetBezierControlPoints3D(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID),ElemID,ilocSide_In=locSideID,SideID_In=iSide)
!        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)                             &
!                                           ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide)     &
!                                           ,SideSlabNormals(1:3,1:3,iSide)                                             &
!                                           ,SideSlabInterVals(1:6,iSide)                                               &
!                                           ,BoundingBoxIsEmpty(iSide)                                                  )
!        ! sanity check
!        ! check flip of master element, has to be zero, because of master side
!        IF(PartElemToSide(E2S_FLIP   ,locSideID,ElemID).NE.0)THEN
!          CALL abort(&
!                __STAMP__&
!                , ' Periodic Side is no master side!')
!        END IF
!        xTest(1:3) = BezierControlPoints3D(1:3,0,0,iSide)
!        xTest      = xTest + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
!        IF(xTest(1)+1e-8.LT.MinMaxGlob(1)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!        IF(xTest(2)+1e-8.LT.MinMaxGlob(2)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!        IF(xTest(3)+1e-8.LT.MinMaxGlob(3)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!        IF(xTest(1)-1e-8.GT.MinMaxGlob(4)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!        IF(xTest(2)-1e-8.GT.MinMaxGlob(5)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!        IF(xTest(3)-1e-8.GT.MinMaxGlob(6)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!      END IF
!      ! the flip has to be set to -1, artificial master side
!      PartElemToSide(E2S_FLIP   ,NBlocSideID,NBElemID) = 0
!      PartElemToSide(E2S_SIDE_ID,NBlocSideID,NBElemID) = newSideID
!      ! rebuild BezierControlPoints3D (simplified version, all BezierControlPoints3D are rebuild)
!      CALL GetBezierControlPoints3D(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,NBElemID),NBElemID,ilocSide_In=NBlocSideID,SideID_In=NewSideID)
!      ! remains equal because of MOVEMENT and MIRRORING of periodic side
!      ! periodic displacement
!      !DO q=0,NGeo
!      !  DO p=0,NGeo
!      !    BezierControlPoints3D(1:3,p,q,newSideID)  = DummyBezierControlPoints3d(1:3,p,q,iSide) &
!      !                                              + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
!      !  END DO ! p=0,NGeo
!      !END DO ! q=0,NGeo
!      !! recompute quark
!      !CALL RotateMasterToSlave(flip,NBlocSideID,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID))
!      DO idir=1,3
!        MinMax(1)=MINVAL(BezierControlPoints3D(iDir,:,:,newSideID))
!        MinMax(2)=MAXVAL(BezierControlPoints3D(iDir,:,:,newSideID))
!        ! this may be required a tolerance due to periodic displacement
!        IF(.NOT.ALMOSTEQUALRELATIVE(MinMax(1),MinMaxGlob(iDir),1e-10))THEN
!          IF(MinMax(1).LT.MinMaxGlob(iDir)) THEN
!            IPWRITE(UNIT_stdOut,*) ' Min-comparison. MinValue, GlobalMin ', MinMax(1),MinMaxGlob(iDir)
!            CALL abort(&
!                __STAMP__&
!                , ' BezierControlPoints3D is moved outside of minvalue of GEO%glob! Direction', iDir)
!          END IF
!        ELSE
!          IF(printBezierControlPointsWarnings)THEN
!            IPWRITE(UNIT_stdOut,*) ' WARNING: Min-comparison. MinValue, GlobalMin ', MinMax(1),MinMaxGlob(iDir)
!          END IF
!        END IF
!        IF(.NOT.ALMOSTEQUALRELATIVE(MinMax(2),MinMaxGlob(iDir+3),1e-10))THEN
!          IF(MinMax(2).GT.MinMaxGlob(iDir+3)) THEN
!            IPWRITE(UNIT_stdOut,*) ' Max-comparison MaxValue, GlobalMax ', MinMax(2),MinMaxGlob(iDir+3)
!            CALL abort(&
!                __STAMP__&
!                , ' BezierControlPoints3D is moved outside of maxvalue of GEO%glob! Direction', iDir)
!          END IF
!        ELSE
!          IF(printBezierControlPointsWarnings)THEN
!            IPWRITE(UNIT_stdOut,*) ' WARNING: Max-comparison MaxValue, GlobalMax ', MinMax(2),MinMaxGlob(iDir+3)
!          END IF
!        END IF
!      END DO

!      ! fill partsidetoelem
!      PartSideToElem(S2E_ELEM_ID,newSideID)=NBElemID
!      PartSideToElem(S2E_NB_ELEM_ID,newSideID)=ElemID
!      PartSideToElem(S2E_FLIP,newSideID)=-1
!      PartSideToElem(S2E_LOC_SIDE_ID,newSideID)=NBlocSideID
!      PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)=locSideID
!      ! mortar type
!      MortarType(1:2,newSideID) = DummyMortarType(1:2,iSide)
!      ! bounding box, etc...
!      CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)                         &
!                                         ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,newSideID) &
!                                         ,SideSlabNormals(1:3,1:3,newSideID)                                         &
!                                         ,SideSlabInterVals(1:6,newSideID)                                           &
!                                         ,BoundingBoxIsEmpty(newSideID)                                              )

!      ! sanity check
!      xTest(1:3) = BezierControlPoints3D(1:3,0,0,newSideID)
!      PVID=SidePeriodicType(newSideID)
!      xTest      = xTest + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
!      IF(xTest(1)+1e-8.LT.MinMaxGlob(1)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
!      IF(xTest(2)+1e-8.LT.MinMaxGlob(2)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
!      IF(xTest(3)+1e-8.LT.MinMaxGlob(3)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
!      IF(xTest(1)-1e-8.GT.MinMaxGlob(4)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
!      IF(xTest(2)-1e-8.GT.MinMaxGlob(5)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
!      IF(xTest(3)-1e-8.GT.MinMaxGlob(6)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
!    END IF
!  END DO ! iSide=1,tmpnSides
!  ! deallocate dummy
!  DEALLOCATE(DummyBezierControlPoints3D)
!  DEALLOCATE(DummySideSlabNormals)
!  DEALLOCATE(DummySideSlabIntervals)
!  DEALLOCATE(DummyBoundingBoxIsEmpty)
!  DEALLOCATE(DummyBC)
!  DEALLOCATE(DummyMortarType)
!  DEALLOCATE(DummyPartSideToElem)
!  DEALLOCATE(DummySidePeriodicType)

!END IF ! nPartPeriodicSides .GT.0
!! reset side-counter
!nPartSides     =nPartPeriodicSides+nSides
!nTotalBCSides =nPartPeriodicSides+nSides

!! sanity check for PartElemToSide
!DO iElem=1,nElems
!  DO ilocSide=1,6
!    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(MortarType(1,SideID).EQ.0)THEN
!      IF(SideID.LE.0)THEN
!        CALL abort(&
!__STAMP__&
!      , ' No Side ID set. critical error!',iElem,REAL(ilocSide))
!      END IF
!    END IF
!  END DO
!END DO ! iElem=1,PP_nElems

!#if USE_MPI
!CALL MPI_BARRIER(MPI_COMM_WORLD,iERROR)
!#endif
!SWRITE(UNIT_StdOut,'(A)') ' Sanity check of duplication successful!'

!END SUBROUTINE DuplicateSlavePeriodicSides


!SUBROUTINE MarkAllBCSides()
!!===================================================================================================================================
!! CAUTION: nTotalBCSides is reset from old value to new current,process-local BCSides
!! The PartBCSideList contains a mapping from the global side list to a local, pure BC side list
!! DG-SideList
!! 1:nBCSides - nInnerSides - nMortarSides - nMPISides
!! DG: periodic sides are no BC sides
!! ParticleTracking treats periodic sides as BC sides and the process needs the actual side,
!! hence it may be required to be duplicated (side at correct position)
!! Particle-Tracking-List before MarkAllBCSides
!! 1:nBCSides - nInnerSides - nSomePeriodicSides - nMortarSides - nMPISides - nMissingPeriodicSides
!! As RefMapping requires only the BC sides, a shorter list is generated over all
!! nTotalBCSides which is NOW smaller than nPartSides or nTotalSides
!! CAUTION:
!! This smaller list is used to build: from 1:nTotalBCSides < nTotalSides and is used for
!! SideNormVec,SideTypes,SideDistance
!! BUT: 1:nTotalSides is STILL used for
!! BezierControlPoints3D, SideSlabInterVals,SideSlabNormals,BoundingBoxIsEmpty
!! and are NOT reshaped yet, hence, the length of the array remains nTotalSides
!! CAUTION/CONTINUOUS:
!! During building of the HALO region, the BezierControlPoints variables are further increased with nTotalSides while the
!! already small arrays increases with nTotalBCSides
!! After building the HALO region, the actual arrays are reshaped and a stored in shorter arrays
!!
!! AND no rule without a break:
!! SidePeriodicType is still on nTotalSides and NOT reshaped
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!USE MOD_Mesh_Vars,               ONLY:nSides
!USE MOD_Particle_Mesh_Vars,      ONLY:PartBCSideList,nTotalSides,nPartPeriodicSides,nTotalBCSides,nPartSides
!USE MOD_Mesh_Vars,               ONLY:BC,nBCSides
!USE MOD_Globals
!!----------------------------------------------------------------------------------------------------------------------------------!
!! insert modules here
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER             :: iSide
!!===================================================================================================================================
!! Note that for DoRefMapping=T: PartBCSideList is increased, due to the periodic sides
!!           for DoRefMapping=F: A new list is created
!SDEALLOCATE(PartBCSideList)
!
!ALLOCATE(PartBCSideList(1:nTotalSides))
!! BC Sides
!PartBCSideList=-1
!nTotalBCSides=0
!DO iSide=1,nBCSides
!  nTotalBCSides=nTotalBCSides+1
!  PartBCSideList(iSide)=nTotalBCSides
!END DO
!
!DO iSide=nBCSides+1,nSides+nPartPeriodicSides
!  IF(BC(iSide).EQ.0) CYCLE
!  nTotalBCSides=nTotalBCSides+1
!  PartBCSideList(iSide)=nTotalBCSides
!END DO ! iSide
!
!! nPartsides
!nPartSides   =nPartPeriodicSides+nSides
!
!END SUBROUTINE MarkAllBCSides


!SUBROUTINE BGMIndexOfElement(ElemID,ElemToBGM)
!!===================================================================================================================================
!! computes the element indices of an given element in the BGM-mesh
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,sVdm_Bezier
!USE MOD_Particle_Surfaces_Vars ,ONLY: sVdm_Bezier
!USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo
!USE MOD_Mesh_Vars              ,ONLY: NGeo
!USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
!USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!INTEGER,INTENT(IN)        :: ElemID
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!INTEGER,INTENT(OUT)       :: ElemToBGM(1:6)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                   :: ilocSide, SideID
!REAL                      :: xmin,xmax,ymin,ymax,zmin,zmax
!REAL                      :: BezierControlPoints3D_tmp(1:3,0:NGeo,0:NGeo)
!!===================================================================================================================================
!
!xmin = HUGE(1.0)
!xmax =-HUGE(1.0)
!ymin = HUGE(1.0)
!ymax =-HUGE(1.0)
!zmin = HUGE(1.0)
!zmax =-HUGE(1.0)
!
!! get min,max of BezierControlPoints of Element
!DO iLocSide = 1,6
!  SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, ElemID)
!  IF(DoRefMapping)THEN
!    IF(SideID.GT.0)THEN
!      IF(PartElemToSide(E2S_FLIP,ilocSide,ElemID).EQ.0)THEN
!        BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
!      ELSE
!        SELECT CASE(ilocSide)
!        CASE(XI_MINUS)
!          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
!        CASE(XI_PLUS)
!          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
!        CASE(ETA_MINUS)
!          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
!        CASE(ETA_PLUS)
!          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
!        CASE(ZETA_MINUS)
!          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
!        CASE(ZETA_PLUS)
!          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
!        END SELECT
!      END IF
!    ELSE
!      SELECT CASE(ilocSide)
!      CASE(XI_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
!      CASE(XI_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
!      CASE(ETA_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
!      CASE(ETA_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
!      CASE(ZETA_MINUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
!      CASE(ZETA_PLUS)
!        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
!      END SELECT
!    END IF
!  ELSE ! pure tracing
!    BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
!  END IF
!  xmin=MIN(xmin,MINVAL(BezierControlPoints3D_tmp(1,:,:)))
!  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D_tmp(1,:,:)))
!  ymin=MIN(ymin,MINVAL(BezierControlPoints3D_tmp(2,:,:)))
!  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D_tmp(2,:,:)))
!  zmin=MIN(zmin,MINVAL(BezierControlPoints3D_tmp(3,:,:)))
!  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D_tmp(3,:,:)))
!END DO ! ilocSide
!
!ElemToBGM(1) = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
!ElemToBGM(2) = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
!ElemToBGM(3) = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
!ElemToBGM(4) = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
!ElemToBGM(5) = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
!ElemToBGM(6) = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
!
!
!END SUBROUTINE BGMIndexOfElement




!SUBROUTINE GetSideOriginAndRadius(nTotalBCSides,SideOrigin,SideRadius)
!!===================================================================================================================================
!! ONLY RefMapping
!! Computes the side origin and radius for each BC Side
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Mesh_Vars              ,ONLY: NGeo
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartBCSideList,nTotalSides
!USE MOD_Basis                  ,ONLY: DeCasteljauInterpolation
!USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!INTEGER,INTENT(IN)       :: nTotalBCSides
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)         :: SideOrigin(1:3,1:nTotalBCSides),SideRadius(1:nTotalBCSides)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                  :: iSide, BCSideID,p,q
!REAL                     :: Xi(1:2), Origin(1:3), Radius, RadiusMax, Vec(1:3)
!!===================================================================================================================================
!
!SideOrigin=0.
!SideRadius=0.
!
!DO iSide=1,nTotalSides
!  BCSideID=PartBCSideList(iSide)
!  IF(BCSideID.LT.1) CYCLE
!  Xi=0.
!  CALL DeCasteljauInterpolation(NGeo,Xi(1:2),BCSideID,Origin)
!  SideOrigin(1:3,BCSideID) = Origin
!  Radius=0.
!  RadiusMax=0.
!  DO q=0,NGeo
!    DO p=0,NGeo
!      Vec(1:3) = BezierControlPoints3D(:,p,q,BCSideID)-Origin
!      Radius=DOT_PRODUCT(Vec,Vec)
!      RadiusMax=MAX(RadiusMax,Radius)
!    END DO ! p=0,NGeo
!  END DO ! q=0,NGeo
!  SideRadius(BCSideID)=SQRT(RadiusMax)
!END DO ! iSide=1,nTotalSides
!
!END SUBROUTINE GetSideOriginAndRadius


!SUBROUTINE GetElemToSideDistance(nTotalBCSides,SideOrigin,SideRadius)
!!===================================================================================================================================
!! computes the distance between each element and it associated sides for DoRefMapping=T
!! only sides for which ElemToSideDistance<lengthPartTrajectory have to be checked during the current tracing step
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Preproc
!USE MOD_Mesh_Vars          ,ONLY: ElemBaryNGeo
!USE MOD_Particle_Mesh_Vars ,ONLY: ElemBC_Shared,ElemRadiusNGeo,BCElem,PartBCSideList,nTotalElems
!USE MOD_Utils              ,ONLY: InsertionSort
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!INTEGER,INTENT(IN)       :: nTotalBCSides
!REAL,INTENT(IN)          :: SideOrigin(1:3,1:nTotalBCSides),SideRadius(1:nTotalBCSides)
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                  :: iElem,ilocSide,SideID,BCSideID
!REAL                     :: Vec(1:3)
!REAL                     :: Origin(1:3)
!!===================================================================================================================================
!
!! loop over all  elements
!DO iElem=1,nTotalElems
!  IF(.NOT.IsTracingBCElem(iElem)) CYCLE
!  ALLOCATE( BCElem(iElem)%ElemToSideDistance(BCElem(iElem)%lastSide) )
!  BCElem(iElem)%ElemToSideDistance(BCElem(iElem)%lastSide)=0.
!  Origin(1:3) = ElemBaryNGeo(1:3,iElem)
!  ! loop over all associated sides
!  DO iLocSide=1,BCElem(iElem)%lastSide
!    SideID=BCElem(iElem)%BCSideID(ilocSide)
!    BCSideID=PartBCSideList(SideID)
!    Vec=Origin - SideOrigin(1:3,BCSideID)
!    BCElem(iElem)%ElemToSideDistance(ilocSide) = SQRT(DOT_PRODUCT(Vec,Vec))-ElemRadiusNGeo(iElem)-SideRadius(BCSideID)
!  END DO ! iLocSide=1,BCElem(iElem)%lastSide
!  ! sort each side distance for each element according to it's distance
!  CALL InsertionSort(BCElem(iElem)%ElemToSideDistance(:),BCElem(iElem)%BCSideID(:),BCElem(iElem)%lastSide)
!END DO ! iElem=1,PP_nElems
!
!END SUBROUTINE GetElemToSideDistance


SUBROUTINE MarkAuxBCElems()
!===================================================================================================================================
! check if auxBCs are inside BoundingBox of Elems
! -- plane: use plane equation f=a1*x+a2*y+a3*z+a4=0 and insert corresponding intervals of box -> fmin and fmax
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemHasAuxBCs,GEO
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


!#if USE_MPI
! !===================================================================================================================================
! !> For each rank an ElemData array 'ElemHaloInfo' is created, which contains information regarding the halo region of each rank
! !> ElemHaloInfo = 0: element not in list
! !>          = 1: local element
! !>          = 2: halo element
! !===================================================================================================================================
! SUBROUTINE SetHaloInfo()
! ! MODULES                                                                                                                          !
! USE MOD_GLobals
! USE MOD_Preproc            ,ONLY: PP_nElems
! USE MOD_Particle_Mesh_Vars ,ONLY: ElemHaloInfoProc
! USE MOD_Particle_MPI_Vars  ,ONLY: PartHaloElemToProc,PartMPI
! USE MOD_Particle_Mesh_Vars ,ONLY: nTotalElems
! USE MOD_IO_HDF5            ,ONLY: AddToElemData,ElementOut
! !----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT NONE
! !----------------------------------------------------------------------------------------------------------------------------------!
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! TYPE tMPIMessage
!   REAL,ALLOCATABLE               :: content(:)            ! Message buffer real
!   LOGICAL,ALLOCATABLE            :: content_log(:)        ! Message buffer logical for BGM
!   INTEGER,ALLOCATABLE            :: content_int(:)        ! Message buffer for integer for adsorption
! END TYPE
!
! TYPE tHaloInfoMPIExchange
!   INTEGER,ALLOCATABLE            :: nHaloElemsSend(:,:)   ! Only MPI neighbors
!   INTEGER,ALLOCATABLE            :: nHaloElemsRecv(:,:)   ! Only MPI neighbors
!   INTEGER                        :: nMPIHaloReceivedElems ! Number of all received particles
!   INTEGER,ALLOCATABLE            :: SendRequest(:,:)      ! Send request message handle 1 - Number, 2-Message
!   INTEGER,ALLOCATABLE            :: RecvRequest(:,:)      ! Receive request message handle,  1 - Number, 2-Message
!   TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)       ! Message, required for particle emission
! END TYPE
!  TYPE (tHaloInfoMPIExchange)     :: HaloInfoMPIExchange   ! Exchange of halo element information between ranks for ElemData output
!
! TYPE(tMPIMessage),ALLOCATABLE    :: HaloInfoRecvBuf(:)    ! HaloInfoRecvBuf with all required types
! TYPE(tMPIMessage),ALLOCATABLE    :: HaloInfoSendBuf(:)    ! HaloInfoSendBuf with all required types
!
! INTEGER,ALLOCATABLE              :: HaloElemTargetProc(:) ! Local rank id for communication
! INTEGER                          :: nSendHaloElems        ! Number of halo elements in HaloInfoMPIExchange%nHaloElemsSend(1,iProc)
! INTEGER                          :: nRecvHaloElems        ! Number of halo elements in HaloInfoMPIExchange%nHaloElemsRecv(1,iProc)
!
!
! CHARACTER(32)                    :: hilf                  ! Auxiliary variable
! INTEGER                          :: yourrank,myelem,iProc,iPos,jPos,messagesize,iElem,ALLOCSTAT,iRank,yourelemid
! INTEGER,PARAMETER                :: HaloInfoCommSize=3
! INTEGER                          :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
! !===================================================================================================================================
! ! Allocate type array for all ranks
! ALLOCATE(ElemHaloInfoProc(0:nProcessors-1))
!
! ! Allocate for local elements: Container with information of my local elements and your halo elements
! DO iRank = 0, nProcessors-1
!   ALLOCATE(ElemHaloInfoProc(iRank)%ElemHaloInfo(PP_nElems))
!   ElemHaloInfoProc(iRank)%ElemHaloInfo = 0 ! Elements that do not belong to the processor and are not halo elements are marked "0"
! END DO ! iRank = 1, nProcessors
!
! ! Mark each local element with its ID
! DO iElem = 1, PP_nElems
!   ElemHaloInfoProc(myrank)%ElemHaloInfo(iElem) = iElem
! END DO ! iElem = 1, PP_nElems
!
! ! Add arrays to ElemData
! DO iRank = 0, nProcessors-1
!   WRITE(UNIT=hilf,FMT='(I0)') iRank ! myrank
!   CALL AddToElemData(ElementOut,'MyRank'//TRIM(hilf)//'_ElemHaloInfo',IntArray=ElemHaloInfoProc(iRank)%ElemHaloInfo)
! END DO ! iRank = 1, nProcessors
!
! ! Allocate halo info arrays
! ALLOCATE( HaloInfoMPIExchange%nHaloElemsSend(2,PartMPI%nMPINeighbors)  &
!         , HaloInfoMPIExchange%nHaloElemsRecv(2,PartMPI%nMPINeighbors)  &
!         , HaloInfoRecvBuf(1:PartMPI%nMPINeighbors)                 &
!         , HaloInfoSendBuf(1:PartMPI%nMPINeighbors)                 &
!         , HaloInfoMPIExchange%SendRequest(2,PartMPI%nMPINeighbors) &
!         , HaloInfoMPIExchange%RecvRequest(2,PartMPI%nMPINeighbors) &
!         , HaloElemTargetProc(PP_nElems+1:nTotalElems)              &
!         , STAT=ALLOCSTAT                                       )
!
! IF (ALLOCSTAT.NE.0) CALL abort(&
!     __STAMP__&
!     ,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)
!
! HaloInfoMPIExchange%nHaloElemsSend=0
! HaloInfoMPIExchange%nHaloElemsRecv=0
!
!
!
!
!
! ! Communicate halo elem info
! !===================================================================================================================================
! ! 1 of 4: SUBROUTINE IRecvNbOfParticles()
! !===================================================================================================================================
! HaloInfoMPIExchange%nHaloElemsRecv=0
! DO iProc=1,PartMPI%nMPINeighbors
!   CALL MPI_IRECV( HaloInfoMPIExchange%nHaloElemsRecv(:,iProc) &
!                 , 2                                           &
!                 , MPI_INTEGER                                 &
!                 , PartMPI%MPINeighbor(iProc)                  &
!                 , 1001                                        &
!                 , PartMPI%COMM                                &
!                 , HaloInfoMPIExchange%RecvRequest(1,iProc)    &
!                 , IERROR )
! END DO ! iProc
!
!
!
!
! !===================================================================================================================================
! ! 2 of 4: SUBROUTINE SendNbOfParticles(doParticle_In)
! !===================================================================================================================================
! ! 1) get number of send particles
! HaloInfoMPIExchange%nHaloElemsSend=0
! HaloElemTargetProc=-1
! DO iElem = PP_nElems+1, nTotalElems
!   ! Send number of halo elements to each proc
!   HaloInfoMPIExchange%nHaloElemsSend(1,PartHaloElemToProc(LOCAL_PROC_ID,iElem)) = &
!       HaloInfoMPIExchange%nHaloElemsSend(1,PartHaloElemToProc(LOCAL_PROC_ID,iElem)) + 1
!   ! Set target iProc of the element for setting the particle message that is sent
!   HaloElemTargetProc(iElem) = PartHaloElemToProc(LOCAL_PROC_ID,iElem)
! END DO ! iElem = PP_nElems+1, nTotalElems
!
! ! 2) send number of particles
! DO iProc=1,PartMPI%nMPINeighbors
!   CALL MPI_ISEND( HaloInfoMPIExchange%nHaloElemsSend(:,iProc) &
!                 , 2                                           &
!                 , MPI_INTEGER                                 &
!                 , PartMPI%MPINeighbor(iProc)                  &
!                 , 1001                                        &
!                 , PartMPI%COMM                                &
!                 , HaloInfoMPIExchange%SendRequest(1,iProc)    &
!                 , IERROR )
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!     __STAMP__&
!     ,' MPI Communication error', IERROR)
! END DO ! iProc
!
!
!
!
!
! !===================================================================================================================================
! ! 3 of 4: SUBROUTINE MPIParticleSend()
! !===================================================================================================================================
! ! 3) Build Message
! DO iProc=1, PartMPI%nMPINeighbors
!   ! allocate SendBuf
!   nSendHaloElems=HaloInfoMPIExchange%nHaloElemsSend(1,iProc)
!   iPos=0
!   MessageSize=nSendHaloElems*HaloInfoCommSize
!
!   ALLOCATE(HaloInfoSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
!   IF (ALLOCSTAT.NE.0) CALL abort(&
!   __STAMP__&
!   ,'  Cannot allocate HaloInfoSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))
!
!   ! fill message
!   DO iElem = PP_nElems+1, nTotalElems
!     ! Element is element with target proc-id equals local proc id
!     IF(HaloElemTargetProc(iElem).NE.iProc) CYCLE
!     ! my rank
!     HaloInfoSendBuf(iProc)%content(1+iPos) = REAL(myrank,KIND=8)
!     jPos=iPos+1
!
!     ! local element ID of new host proc: PEM%GlobalElemID(PartID)
!     HaloInfoSendBuf(iProc)%content(    1+jPos)    = REAL(PartHaloElemToProc(NATIVE_ELEM_ID,iElem),KIND=8)
!     jPos=jPos+1
!
!     ! local element ID of old proc (halo elem id)
!     HaloInfoSendBuf(iProc)%content(    1+jPos)    = REAL(iElem,KIND=8)
!     jPos=jPos+1
!
!     IF(MOD(jPos,HaloInfoCommSize).NE.0) THEN
!       IPWRITE(UNIT_stdOut,*)  'HaloInfoCommSize',HaloInfoCommSize
!       IPWRITE(UNIT_stdOut,*)  'jPos',jPos
!       CALL Abort(&
!           __STAMP__&
!           ,' CalcHaloInfo: wrong sending message size!')
!     END IF
!     iPos=iPos+HaloInfoCommSize
!   END DO  ! iElem = PP_nElems+1, nTotalElems
! END DO ! iProc
!
!
!
! ! 4) Finish Received number of halo elements
! DO iProc=1,PartMPI%nMPINeighbors
!   CALL MPI_WAIT(HaloInfoMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!     __STAMP__&
!     ,' MPI Communication error', IERROR)
!   CALL MPI_WAIT(HaloInfoMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!     __STAMP__&
!     ,' MPI Communication error', IERROR)
! END DO ! iProc
!
! ! total number of received particles: add up number of all received ranks
! HaloInfoMPIExchange%nMPIHaloReceivedElems=SUM(HaloInfoMPIExchange%nHaloElemsRecv(1,:))
!
!
!
! ! 5) Allocate received buffer and open MPI_IRECV
! DO iProc=1,PartMPI%nMPINeighbors
!   nRecvHaloElems=HaloInfoMPIExchange%nHaloElemsRecv(1,iProc)
!   MessageSize=nRecvHaloElems*HaloInfoCommSize
!   ALLOCATE(HaloInfoRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
!   IF (ALLOCSTAT.NE.0) THEN
!     IPWRITE(*,*) 'sum of total received particles            ', SUM(HaloInfoMPIExchange%nHaloElemsRecv(1,:))
!     IPWRITE(*,*) 'sum of total received deposition particles ', SUM(HaloInfoMPIExchange%nHaloElemsRecv(2,:))
!     CALL abort(&
!     __STAMP__&
!     ,'  Cannot allocate HaloInfoRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
!   END IF
!   CALL MPI_IRECV( HaloInfoRecvBuf(iProc)%content                                 &
!                 , MessageSize                                                &
!                 , MPI_DOUBLE_PRECISION                                       &
!                 , PartMPI%MPINeighbor(iProc)                                 &
!                 , 1002                                                       &
!                 , PartMPI%COMM                                               &
!                 , HaloInfoMPIExchange%RecvRequest(2,iProc)                       &
!                 , IERROR )
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!     __STAMP__&
!     ,' MPI Communication error', IERROR)
! END DO ! iProc
!
! ! 6) Send halo elements
! DO iProc=1,PartMPI%nMPINeighbors
!   nSendHaloElems = HaloInfoMPIExchange%nHaloElemsSend(1,iProc)
!   MessageSize    = nSendHaloElems*HaloInfoCommSize
!   CALL MPI_ISEND( HaloInfoSendBuf(iProc)%content                             &
!                 , MessageSize                                                &
!                 , MPI_DOUBLE_PRECISION                                       &
!                 , PartMPI%MPINeighbor(iProc)                                 &
!                 , 1002                                                       &
!                 , PartMPI%COMM                                               &
!                 , HaloInfoMPIExchange%SendRequest(2,iProc)                       &
!                 , IERROR )
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!     __STAMP__&
!     ,' MPI Communication error', IERROR)
! END DO ! iProc
!
!
!
!
!
!
! !===================================================================================================================================
! ! 4 of 4: SUBROUTINE MPIParticleRecv()
! !===================================================================================================================================
! DO iProc=1,PartMPI%nMPINeighbors
!   CALL MPI_WAIT(HaloInfoMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!     __STAMP__&
!     ,' MPI Communication error', IERROR)
! END DO ! iProc
!
! DO iProc=1,PartMPI%nMPINeighbors
!   nRecvHaloElems=HaloInfoMPIExchange%nHaloElemsRecv(1,iProc)
!   MessageSize=nRecvHaloElems*HaloInfoCommSize
!   ! finish communication with iproc
!   CALL MPI_WAIT(HaloInfoMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
!   ! Evaluate the received data and assign the halo information to the local elements
!   DO iPos=0,MessageSize-1,HaloInfoCommSize
!     IF(nRecvHaloElems.EQ.0) EXIT
!     yourrank   = INT(HaloInfoRecvBuf(iProc)%content( 1+iPos),KIND=4)
!     jPos=iPos+1
!
!     myelem     = INT(HaloInfoRecvBuf(iProc)%content( 1+jPos),KIND=4)
!     jPos=jPos+1
!
!     yourelemid     = INT(HaloInfoRecvBuf(iProc)%content( 1+jPos),KIND=4)
!     jPos=jPos+1
!
!     IF(MOD(jPos,HaloInfoCommSize).NE.0)THEN
!       IPWRITE(UNIT_stdOut,*)  'HaloInfoCommSize',HaloInfoCommSize
!       IPWRITE(UNIT_stdOut,*)  'jPos',jPos
!       CALL Abort(&
!           __STAMP__&
!           ,' HaloInfoCommSize-wrong receiving message size!')
!     END IF
!
!     ! Set halo info: halo elements are marked with "-yourelemid" (negative ElemID of the halo region element)
!     ElemHaloInfoProc(yourrank)%ElemHaloInfo(myelem) = -yourelemid
!   END DO
! END DO ! iProc
!
!
! ! deallocate send,receive buffer
! DO iProc=1,PartMPI%nMPINeighbors
!   SDEALLOCATE(HaloInfoRecvBuf(iProc)%content)
!   SDEALLOCATE(HaloInfoSendBuf(iProc)%content)
! END DO ! iProc
!
!
! ! De-allocate halo info arrays
! SDEALLOCATE(HaloInfoMPIExchange%nHaloElemsSend)
! SDEALLOCATE(HaloInfoMPIExchange%nHaloElemsRecv)
! SDEALLOCATE(HaloInfoRecvBuf)
! SDEALLOCATE(HaloInfoSendBuf)
! SDEALLOCATE(HaloInfoMPIExchange%SendRequest)
! SDEALLOCATE(HaloInfoMPIExchange%RecvRequest)
! SDEALLOCATE(HaloElemTargetProc)
!
! END SUBROUTINE SetHaloInfo
!#endif /*USE_MPI*/

END MODULE MOD_Particle_Mesh
