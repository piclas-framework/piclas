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
PUBLIC
SAVE
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

INTERFACE InitFIBGM
  MODULE PROCEDURE InitFIBGM
END INTERFACE

INTERFACE SingleParticleToExactElement
  MODULE PROCEDURE SingleParticleToExactElement
END INTERFACE

INTERFACE SingleParticleToExactElementNoMap
  MODULE PROCEDURE SingleParticleToExactElementNoMap
END INTERFACE

INTERFACE MapRegionToElem
  MODULE PROCEDURE MapRegionToElem
END INTERFACE

INTERFACE PointToExactElement
  MODULE PROCEDURE PointToExactElement
END INTERFACE

INTERFACE BuildElementBasis
  MODULE PROCEDURE BuildElementBasis
END INTERFACE

INTERFACE CountPartsPerElem
  MODULE PROCEDURE CountPartsPerElem
END INTERFACE

!INTERFACE CheckIfCurvedElem
!  MODULE PROCEDURE CheckIfCurvedElem
!END INTERFACE

INTERFACE InitElemBoundingBox
  MODULE PROCEDURE InitElemBoundingBox
END INTERFACE

INTERFACE InsideElemBoundingBox
  MODULE PROCEDURE InsideElemBoundingBox
END INTERFACE

INTERFACE GetElemAndSideType
  MODULE PROCEDURE GetElemAndSideType
END INTERFACE

INTERFACE PartInElemCheck
  MODULE PROCEDURE PartInElemCheck
END INTERFACE

INTERFACE ParticleInsideQuad3D
  MODULE PROCEDURE ParticleInsideQuad3D
END INTERFACE

INTERFACE MarkAuxBCElems
  MODULE PROCEDURE MarkAuxBCElems
END INTERFACE

INTERFACE BoundsOfElement
  MODULE PROCEDURE BoundsOfElement
END INTERFACE

PUBLIC::CountPartsPerElem
PUBLIC::BuildElementBasis,CheckIfCurvedElem
PUBLIC::MapRegionToElem,PointToExactElement
PUBLIC::InitParticleMesh,FinalizeParticleMesh, InitFIBGM, SingleParticleToExactElement, SingleParticleToExactElementNoMap
PUBLIC::InsideElemBoundingBox
PUBLIC::PartInElemCheck
PUBLIC::ParticleInsideQuad3D
PUBLIC::InitParticleGeometry
PUBLIC::MarkAuxBCElems
PUBLIC::BoundsOfElement
!===================================================================================================================================
!
PUBLIC::DefineParametersParticleMesh
CONTAINS

!==================================================================================================================================
!> Define parameters for particle tracking
!==================================================================================================================================
SUBROUTINE DefineParametersParticleMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('Tracking')

CALL prms%CreateLogicalOption( 'DoRefMapping'&
  , 'Refmapping [T] or Tracing [F] algorithms are used for tracking of particles.'&
  , '.TRUE.')

CALL prms%CreateLogicalOption( 'TriaTracking'&
  , 'Using Triangle-aproximation [T] or (bi-)linear and bezier (curved) description [F] of sides for tracing algorithms.'//&
  ' Currently flag is only used in DSMC timediscs. Requires DoRefMapping=F.'&
  ,'.FALSE.')
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
CALL prms%CreateRealOption(    'BezierEpsilonBilinear'&
    , ' Bi-linear tolerance for the bi-linear - planar decision.' , '1e-6')
CALL prms%CreateIntOption(     'BezierElevation'&
  , ' Use BezierElevation>0 to tighten the bounding box. Typical values>10','0')
CALL prms%CreateIntOption(     'BezierSampleN'&
  , 'TODO-DEFINE-PARAMETER\n'//&
    'Default value: NGeo equidistant sampling of bezier surface for emission','0')


! Background mesh init variables
CALL prms%CreateRealArrayOption('Part-FIBGMdeltas'&
  , 'Define the deltas for the cartesian Fast-Init-Background-Mesh.'//&
  ' They should be of the similar size as the smallest cells of the used mesh for simulation.'&
  , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('Part-FactorFIBGM'&
  , 'Factor with which the background mesh will be scaled.'&
  , '1. , 1. , 1.')
CALL prms%CreateLogicalOption( 'printMPINeighborWarnings'&
    ,  ' Print warning if the MPI-Halo-region between to procs are not overlapping. Only one proc find the other in halo ' &
    ,'.FALSE.')
CALL prms%CreateLogicalOption( 'printBezierControlPointsWarnings'&
    ,  ' Print warning if MINVAL(BezierControlPoints3d(iDir,:,:,newSideID)) and global boundaries are too close ' &
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
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars, ONLY:BezierEpsilonBilinear,BezierElevation,BezierControlPoints3DElevated
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping,MeasureTrackTime,FastPeriodic,CountNbOfLostParts,nLostParts,CartesianPeriodic
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking, WriteTriaDebugMesh
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars, ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
USE MOD_Mesh_Vars,              ONLY:nElems,nSides,nNodes,SideToElem,ElemToSide,NGeo,NGeoElevated,OffSetElem,ElemToElemGlob
USE MOD_ReadInTools,            ONLY:GETREAL,GETINT,GETLOGICAL,GetRealArray
USE MOD_Particle_Surfaces_Vars, ONLY:BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux,WriteTriaSurfaceFluxDebugMesh
USE MOD_Mesh_Vars,              ONLY:useCurveds,NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ALLOCSTAT,RefMappingGuessProposal
INTEGER           :: iElem, ilocSide,iSide,iSample,ElemIDGlob
CHARACTER(LEN=2)  :: hilf
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH ...'
IF(ParticleMeshInitIsDone) CALL abort(&
__STAMP__&
, ' Particle-Mesh is already initialized.')
! allocate and duplicate partElemToside
nTotalSides=nSides
nTotalBCSides=nSides
nTotalElems=nElems
nTotalNodes=nNodes
ALLOCATE(PartElemToSide(1:2,1:6,1:nTotalSides)    &
        ,PartSideToElem(1:5,1:nTotalSides)        &
        ,PartElemToElemGlob(1:4,1:6,1:nTotalElems)&
        ,STAT=ALLOCSTAT                      )
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate particle mesh vars!')
! nullify
PartElemToSide=-1
PartSideToElem=-1
PartElemToElemGlob=-1

DoRefMapping       = GETLOGICAL('DoRefMapping',".TRUE.")
TriaTracking       = GETLOGICAL('TriaTracking','.FALSE.')

IF ((DoRefMapping.OR.UseCurveds.OR.(NGeo.GT.1)).AND.(TriaTracking)) THEN
  CALL abort(&
__STAMP__&
,'DoRefMapping=T .OR. UseCurveds=T .OR. NGEO>1! Not possible with TriaTracking=T at the same time!')
ELSE IF (TriaTracking) THEN
  WriteTriaDebugMesh = GETLOGICAL('Write-Tria-DebugMesh','.FALSE.')
ELSE
  WriteTriaDebugMesh = .FALSE.
END IF
CountNbOfLostParts = GETLOGICAL('CountNbOfLostParts',".FALSE.")
nLostParts         = 0 

#ifdef CODE_ANALYZE
PARTOUT            = GETINT('PartOut','0')
MPIRankOut         = GETINT('MPIRankOut','0')
#endif /*CODE_ANALYZE*/

!IF(.NOT.DoRefMapping) THEN
!  SDEALLOCATE(nTracksPerElem)
!END IF
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
WRITE(hilf,'(I2.2)') RefMappingGuessProposal
RefMappingGuess = GETINT('RefMappingGuess',hilf)
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

BezierEpsilonBilinear = GETREAL('BezierEpsilonBilinear','1e-6')

BezierElevation = GETINT('BezierElevation','0')
NGeoElevated    = NGeo + BezierElevation
SDEALLOCATE(BezierControlPoints3DElevated)
ALLOCATE(BezierControlPoints3DElevated(1:3,0:NGeo+BezierElevation,0:NGeo+BezierElevation,1:nSides) &
        ,STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate BezierControlPoints3DElevated!')
BezierControlPoints3DElevated=0.

! BezierAreaSample stuff:
WRITE(hilf,'(L1)') TriaTracking
TriaSurfaceFlux = GETLOGICAL('TriaSurfaceFlux',TRIM(hilf))
IF (TriaSurfaceFlux) THEN
  BezierSampleN = 1
  SurfFluxSideSize=(/1,2/)
  WriteTriaSurfaceFluxDebugMesh = GETLOGICAL('Write-TriaSurfaceFlux-DebugMesh','.FALSE.')
ELSE
  WRITE(hilf,'(I2.2)') NGeo
  BezierSampleN = GETINT('BezierSampleN',hilf)
  WriteTriaSurfaceFluxDebugMesh=.FALSE.
  SurfFluxSideSize=BezierSampleN
  ALLOCATE(BezierSampleXi(0:BezierSampleN))!,STAT=ALLOCSTAT)
  DO iSample=0,BezierSampleN
    BezierSampleXi(iSample)=-1.+2.0/BezierSampleN*iSample
  END DO
END IF

! copy
DO iElem=1,PP_nElems
  DO iLocSide=1,6
    PartElemToSide(:,iLocSide,iElem)=ElemToSide(:,iLocSide,iElem)
  END DO 
  ElemIDGlob=OffSetElem+iElem
  PartElemToElemGlob(1:4,1:6,iElem)=ElemToElemGlob(1:4,1:6,ElemIDGlob)
END DO
DO iSide=1,nSides
  PartSideToElem(:,iSide)=SideToElem(:,iSide)
END DO 


ParticleMeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
! Subroutine for particle geometry initialization (GEO container)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nElems, nNodes
USE MOD_Mesh_Vars              ,ONLY: Elems, offsetElem, ElemToSide
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Tracking_Vars ,ONLY: WriteTriaDebugMesh
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem, iLocSide, iNode, jNode
INTEGER            :: nStart, NodeNum
INTEGER            :: ALLOCSTAT
INTEGER            :: NodeMap(4,6),nSides
REAL               :: A(3,3),detcon
REAL,ALLOCATABLE   :: Coords(:,:,:,:)
CHARACTER(32)      :: hilf
CHARACTER(LEN=255) :: FileString
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
         GEO%ConcaveElemSide(1:6,1:nElems), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(__STAMP__&
 ,'ERROR in InitParticleGeometry: Cannot allocate GEO%... stuff!')
END IF

!ALLOCATE(GEO%ElemToNodeIDGlobal(1:8,1:nElems))

GEO%ElemToNodeID(:,:)=0
GEO%ElemSideNodeID(:,:,:)=0
GEO%NodeCoords(:,:)=0.
GEO%ConcaveElemSide(:,:)=.FALSE.
iNode=0
DO iElem=1,nElems
  DO jNode=1,8
    Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=0
  END DO
END DO
DO iElem=1,nElems
  !--- Save corners of sides
  DO jNode=1,8
    IF (Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID.EQ.0) THEN
      iNode=iNode+1
      Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=iNode
      GEO%NodeCoords(1:3,iNode)=Elems(iElem+offsetElem)%ep%node(jNode)%np%x(1:3)
    END IF
    GEO%ElemToNodeID(jNode,iElem)=Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID
    !GEO%ElemToNodeIDGlobal(jNode,iElem) = Elems(iElem+offsetElem)%ep%node(jNode)%np%ind
  END DO
END DO

DO iElem=1,nElems
  DO iLocSide=1,6
    nStart=MAX(0,ElemToSide(E2S_FLIP,iLocSide,iElem)-1)
    GEO%ElemSideNodeID(1:4,iLocSide,iElem)=(/Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart  ,4)+1,iLocSide))%np%NodeID,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+1,4)+1,iLocSide))%np%NodeID,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+2,4)+1,iLocSide))%np%NodeID,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+3,4)+1,iLocSide))%np%NodeID/)
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

!-- write debug-mesh
IF (WriteTriaDebugMesh) THEN
  nSides=6
  WRITE(UNIT=hilf,FMT='(I4.4)') myRank
  FileString='TRIA-DebugMesh_PROC'//TRIM(hilf)//'.vtu'
  ALLOCATE(Coords(1:3,1:4,1:nSides,1:nElems))
  DO iElem = 1,nElems ; DO iLocSide = 1,nSides ; DO iNode = 1,4
    Coords(:,iNode,iLocSide,iElem)=GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,iLocSide,iElem))
  END DO ; END DO ; END DO
  CALL WriteTriaDataToVTK(nSides,nElems,Coords(1:3,1:4,1:6,1:nElems),FileString)
  SDEALLOCATE(Coords)
END IF !WriteTriaDebugMesh

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


SUBROUTINE FinalizeParticleMesh()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars             , ONLY: nElems, nNodes
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Tracking_Vars, ONLY: Distance,ListDistance
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

SDEALLOCATE(PartElemToSide)
SDEALLOCATE(PartSideToElem)
SDEALLOCATE(PartElemToElemGlob)
SDEALLOCATE(PartElemToElemAndSide)
SDEALLOCATE(PartBCSideList)
SDEALLOCATE(SidePeriodicType)
SDEALLOCATE(SidePeriodicDisplacement)
SDEALLOCATE(IsTracingBCElem)
SDEALLOCATE(TracingBCInnerSides)
SDEALLOCATE(TracingBCTotalSides)
SDEALLOCATE(ElemType)
SDEALLOCATE(GEO%PeriodicVectors)
SDEALLOCATE(GEO%PeriodicVectorsLength)
SDEALLOCATE(GEO%FIBGM)
SDEALLOCATE(GEO%Volume)
SDEALLOCATE(GEO%CharLength)
SDEALLOCATE(GEO%DeltaEvMPF)
SDEALLOCATE(GEO%ElemToFIBGM)
SDEALLOCATE(GEO%TFIBGM)

SDEALLOCATE(GEO%ElemToNodeID)
SDEALLOCATE(GEO%ElemSideNodeID)
!SDEALLOCATE(GEO%ElemToNodeIDGlobal)
SDEALLOCATE(GEO%NodeCoords)
SDEALLOCATE(GEO%ConcaveElemSide)
SDEALLOCATE(GEO%ElemsOnNode)
SDEALLOCATE(GEO%NeighNodesOnNode)
SDEALLOCATE(GEO%NumNeighborElems)
IF (ALLOCATED(GEO%ElemToNeighElems)) THEN
  DO iElem=1,nElems
    SDEALLOCATE(GEO%ElemToNeighElems(iElem)%ElemID)
  END DO
END IF
SDEALLOCATE(GEO%ElemToNeighElems)
IF (ALLOCATED(GEO%NodeToElem)) THEN
  DO iNode=1,nNodes
    SDEALLOCATE(GEO%NodeToElem(iNode)%ElemID)
  END DO
END IF
SDEALLOCATE(GEO%NodeToElem)
IF (ALLOCATED(GEO%NodeToNeighNode)) THEN
  DO iNode=1,nNodes
    SDEALLOCATE(GEO%NodeToNeighNode(iNode)%ElemID)
  END DO
END IF
SDEALLOCATE(GEO%NodeToNeighNode)

SDEALLOCATE(BCElem)
SDEALLOCATE(XiEtaZetaBasis)
SDEALLOCATE(slenXiEtaZetaBasis)
SDEALLOCATE(ElemRadiusNGeo)
SDEALLOCATE(ElemRadius2NGeo)
SDEALLOCATE(EpsOneCell)
SDEALLOCATE(Distance)
SDEALLOCATE(ListDistance)
SDEALLOCATE(isTracingTrouble)
SDEALLOCATE(ElemTolerance)
SDEALLOCATE(ElemToGlobalElemID)

ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh


SUBROUTINE SingleParticleToExactElement(iPart,doHalo,initFix,doRelocate)
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,               ONLY:PartState,PEM,PDM,PartPosRef,KeepWallParticles
USE MOD_Particle_Mesh_Vars,          ONLY:Geo
USE MOD_Particle_Tracking_Vars,      ONLY:DoRefMapping,TriaTracking
USE MOD_Particle_Mesh_Vars,          ONLY:epsOneCell,IsTracingBCElem,ElemRadius2NGeo
USE MOD_Eval_xyz,                    ONLY:GetPositionInRefElem
USE MOD_Utils,                       ONLY:InsertionSort !BubbleSortID
USE MOD_Particle_Tracking_Vars,      ONLY:DoRefMapping,Distance,ListDistance
USE MOD_Particle_Boundary_Condition, ONLY:PARTSWITCHELEMENT
USE MOD_Particle_MPI_Vars,           ONLY:PartHaloElemToProc
USE MOD_Mesh_Vars,                   ONLY:ElemToSide,BC,ElemBaryNGeo
USE MOD_Particle_MPI_Vars,           ONLY:SafetyFactor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
LOGICAL,INTENT(IN)                :: doHalo
LOGICAL,INTENT(IN)                :: initFix
LOGICAL,INTENT(IN)                :: doRelocate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                           :: InElementCheck,ParticleFound                                
REAL                              :: xi(1:3),Distance2,Det(6,2)
#ifdef MPI
INTEGER                           :: XiDir,locSideID,flip,SideID
REAL                              :: locXi,locEta,tmpXi
INTEGER                           :: Moved(2)
#endif /*MPI*/
!===================================================================================================================================
ParticleFound = .FALSE.
IF (KeepWallParticles) THEN
  IF (PDM%ParticleAtWall(iPart)) THEN
    PEM%Element(iPart) = PEM%lastElement(iPart)
    ParticleFound = .TRUE.
    RETURN 
  END IF
END IF

IF(DoHALO)THEN
  IF ( (PartState(iPart,1).LT.GEO%xminglob).OR.(PartState(iPart,1).GT.GEO%xmaxglob).OR. &
       (PartState(iPart,2).LT.GEO%yminglob).OR.(PartState(iPart,2).GT.GEO%ymaxglob).OR. &
       (PartState(iPart,3).LT.GEO%zminglob).OR.(PartState(iPart,3).GT.GEO%zmaxglob)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
ELSE
  IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
       (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
       (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
END IF

!IF (.NOT.DoRelocate) THEN
  !IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
       !(PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
       !(PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
     !PDM%ParticleInside(iPart) = .FALSE.
     !RETURN
  !END IF
!END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
CellX = MAX(MIN(GEO%TFIBGMimax,CellX),GEO%TFIBGMimin)
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MAX(MIN(GEO%TFIBGMjmax,CellY),GEO%TFIBGMjmin)
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MAX(MIN(GEO%TFIBGMkmax,CellZ),GEO%TFIBGMkmin)


IF (TriaTracking) THEN
  !--- check all cells associated with this background mesh cell
  DO iBGMElem = 1, GEO%FIBGM(CellX,CellY,CellZ)%nElem
    ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
    CALL ParticleInsideQuad3D(PartState(iPart,1:3),ElemID,InElementCheck,Det)
    IF (InElementCheck) THEN
       PEM%Element(iPart) = ElemID
       ParticleFound = .TRUE.
       EXIT
    END IF
  END DO
  IF (.NOT.ParticleFound) THEN
    PDM%ParticleInside(iPart) = .FALSE.
  END IF
  RETURN
END IF

!--- check all cells associated with this beckground mesh cell
nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem

! get closest element barycenter
Distance=-1.
ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  Distance2=(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
           +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
           +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) 
  IF(Distance2.GT.ElemRadius2NGeo(ElemID))THEN
    Distance(iBGMElem)=-1.
  ELSE
    Distance(iBGMElem)=Distance2
  END IF
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

IF(ALMOSTEQUAL(MAXVAL(Distance),-1.))THEN
  PDM%ParticleInside(iPart) = .FALSE.
  IF(DoRelocate) CALL abort(&
  __STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfoOpt=SafetyFactor)
  RETURN
END IF

!CALL BubbleSortID(Distance,ListDistance,nBGMElems)
IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element  
DO iBGMElem=1,nBGMElems
  IF(ALMOSTEQUAL(Distance(iBGMElem),-1.))CYCLE
  ElemID=ListDistance(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  IF(IsTracingBCElem(ElemID))THEN
    CALL PartInElemCheck(PartState(iPart,1:3),iPart,ElemID,InElementCheck)
    IF(.NOT.InElementCheck) CYCLE
  END IF

  CALL GetPositionInRefElem(PartState(iPart,1:3),xi,ElemID)
  IF(MAXVAL(ABS(Xi)).LT.epsOneCell(ElemID))THEN ! particle outside
    IF(.NOT.InitFix)THEN
      InElementCheck=.TRUE.
    ELSE
     InElementCheck=.TRUE.
     ! inelementcheck can only be set to false in the following part
#ifdef MPI
!     ! check if xi is larger than unity, than the
!     ! particle is found at least twice
     IF(MAXVAL(ABS(Xi)).GT.0.99999999)THEN ! particle possible outside
       XiDir = MAXLOC(ABS(Xi),1)
       ! now, get neighbor-side id
       SELECT CASE(XiDir)
       CASE(1) ! Xi
         IF(Xi(XiDir).GT.0)THEN
           ! XI_PLUS
           locSideID=XI_PLUS
           locXi=Xi(3)
           locEta=Xi(2)
         ELSE
           ! XI_MINUS
           locSideID=XI_MINUS
           locXi=Xi(2)
           locEta=Xi(3)
         END IF
       CASE(2) ! Eta
         IF(Xi(XiDir).GT.0)THEN
           locSideID=ETA_PLUS
           locXi=-Xi(1)
           locEta=Xi(3)
         ELSE
           locSideID=ETA_MINUS
           locXi=Xi(1)
           locEta=Xi(3)
         END IF
       CASE(3) ! Zeta
         IF(Xi(XiDir).GT.0)THEN
           locSideID=ZETA_PLUS
           locXi =Xi(1)
           locEta=Xi(2)
         ELSE
           locSideID=ZETA_MINUS
           locXi=Xi(2)
           locEta=Xi(1)
         END IF
       CASE DEFAULT
         CALL abort(&
__STAMP__&
, ' Error in  mesh-connectivity!')
       END SELECT
       ! get flip and rotate xi and eta into side-master system
       flip     =ElemToSide(E2S_FLIP,locSideID,ElemID)
       SideID   =ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
       SELECT CASE(Flip)
       CASE(1) ! slave side, SideID=q,jSide=p
         tmpXi=locEta
         locEta=locXi
         locXi=tmpXi
       CASE(2) ! slave side, SideID=N-p,jSide=q
         locXi=-locXi
         locEta=locEta
       CASE(3) ! slave side, SideID=N-q,jSide=N-p
         tmpXi =-locEta
         locEta=-locXi
         locXi=tmpXi
       CASE(4) ! slave side, SideID=p,jSide=N-q
         locXi =locXi
         locEta=-locEta
       END SELECT
       IF(BC(SideID).GT.0)THEN
         InElementCheck=.FALSE.
       ELSE
!         ! check if neighbor element is an mpi-element and if yes,
!         ! only take the particle if I am the lower rank
         Moved = PARTSWITCHELEMENT(locxi,loceta,locSideID,SideID,ElemID)
         IF(Moved(1).GT.PP_nElems)THEN
           IF(PartHaloElemToProc(NATIVE_PROC_ID,Moved(1)).LT.MyRank)THEN
             InElementCheck=.FALSE.
           END IF
         END IF
       END IF
     END IF
#endif /*MPI*/      
    END IF
  ELSE ! particle at face,edge or node, check most possible point
    InElementCheck=.FALSE.
  END IF
  IF (InElementCheck) THEN 
    PEM%Element(iPart) = ElemID
    IF(DoRefMapping) PartPosRef(1:3,iPart) = Xi
    ParticleFound = .TRUE.
    EXIT
  END IF
END DO ! iBGMElem

! particle not found
IF (.NOT.ParticleFound) THEN
  IF(DoRelocate) CALL abort(&
  __STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfoOpt=SafetyFactor)
  PDM%ParticleInside(iPart) = .FALSE.
END IF

END SUBROUTINE SingleParticleToExactElement


SUBROUTINE SingleParticleToExactElementNoMap(iPart,doHALO,doRelocate)
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,          ONLY:PartState,PEM,PDM
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:Geo
USE MOD_Utils,                  ONLY:InsertionSort !BubbleSortID
USE MOD_Particle_Tracking_Vars, ONLY:Distance,ListDistance
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
USE MOD_Particle_MPI_Vars,      ONLY:SafetyFactor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
LOGICAL,INTENT(IN)                :: doHalo
LOGICAL,INTENT(IN)                :: doRelocate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                           :: ParticleFound,InElementCheck
REAL                              :: Distance2
!===================================================================================================================================

ParticleFound = .FALSE.
IF(DoHALO)THEN
  IF ( (PartState(iPart,1).LT.GEO%xminglob).OR.(PartState(iPart,1).GT.GEO%xmaxglob).OR. &
       (PartState(iPart,2).LT.GEO%yminglob).OR.(PartState(iPart,2).GT.GEO%ymaxglob).OR. &
       (PartState(iPart,3).LT.GEO%zminglob).OR.(PartState(iPart,3).GT.GEO%zmaxglob)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
ELSE
  IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
       (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
       (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
CellX = MAX(MIN(GEO%TFIBGMimax,CellX),GEO%TFIBGMimin)
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MAX(MIN(GEO%TFIBGMjmax,CellY),GEO%TFIBGMjmin)
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MAX(MIN(GEO%TFIBGMkmax,CellZ),GEO%TFIBGMkmin)

!--- check all cells associated with this beckground mesh cell

nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem


! get closest element barycenter
Distance=-1.

ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  Distance2=(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
           +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
           +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) 
  IF(Distance2.GT.ElemRadius2NGeo(ElemID))THEN
    Distance(iBGMElem)=-1.
  ELSE
    Distance(iBGMElem)=Distance2
  END IF
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

IF(ALMOSTEQUAL(MAXVAL(Distance),-1.))THEN
  PDM%ParticleInside(iPart) = .FALSE.
  IF(DoRelocate)THEN
    IPWRITE(UNIT_StdOut,*) 'Position',PartState(iPart,1:3)
    CALL abort(&
  __STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfoOpt=SafetyFactor)
  END IF
  RETURN
END IF

!CALL BubbleSortID(Distance,ListDistance,nBGMElems)
IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)
! loop through sorted list and start by closest element  

DO iBGMElem=1,nBGMElems
  IF(ALMOSTEQUAL(Distance(iBGMElem),-1.))CYCLE
  ElemID=ListDistance(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  CALL PartInElemCheck(PartState(iPart,1:3),iPart,ElemID,InElementCheck)

  IF(InElementCheck)THEN
    ! no intersection found and particle is in final element
    PEM%Element(iPart) = ElemID
    ParticleFound=.TRUE.
    EXIT
  END IF
END DO ! iBGMElem

! particle not found
IF (.NOT.ParticleFound) THEN
!  IF(DoRelocate) CALL abort(&
!  __STAMP__&
!  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
!  RealInfoOpt=SafetyFactor)
  PDM%ParticleInside(iPart) = .FALSE.
END IF
END SUBROUTINE SingleParticleToExactElementNoMap


SUBROUTINE PartInElemCheck(PartPos_In,PartID,ElemID,FoundInElem,IntersectPoint_Opt,&
#ifdef CODE_ANALYZE
    Sanity_Opt,Tol_Opt,CodeAnalyze_Opt)
#else
       Tol_Opt)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
! Checks if particle is in Element
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
USE MOD_Particle_Surfaces_Vars, ONLY:SideType,SideNormVec
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartBCSideList
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Intersection,  ONLY:ComputePlanarRectIntersection
USE MOD_Particle_Intersection,  ONLY:ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection,  ONLY:ComputeBiLinearIntersection
USE MOD_Particle_Intersection,  ONLY:ComputeCurvedIntersection
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Globals,                ONLY:MyRank,UNIT_stdout
USE MOD_Mesh_Vars,              ONLY:NGeo
USE MOD_Particle_Tracking_Vars, ONLY:PartOut,MPIRankOut
USE MOD_Particle_Surfaces,      ONLY:OutputBezierControlPoints
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3d
USE MOD_Particle_Intersection,  ONLY:OutputTrajectory
#endif /*CODE_ANALYZE*/
USE MOD_Particle_Vars,          ONLY:LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: ElemID,PartID
REAL,INTENT(IN)                          :: PartPos_In(1:3)
#ifdef CODE_ANALYZE
LOGICAL,INTENT(IN),OPTIONAL              :: CodeAnalyze_Opt
LOGICAL,INTENT(IN),OPTIONAL              :: Sanity_Opt
#endif /*CODE_ANALYZE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                      :: FoundInElem
REAL,INTENT(OUT),OPTIONAL                :: IntersectPoint_Opt(1:3)
REAL,INTENT(OUT),OPTIONAL                :: Tol_Opt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: ilocSide,flip,SideID,BCSideID
REAL                                     :: PartTrajectory(1:3),NormVec(1:3)
REAL                                     :: lengthPartTrajectory,PartPos(1:3),LastPosTmp(1:3)
LOGICAL                                  :: isHit
REAL                                     :: alpha,eta,xi,IntersectPoint(1:3)
!===================================================================================================================================

IF(PRESENT(tol_Opt)) tol_Opt=-1.
! virtual move to element barycenter
LastPosTmp(1:3) =LastPartPos(PartID,1:3)
LastPartPos(PartID,1:3) =ElemBaryNGeo(1:3,ElemID)
PartPos(1:3) =PartPos_In(1:3)

PartTrajectory=PartPos - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )


#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      IPWRITE(UNIT_stdout,*) ' --------------------------------------------- '
      IPWRITE(UNIT_stdout,*) ' PartInElemCheck '
      CALL OutputTrajectory(PartID,PartPos,PartTrajectory,lengthPartTrajectory)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

IF(ALMOSTZERO(lengthPartTrajectory))THEN
  FoundInElem =.TRUE.
  LastPartPos(PartID,1:3) = LastPosTmp(1:3) 
  ! bugfix by Tilman
  RETURN
END IF
PartTrajectory=PartTrajectory/lengthPartTrajectory
isHit=.FALSE.
alpha=-1.
DO ilocSide=1,6

  !SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
  SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
  flip  = PartElemToSide(E2S_FLIP,ilocSide,ElemID)
  IF(DoRefMapping)THEN
    IF(SideID.LT.1) CYCLE
    BCSideID=SideID
    SideID=PartBCSideList(BCSideID)
    IF(SideID.LT.1) CYCLE
  END IF

  SELECT CASE(SideType(SideID))
  CASE(PLANAR_RECT)
    CALL ComputePlanarRectIntersection(ishit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta ,PartID,flip,SideID)
  CASE(PLANAR_CURVED)
    CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha,xi,eta,PartID,flip,SideID)
  CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha &
                                                                                       ,xi      &
                                                                                       ,eta      &
                                                                                       ,PartID,SideID &
                                                                                       ,ElemCheck_Opt=.TRUE.)
  CASE(CURVED)
    CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha,xi,eta,PartID,SideID,ElemCheck_Opt=.TRUE.)
  END SELECT

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(15("="))')
      WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (PartInElemCheck): '
      WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,'| Hit: ',isHit
      WRITE(UNIT_stdout,'(2(A,G0))')  '     | LengthPT: ',LengthPartTrajectory,' | Alpha: ',Alpha
      WRITE(UNIT_stdout,'(A,2(X,G0))') '     | Intersection xi/eta: ',xi,eta
    END IF
  END IF
  IF(PRESENT(Sanity_Opt))THEN
    IF(Sanity_Opt)THEN
      IF(alpha.GT.-1)THEN
        ! alpha is going from barycenter to point
        ! here, the tolerance for the ratio alpha/LengthPartTrajectory for tracing with element-corners is determined.
        IF(PRESENT(tol_Opt)) tol_Opt=MAX(ABS(1.-alpha/LengthPartTrajectory),tol_Opt)
        ! mark element as trouble element if rel. tol from alpha/LengthPartTrajectory to 1 > 1e-4
        ! tolerance 1e-4 is from ANSA_BOX grid (experimental, arbitrary)
        IF(ALMOSTEQUALRELATIVE(alpha/LengthPartTrajectory,1.,0.002)) THEN
          alpha=-1
        ELSE
          print*,'alpha',alpha,LengthPartTrajectory,alpha/LengthPartTrajectory,ABS(1.-alpha/LengthPartTrajectory),tol_Opt
        END IF
      END IF
    END IF
  END IF
  ! Dirty fix for PartInElemCheck if Lastpartpos is almost on side (tolerance issues) 
  IF(PRESENT(CodeAnalyze_Opt))THEN
    IF(CodeAnalyze_Opt)THEN
      IF((alpha)/LengthPartTrajectory.GT.0.9)THEN
        alpha = -1.0
      END IF
    END IF
  END IF
#endif /*CODE_ANALYZE*/
  IF(alpha.GT.-1)THEN
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      NormVec=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
    END SELECT 
    IF(flip.NE.0) NormVec=-NormVec
    IntersectPoint=LastPartPos(PartID,1:3)+alpha*PartTrajectory

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,*) '     | alpha          ',alpha
      WRITE(UNIT_stdout,*) '     | Normal vector  ',NormVec
      WRITE(UNIT_stdout,*) '     | PartTrajectory ',PartTrajectory
      WRITE(UNIT_stdout,*) '     | Dotprod        ',DOT_PRODUCT(NormVec,PartTrajectory)
      WRITE(UNIT_stdout,*) '     | Point 2        ', LastPartPos(PartID,1:3)+alpha*PartTrajectory+NormVec
      WRITE(UNIT_stdout,*) '     | Beziercontrolpoints3d-x'
      CALL OutputBezierControlPoints(BezierControlPoints3D_in=BezierControlPoints3D(1:3,:,:,SideID))
    END IF
  END IF
#endif /*CODE_ANALYZE*/

    IF(DOT_PRODUCT(NormVec,PartTrajectory).LT.0.)THEN
      alpha=-1.0
    ELSE
      EXIT
    END IF
    ! PO: should now be obsolete
    !IF(DoRefMapping)THEN
    !  IF(DOT_PRODUCT(NormVec,PartState(PartID,4:6)).LT.0.) alpha=-1.0
    !END IF ! DoRefMapping
  END IF
END DO ! ilocSide
FoundInElem=.TRUE.
IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt=0.
IF(alpha.GT.-1) THEN
  FoundInElem=.FALSE.
  IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt=IntersectPoint
END IF
LastPartPos(PartID,1:3) = LastPosTmp(1:3) 

END SUBROUTINE PartInElemCheck


SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det)
!===================================================================================================================================
! checks if particle is inside of linear element with triangulated faces
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,  ONLY : GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   ,INTENT(OUT)           :: Det(6,2)
LOGICAL,INTENT(OUT)           :: InElementCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, NodeNum
LOGICAL                       :: PosCheck, NegCheck
REAL                          :: A(1:3,1:4), cross(3)
!===================================================================================================================================
  InElementCheck = .TRUE.
  DO iLocSide = 1,6                 ! for all 6 sides of the element
     !--- initialize flags for side checks
     PosCheck = .FALSE.
     NegCheck = .FALSE.
     !--- A = vector from particle to node coords
     DO NodeNum = 1,4
       A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,ElemID)) - PartStateLoc(1:3)
     END DO

     !--- compute cross product for vector 1 and 3
     cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
     cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
     cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

     !--- negative determinant of triangle 1 (points 1,3,2):
     Det(iLocSide,1) = cross(1) * A(1,2) + &
                       cross(2) * A(2,2) + &
                       cross(3) * A(3,2)
     Det(iLocSide,1) = -det(iLocSide,1)
     !--- determinant of triangle 2 (points 1,3,4):
     Det(iLocSide,2) = cross(1) * A(1,4) + &
                       cross(2) * A(2,4) + &
                       cross(3) * A(3,4)
     IF (Det(iLocSide,1).LT.0) THEN
       NegCheck = .TRUE.
     ELSE
       PosCheck = .TRUE.
     END IF
     IF (Det(iLocSide,2).LT.0) THEN
       NegCheck = .TRUE.
     ELSE
       PosCheck = .TRUE.
     END IF

     !--- final determination whether particle is in element
     IF (GEO%ConcaveElemSide(iLocSide,ElemID)) THEN
       IF (.NOT.PosCheck) InElementCheck = .FALSE.
     ELSE
       IF (NegCheck) InElementCheck = .FALSE.
     END IF
  END DO
  RETURN
END SUBROUTINE ParticleInsideQuad3D


SUBROUTINE InitFIBGM()
!===================================================================================================================================
! Build Fast-Init-Background-Mesh.
! The BGM is a cartesian mesh for easier locating of particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GetRealArray,GetLogical
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,nTotalElems,nTotalBCSides, FindNeighbourElems
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
#ifdef MPI
USE MOD_Particle_MPI           ,ONLY: InitHALOMesh
USE MOD_Particle_MPI_Vars      ,ONLY: printMPINeighborWarnings,printBezierControlPointsWarnings
#endif /*MPI*/
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
USE MOD_PICDepo_Vars           ,ONLY: ElemRadius2_sf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: StartT,EndT
INTEGER                  :: iElem,ElemToBGM(1:6,1:PP_nElems)
INTEGER,ALLOCATABLE      :: HaloElemToBGM(:,:)
REAL,ALLOCATABLE         :: SideOrigin(:,:), SideRadius(:) 
!=================================================================================================================================

SWRITE(UNIT_StdOut,'(66("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT ELEMENT BASIS...' 
!! Read parameter for FastInitBackgroundMesh (FIBGM)
GEO%FIBGMdeltas(1:3) = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
GEO%FactorFIBGM(1:3) = GETREALARRAY('Part-FactorFIBGM',3,'1. , 1. , 1.')
GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

StartT=PICLASTIME()
ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:nTotalElems) &
        ,slenXiEtaZetaBasis(1:6,1:nTotalElems) &
        ,ElemRadiusNGeo(1:nTotalElems)         &
        ,ElemRadius2NGeo(1:nTotalElems)        )
CALL BuildElementBasis()
EndT=PICLASTIME()
IF(PartMPI%MPIROOT)THEN
  WRITE(UNIT_stdOut,'(A,F12.3,A)',ADVANCE='YES')' INIT ELEMENT-BASIS TOOK          [',EndT-StartT,'s]'
END IF

SWRITE(UNIT_StdOut,'(66("-"))')
StartT=PICLASTIME()
! get new min max
SWRITE(UNIT_stdOut,'(A)')' Getting FIBGM-minmax ...' 
CALL GetFIBGMminmax()
! sort elem in bgm cells
SWRITE(UNIT_stdOut,'(A)')' Getting element range in FIBGM ...' 
DO iElem=1,PP_nElems
  CALL BGMIndexOfElement(iElem,ElemToBGM(1:6,iElem)) 
END DO ! iElem = nElems+1,nTotalElems
SWRITE(UNIT_stdOut,'(A)')' Building FIBGM ...' 
CALL GetFIBGM(ElemToBGM)
EndT=PICLASTIME()
IF(PartMPI%MPIROOT)THEN
  WRITE(UNIT_stdOut,'(A,F12.3,A)',ADVANCE='YES')' Init FIBGM took                  [',EndT-StartT,'s]'
END IF
SWRITE(UNIT_StdOut,'(66("-"))')

CALL DuplicateSlavePeriodicSides()
! CAUTION: in MarkAllBCSides, a counter is reset for refmapping
CALL MarkAllBCSides()
! get elem and side types
CALL GetElemAndSideType()

StartT=PICLASTIME()
#ifdef MPI
SWRITE(UNIT_stdOut,'(A)')' INIT HALO REGION...' 
!CALL Initialize()  ! Initialize parallel environment for particle exchange between MPI domains
printMPINeighborWarnings = GETLOGICAL('printMPINeighborWarnings','.FALSE.')
printBezierControlPointsWarnings = GETLOGICAL('printBezierControlPointsWarnings','.FALSE.')
CALL InitHaloMesh()
! HALO mesh and region build. Unfortunately, the local FIBGM has to be extended to include the HALO elements :(
! rebuild is a local operation
#endif /*MPI*/

IF(nTotalElems.GT.PP_nElems)THEN
  ALLOCATE(HaloElemToBGM(1:6,PP_nElems+1:nTotalElems))
  DO iElem=PP_nElems+1,nTotalElems
    CALL BGMIndexOfElement(iElem,HaloElemToBGM(1:6,iElem)) 
  END DO ! iElem = nElems+1,nTotalElems
  CALL AddHALOCellsToFIBGM(ElemToBGM,HaloElemToBGM)
  DEALLOCATE(HaloElemToBGM)
ELSE
  CALL AddHALOCellsToFIBGM(ElemToBGM)
END IF

EndT=PICLASTIME()
IF(PartMPI%MPIROOT)THEN
   WRITE(UNIT_stdOut,'(A,F8.3,A)',ADVANCE='YES')' Construction of halo region took [',EndT-StartT,'s]'
END IF

IF(DoRefMapping)THEN
  IF(PartMPI%MPIROOT)THEN
     WRITE(UNIT_stdOut,'(A)') ' Reshaping arrays to reduced list...'
  END IF
  ! remove inner BezierControlPoints3D and SlabNormals, usw.
  CALL ReshapeBezierSides()
  ! compute side origin and radius for all sides in PartBCSideList
  IF(PartMPI%MPIROOT)THEN
     WRITE(UNIT_stdOut,'(A)') ' GetSideOrigin and Radius..'
  END IF
  ! remove inner BezierControlPoints3D and SlabNormals, usw.
  ALLOCATE( SideOrigin(1:3,1:nTotalBCSides) &
          , SideRadius(    1:nTotalBCSides) )
  CALL GetSideOriginAndRadius(nTotalBCSides,SideOrigin,SideRadius)
END IF

! get BCElem mapping, epsOnCell and calculate number of different elements and sides
IF(PartMPI%MPIROOT)THEN
   WRITE(UNIT_stdOut,'(A)') ' GetBCElemMap ...'
END IF
CALL GetBCElemMap()
IF(PartMPI%MPIROOT)THEN
   WRITE(UNIT_stdOut,'(A)') ' CaclElemAndSideNum ...'
END IF
CALL CalcElemAndSideNum()
! get basevectors for (bi-)linear sides
IF(PartMPI%MPIROOT)THEN
   WRITE(UNIT_stdOut,'(A)') ' LinearSideBaseVectors ...'
END IF
CALL GetLinearSideBaseVectors()
IF(PartMPI%MPIROOT)THEN
   WRITE(UNIT_stdOut,'(A)') ' Elem-Connectivity ...'
END IF
! check connectivity of particle mesh
CALL ElemConnectivity()

#if (PP_TimeDiscMethod==1001)
FindNeighbourElems = .TRUE.
#endif
IF (FindNeighbourElems) THEN
  ! build node conectivity of particle mesh
  IF(PartMPI%MPIROOT)THEN
     WRITE(UNIT_stdOut,'(A)') ' Node-Neighbourhood ...'
  END IF
  CALL NodeNeighbourhood()
END IF

SDEALLOCATE(XiEtaZetaBasis)
SDEALLOCATE(slenXiEtaZetaBasis)
SDEALLOCATE(ElemRadiusNGeo)
SDEALLOCATE(ElemRadius2NGeo)
ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:nTotalElems) &
        ,slenXiEtaZetaBasis(1:6,1:nTotalElems) &
        ,ElemRadiusNGeo(1:nTotalElems)         &
        ,ElemRadius2NGeo(1:nTotalElems)        )
SWRITE(UNIT_stdOut,'(A)')' BUILD ElementBasis ...'
SDEALLOCATE(ElemRadius2_sf) ! deallocate when using LB (it would be allocated twice because the call is executed twice)
CALL BuildElementBasis()
SWRITE(UNIT_stdOut,'(A)')' BUILD ElementBasis DONE!'
IF(DoRefMapping) THEN
  ! compute distance between each side associated with  the element and its origin
  CALL GetElemToSideDistance(nTotalBCSides,SideOrigin,SideRadius)
  DEALLOCATE( SideOrigin, SideRadius)
END IF
SWRITE(UNIT_stdOut,'(A)')' ... DONE!' 
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitFIBGM


SUBROUTINE GetFIBGM(ElemToBGM)
!===================================================================================================================================
! build local FIBGM mesh for process local FIBGM mesh including HALO region
! mode 1: build local BGM and interconnections with other processes
! mode 2: rebuild BGM including HALO region
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Partilce_Periodic_BC ,ONLY: InitPeriodicBC
USE MOD_Particle_Mesh_Vars   ,ONLY: GEO
USE MOD_PICDepo              ,ONLY: InitializeDeposition
USE MOD_Particle_MPI_Vars    ,ONLY: SafetyFactor,halo_eps_velo,halo_eps,halo_eps2
#ifndef PP_HDG
USE MOD_CalcTimeStep         ,ONLY: CalcTimeStep
#endif /*PP_HDG*/
USE MOD_Equation_Vars        ,ONLY: c
USE MOD_Particle_Vars        ,ONLY: manualtimestep
#if (PP_TimeDiscMethod==201)
USE MOD_Particle_Vars        ,ONLY: dt_part_ratio
#endif
USE MOD_ChangeBasis          ,ONLY: ChangeBasis2D
#ifdef MPI
USE MOD_Particle_MPI         ,ONLY: InitHALOMesh
USE MOD_Particle_Mesh_Vars   ,ONLY: FIBGMCellPadding
USE MOD_PICDepo_Vars         ,ONLY: DepositionType, r_sf
USE MOD_Particle_MPI_Vars    ,ONLY: PartMPI
USE MOD_Particle_Mesh_Vars   ,ONLY: NbrOfCases,casematrix
#endif /*MPI*/
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars        ,ONLY: RK_c,nRKStages
#endif
USE MOD_ReadInTools          ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN)    :: mode
INTEGER,INTENT(IN)     :: ElemToBGM(1:6,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                  :: localXmin,localXmax,localymin,localymax,localzmin,localzmax
INTEGER                          :: BGMimin,BGMimax,BGMjmin,BGMjmax,BGMkmin,BGMkmax
!REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                          :: iBGM,jBGM,kBGM,iElem
INTEGER                          :: BGMCellXmax,BGMCellXmin
INTEGER                          :: BGMCellYmax,BGMCellYmin
INTEGER                          :: BGMCellZmax,BGMCellZmin
INTEGER                          :: ALLOCSTAT
INTEGER                          :: iProc
REAL                             :: deltaT
REAL                             :: globalDiag
#ifdef MPI
INTEGER                          :: ii,jj,kk,i,j
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
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
INTEGER                          :: iStage
#endif
!===================================================================================================================================

! zeros
#ifdef MPI
ii=0
jj=0
kk=0
#endif /*MPI*/


#ifdef MPI
  ! allocate and initialize MPINeighbor
  ALLOCATE(PartMPI%isMPINeighbor(0:PartMPI%nProcs-1))
  PartMPI%isMPINeighbor(:) = .FALSE.
  PartMPI%nMPINeighbors=0
#endif   


  CALL InitPeriodicBC()
  ! reduce beziercontrolpoints to boundary sides
  !IF(DoRefMapping) CALL ReshapeBezierSides()
  !CALL InitializeInterpolation() ! not any more required ! has to be called earliear
  CALL InitializeDeposition()     ! has to remain here, because domain can have changed
  !CALL InitPIC()                 ! does not depend on domain

! deallocate stuff // required for dynamic load balance
#ifdef MPI
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
#endif /*MPI*/

!--- compute number of background cells in each direction
!BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)+1.00001)
!BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)+0.99999)
!BGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)+1.00001)
!BGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)+0.99999)
!BGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)+1.00001)
!BGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)+0.99999)

! now fail safe, enlarge the BGM grid for safety reasons
BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))-1
BGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
BGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))-1
BGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
BGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))-1

!--- JN: For MPI communication, information also about the neighboring FIBGM cells is needed
!--- AS: shouldn't we add up here the nPaddingCells? 
!--- TS: What we need to do is increase the BGM area for shape_function ONLY
!        Reason: if a particle moves outside the domain, there still needs to be a
!                BGM with an associated ShapeProc at the particle position
!        Particle may only move c*dt*Safetyfactor.
!--- PO: modified for curved and shape-function influence
!        c*dt*SafetyFactor+r_cutoff
IF (ManualTimeStep.EQ.0.0) THEN
#ifndef PP_HDG
  deltaT=CALCTIMESTEP()
#else
   CALL abort(&
__STAMP__&
, 'ManualTimeStep.EQ.0.0 -> ManualTimeStep is not defined correctly! Particles-ManualTimeStep = ',RealInfoOpt=ManualTimeStep)
#endif /*PP_HDG*/
ELSE
  deltaT=ManualTimeStep
END IF
IF (halo_eps_velo.EQ.0) halo_eps_velo = c
#if (PP_TimeDiscMethod==4 || PP_TimeDiscMethod==200 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==1000)
IF (halo_eps_velo.EQ.c) THEN
   CALL abort(&
__STAMP__&
, 'halo_eps_velo.EQ.c -> Halo Eps Velocity for MPI not defined')
END IF
#endif
#if (PP_TimeDiscMethod==201)
deltaT=CALCTIMESTEP()
halo_eps = c*deltaT*SafetyFactor*max(dt_part_ratio,1.0)
#elif (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
halo_eps = RK_c(2)
DO iStage=2,nRKStages-1
  halo_eps = MAX(halo_eps,RK_c(iStage+1)-RK_c(iStage))
END DO
halo_eps = MAX(halo_eps,1.-RK_c(nRKStages))
SWRITE(UNIT_stdOut,'(A38,E24.12)') ' |                 max. RKdtFrac  |    ',halo_eps 
halo_eps = halo_eps*halo_eps_velo*deltaT*SafetyFactor !dt multiplied with maximum RKdtFrac
#else
halo_eps = halo_eps_velo*deltaT*SafetyFactor ! for RK too large
#endif

! limit halo_eps to diagonal of bounding box
globalDiag = SQRT( (GEO%xmaxglob-GEO%xminglob)**2 & 
                 + (GEO%ymaxglob-GEO%yminglob)**2 & 
                 + (GEO%zmaxglob-GEO%zminglob)**2 ) 
IF(halo_eps.GT.globalDiag)THEN
  CALL PrintOption('unlimited halo distance','CALCUL.',RealOpt=halo_eps)
  SWRITE(UNIT_stdOut,'(A38)') ' |   limitation of halo distance  |    '
  halo_eps=globalDiag
END IF

halo_eps2=halo_eps*halo_eps
CALL PrintOption('halo distance','CALCUL.',RealOpt=halo_eps)


#ifdef MPI
IF ((DepositionType.EQ.'shape_function')             &
.OR.(DepositionType.EQ.'shape_function_cylindrical') &
.OR.(DepositionType.EQ.'shape_function_spherical')   &
.OR.(DepositionType.EQ.'shape_function_simple')      &
.OR.(DepositionType.EQ.'shape_function_1d')          )THEN
  ! and changed, tooo
  BGMimax = INT((MIN(GEO%xmax+halo_eps,GEO%xmaxglob)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
  BGMimin = INT((MAX(GEO%xmin-halo_eps,GEO%xminglob)-GEO%xminglob)/GEO%FIBGMdeltas(1))-1
  BGMjmax = INT((MIN(GEO%ymax+halo_eps,GEO%ymaxglob)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
  BGMjmin = INT((MAX(GEO%ymin-halo_eps,GEO%yminglob)-GEO%yminglob)/GEO%FIBGMdeltas(2))-1
  BGMkmax = INT((MIN(GEO%zmax+halo_eps,GEO%zmaxglob)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
  BGMkmin = INT((MAX(GEO%zmin-halo_eps,GEO%zminglob)-GEO%zminglob)/GEO%FIBGMdeltas(3))-1
END IF
#endif

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
  CALL abort(&
__STAMP__&
, 'Problem allocating GEO%FIBGM!' )
END IF

! null number of element per BGM cell
DO kBGM = BGMkmin,BGMkmax
   DO jBGM = BGMjmin,BGMjmax
     DO iBGM = BGMimin,BGMimax
         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
      END DO
   END DO
END DO

!--- compute number of elements in each background cell
DO iElem=1,PP_nElems
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = ElemToBGM(1,iElem)
  BGMCellXmax = ElemToBGM(2,iElem)
  BGMCellYmin = ElemToBGM(3,iElem)
  BGMCellYmax = ElemToBGM(4,iElem)
  BGMCellZmin = ElemToBGM(5,iElem)
  BGMCellZmax = ElemToBGM(6,iElem)
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
DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
      IF(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem.EQ.0) CYCLE
      ALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%FIBGM(iBGM,jBGM,kBGM)%nElem))
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!--- map elements to background cells
DO iElem=1,PP_nElems
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = ElemToBGM(1,iElem)
  BGMCellXmax = ElemToBGM(2,iElem)
  BGMCellYmin = ElemToBGM(3,iElem)
  BGMCellYmax = ElemToBGM(4,iElem)
  BGMCellZmin = ElemToBGM(5,iElem)
  BGMCellZmax = ElemToBGM(6,iElem)
  ! add current Element to BGM-Elem
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1    
        GEO%FIBGM(iBGM,jBGM,kBGM)%Element(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem


!IF(mode.EQ.2) RETURN
SWRITE(UNIT_stdOut,'(A)')' Building MPI-FIBGM ...' 
#ifdef MPI
!--- MPI stuff for background mesh (FastinitBGM)
BGMCells=0 
ALLOCATE(BGMCellsArray(1:(BGMimax-BGMimin+1)*(BGMjmax-BGMjmin+1)*(BGMkmax-BGMkmin+1)*3))
!Count BGMCells with Elements inside and save their indices in BGMCellsArray
DO kBGM=BGMkmin, BGMkmax
  DO jBGM=BGMjmin, BGMjmax
    DO iBGM=BGMimin, BGMimax  
      IF (GEO%FIBGM(iBGM,jBGM,kBGM)%nElem .GT. 0) THEN
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
!Gather indices of every Procs' Cells
CALL MPI_ALLGATHERV(BGMCellsArray(1:BGMCells*3), BGMCells*3, MPI_INTEGER, GlobalBGMCellsArray, &    
                   & NbrOfBGMCells(0:PartMPI%nProcs-1)*3, Displacement, MPI_INTEGER, PartMPI%COMM, IERROR)

!--- JN: first: count required array size for ReducedBGMArray
!--- TS: Define padding stencil (max of halo and shape padding)
!        Reason: This padding is used to build the ReducedBGM, so any information 
!                outside this region is lost 
IF (GEO%nPeriodicVectors.GT.0) THEN  !Periodic (can't be done below because ReducedBGMArray is sorted by proc)
  FIBGMCellPadding(1:3)=1
  IF(.NOT.GEO%directions(1)) FIBGMCellPadding(1) = INT(halo_eps/GEO%FIBGMdeltas(1))+1
  IF(.NOT.GEO%directions(2)) FIBGMCellPadding(2) = INT(halo_eps/GEO%FIBGMdeltas(2))+1
  IF(.NOT.GEO%directions(3)) FIBGMCellPadding(3) = INT(halo_eps/GEO%FIBGMdeltas(3))+1
ELSE
  FIBGMCellPadding(1:3) = INT(halo_eps/GEO%FIBGMdeltas(1:3))+1
END IF
! halo region already included in BGM
!FIBGMCellPadding(1:3) = 0
nShapePaddingX = 0
nShapePaddingY = 0
nShapePaddingZ = 0
IF ((DepositionType.EQ.'shape_function')             &
.OR.(DepositionType.EQ.'shape_function_cylindrical') &
.OR.(DepositionType.EQ.'shape_function_simple')      &
.OR.(DepositionType.EQ.'shape_function_spherical')   &
.OR.(DepositionType.EQ.'shape_function_1d')          )THEN
  nShapePaddingX = INT(r_sf/GEO%FIBGMdeltas(1)+0.9999999)
  nShapePaddingY = INT(r_sf/GEO%FIBGMdeltas(2)+0.9999999)
  nShapePaddingZ = INT(r_sf/GEO%FIBGMdeltas(3)+0.9999999)
  !IPWRITE(*,*) 'nShapePaddingX',nShapePaddingX
  !IPWRITE(*,*) 'nShapePaddingY',nShapePaddingY
  !IPWRITE(*,*) 'nShapePaddingZ',nShapePaddingZ
 ! IF(mode.EQ.2) THEN
 !   IF((nShapePaddingX.EQ.0)    &
 !     .OR.(nShapePaddingY.EQ.0) &
 !     .OR.(nShapePaddingZ.EQ.0))THEN 
 !       CALL abort(__STAMP__&
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

! now create a temporary array in which for all BGM Cells + ShapePadding the processes are saved 
! reason: this way, the ReducedBGM List only needs to be searched once and not once for each BGM Cell+Stencil

! first count the maximum number of procs that exist within each BGM cell (inkl. Shape Padding region)
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
          !  IPWRITE(UNIT_stdOut,*) ' Warning, something wrong with halo region'
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
!           CALL abort(__STAMP__&
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



SUBROUTINE AddHALOCellsToFIBGM(ElemToBGM,HaloElemToBGM)
!===================================================================================================================================
! remap all elements including halo-elements into FIBGM
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,            ONLY : UNIT_StdOut
USE MOD_ChangeBasis,                        ONLY:ChangeBasis2D
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO,nTotalElems
USE MOD_Particle_MPI_Vars,                  ONLY:PartMPI
USE MOD_Particle_Tracking_Vars,             ONLY:Distance,ListDistance
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN)    :: mode
INTEGER,INTENT(IN)               :: ElemToBGM(1:6,1:PP_nElems)
INTEGER,INTENT(IN),OPTIONAL      :: HaloElemToBGM(1:6,PP_nElems+1:nTotalElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: BGMimin,BGMimax,BGMjmin,BGMjmax,BGMkmin,BGMkmax,Allocstat
REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                          :: iBGM,jBGM,kBGM,iElem
INTEGER                          :: BGMCellXmax,BGMCellXmin
INTEGER                          :: BGMCellYmax,BGMCellYmin
INTEGER                          :: BGMCellZmax,BGMCellZmin
LOGICAL, ALLOCATABLE             :: ElementFound(:)
INTEGER                          :: maxnBGMElems
!===================================================================================================================================


! current min,max
BGMimax=GEO%FIBGMimax
BGMimin=GEO%FIBGMimin
BGMjmax=GEO%FIBGMjmax
BGMjmin=GEO%FIBGMjmin
BGMkmax=GEO%FIBGMkmax
BGMkmin=GEO%FIBGMkmin

GEO%TFIBGMimax =GEO%FIBGMimax
GEO%TFIBGMimin =GEO%FIBGMimin
GEO%TFIBGMjmax =GEO%FIBGMjmax
GEO%TFIBGMjmin =GEO%FIBGMjmin
GEO%TFIBGMkmax =GEO%FIBGMkmax
GEO%TFIBGMkmin =GEO%FIBGMkmin

BGMCellXmax = BGMimax
BGMCellXmin = BGMimin
BGMCellYmax = BGMjmax
BGMCellYmin = BGMjmin
BGMCellZmax = BGMkmax
BGMCellZmin = BGMkmin


DO iElem=1,nTotalElems
  IF(iElem.LE.PP_nElems)THEN
    BGMCellXmin = ElemToBGM(1,iElem)
    BGMCellXmax = ElemToBGM(2,iElem)
    BGMCellYmin = ElemToBGM(3,iElem)
    BGMCellYmax = ElemToBGM(4,iElem)
    BGMCellZmin = ElemToBGM(5,iElem)
    BGMCellZmax = ElemToBGM(6,iElem)
  ELSE
    IF(.NOT.GEO%directions(1)) BGMCellXmin = HaloElemToBGM(1,iElem)
    IF(.NOT.GEO%directions(1)) BGMCellXmax = HaloElemToBGM(2,iElem)
    IF(.NOT.GEO%directions(2)) BGMCellYmin = HaloElemToBGM(3,iElem)
    IF(.NOT.GEO%directions(2)) BGMCellYmax = HaloElemToBGM(4,iElem)
    IF(.NOT.GEO%directions(3)) BGMCellZmin = HaloElemToBGM(5,iElem)
    IF(.NOT.GEO%directions(3)) BGMCellZmax = HaloElemToBGM(6,iElem)
  END IF

  BGMimin=MIN(BGMimin,BGMCellXmin)
  BGMimax=MAX(BGMimax,BGMCellXmax)
  BGMjmin=MIN(BGMjmin,BGMCellYmin)
  BGMjmax=MAX(BGMjmax,BGMCellYmax)
  BGMkmin=MIN(BGMkmin,BGMCellZmin)
  BGMkmax=MAX(BGMkmax,BGMCellZmax)

END DO ! iElem = nElems+1,nTotalElems

GEO%TFIBGMimax =BGMimax
GEO%TFIBGMimin =BGMimin
GEO%TFIBGMjmax =BGMjmax
GEO%TFIBGMjmin =BGMjmin
GEO%TFIBGMkmax =BGMkmax
GEO%TFIBGMkmin =BGMkmin

ALLOCATE(GEO%TFIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
__STAMP__&
,' ERROR in AddElemsToTFIBGM: Cannot allocate GEO%TFIBGM!')
END IF

ALLOCATE( ElementFound(1:nTotalElems) )
ElementFound = .FALSE.

! null number of elements per BGM-Cell
DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
       GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM


!--- compute number of elements in each background cell
DO iElem=1,PP_nElems
  !--- find minimum and maximum BGM cell for current element
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = MIN(MAX(ElemToBGM(1,iElem),BGMimin),BGMimax)
  BGMCellXmax = MAX(MIN(ElemToBGM(2,iElem),BGMimax),BGMimin)
  BGMCellYmin = MIN(MAX(ElemToBGM(3,iElem),BGMjmin),BGMjmax)
  BGMCellYmax = MAX(MIN(ElemToBGM(4,iElem),BGMjmax),BGMjmin)
  BGMCellZmin = MIN(MAX(ElemToBGM(5,iElem),BGMkmin),BGMkmax)
  BGMCellZmax = MAX(MIN(ElemToBGM(6,iElem),BGMkmax),BGMkmin)
  ! add ecurrent element to number of BGM-elems
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
         GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1
         ElementFound(iElem) = .TRUE.
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

DO iElem=PP_nElems+1,nTotalElems
  !--- find minimum and maximum BGM cell for current element
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = MIN(MAX(HaloElemToBGM(1,iElem),BGMimin),BGMimax)
  BGMCellXmax = MAX(MIN(HaloElemToBGM(2,iElem),BGMimax),BGMimin)
  BGMCellYmin = MIN(MAX(HaloElemToBGM(3,iElem),BGMjmin),BGMjmax)
  BGMCellYmax = MAX(MIN(HaloElemToBGM(4,iElem),BGMjmax),BGMjmin)
  BGMCellZmin = MIN(MAX(HaloElemToBGM(5,iElem),BGMkmin),BGMkmax)
  BGMCellZmax = MAX(MIN(HaloElemToBGM(6,iElem),BGMkmax),BGMkmin)
  ! add ecurrent element to number of BGM-elems
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1
        ElementFound(iElem) = .TRUE.
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem


!--- allocate mapping variable and clean number for mapping (below)
DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
      IF(GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem.EQ.0) CYCLE
      ALLOCATE(GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem))
      GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!--- map elements to background cells
DO iElem=1,PP_nElems
  !--- find minimum and maximum BGM cell for current element
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = MIN(MAX(ElemToBGM(1,iElem),BGMimin),BGMimax)
  BGMCellXmax = MAX(MIN(ElemToBGM(2,iElem),BGMimax),BGMimin)
  BGMCellYmin = MIN(MAX(ElemToBGM(3,iElem),BGMjmin),BGMjmax)
  BGMCellYmax = MAX(MIN(ElemToBGM(4,iElem),BGMjmax),BGMjmin)
  BGMCellZmin = MIN(MAX(ElemToBGM(5,iElem),BGMkmin),BGMkmax)
  BGMCellZmax = MAX(MIN(ElemToBGM(6,iElem),BGMkmax),BGMkmin)

  ! add current Element to BGM-Elem
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
        GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1    
        GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem
DO iElem=PP_nElems+1,nTotalElems
  !--- find minimum and maximum BGM cell for current element
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = MIN(MAX(HaloElemToBGM(1,iElem),BGMimin),BGMimax)
  BGMCellXmax = MAX(MIN(HaloElemToBGM(2,iElem),BGMimax),BGMimin)
  BGMCellYmin = MIN(MAX(HaloElemToBGM(3,iElem),BGMjmin),BGMjmax)
  BGMCellYmax = MAX(MIN(HaloElemToBGM(4,iElem),BGMjmax),BGMjmin)
  BGMCellZmin = MIN(MAX(HaloElemToBGM(5,iElem),BGMkmin),BGMkmax)
  BGMCellZmax = MAX(MIN(HaloElemToBGM(6,iElem),BGMkmax),BGMkmin)

  ! add current Element to BGM-Elem
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
        GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem = GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem + 1    
        GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem


DO iElem=1,PP_nElems
  IF(.NOT.ElementFound(iElem))THEN
    !--- find minimum and maximum BGM cell for current element
    ! here fancy stuff, because element could be wide out of element range
    BGMCellXmin = ElemToBGM(1,iElem)
    BGMCellXmax = ElemToBGM(2,iElem)
    BGMCellYmin = ElemToBGM(3,iElem)
    BGMCellYmax = ElemToBGM(4,iElem)
    BGMCellZmin = ElemToBGM(5,iElem)
    BGMCellZmax = ElemToBGM(6,iElem)

    IPWRITE(UNIT_stdOut,*) ' TFIBGM , iElem'
    IPWRITE(UNIT_stdOut,*) 'xmin',GEO%xmin
    IPWRITE(UNIT_stdOut,*) 'xmax',GEO%xmax
    IPWRITE(UNIT_stdOut,*) 'ymin',GEO%ymin
    IPWRITE(UNIT_stdOut,*) 'ymax',GEO%ymax
    IPWRITE(UNIT_stdOut,*) 'zmin',GEO%zmin
    IPWRITE(UNIT_stdOut,*) 'zmax',GEO%zmax
    IPWRITE(UNIT_stdOut,*) ' BGM , iBGM'
    IPWRITE(UNIT_stdOut,*) 'xmin', BGMimin,BGMCellXmin
    IPWRITE(UNIT_stdOut,*) 'xmax', BGMimax,BGMCellXmax
    IPWRITE(UNIT_stdOut,*) 'ymin', BGMjmin,BGMCellYmin
    IPWRITE(UNIT_stdOut,*) 'ymax', BGMjmax,BGMCellYmax
    IPWRITE(UNIT_stdOut,*) 'zmin', BGMkmin,BGMCellYmin
    IPWRITE(UNIT_stdOut,*) 'zmax', BGMkmax,BGMCellYmax
    CALL abort(&
__STAMP__&
,' Element not located in FIBGM! iElem, myRank',iElem,REAL(PartMPI%MyRank))
  END IF
END DO ! iElem

DO iElem=PP_nElems+1,nTotalElems
  IF(.NOT.ElementFound(iElem))THEN
    !--- find minimum and maximum BGM cell for current element
    ! here fancy stuff, because element could be wide out of element range
    BGMCellXmin = HaloElemToBGM(1,iElem)
    BGMCellXmax = HaloElemToBGM(2,iElem)
    BGMCellYmin = HaloElemToBGM(3,iElem)
    BGMCellYmax = HaloElemToBGM(4,iElem)
    BGMCellZmin = HaloElemToBGM(5,iElem)
    BGMCellZmax = HaloElemToBGM(6,iElem)

    IPWRITE(UNIT_stdOut,*) ' TFIBGM , iElem'
    IPWRITE(UNIT_stdOut,*) 'xmin',GEO%xmin,xmin
    IPWRITE(UNIT_stdOut,*) 'xmax',GEO%xmax,xmax
    IPWRITE(UNIT_stdOut,*) 'ymin',GEO%ymin,ymin
    IPWRITE(UNIT_stdOut,*) 'ymax',GEO%ymax,ymax
    IPWRITE(UNIT_stdOut,*) 'zmin',GEO%zmin,zmin
    IPWRITE(UNIT_stdOut,*) 'zmax',GEO%zmax,zmax
    IPWRITE(UNIT_stdOut,*) ' BGM , iBGM'
    IPWRITE(UNIT_stdOut,*) 'xmin', BGMimin,BGMCellXmin
    IPWRITE(UNIT_stdOut,*) 'xmax', BGMimax,BGMCellXmax
    IPWRITE(UNIT_stdOut,*) 'ymin', BGMjmin,BGMCellYmin
    IPWRITE(UNIT_stdOut,*) 'ymax', BGMjmax,BGMCellYmax
    IPWRITE(UNIT_stdOut,*) 'zmin', BGMkmin,BGMCellYmin
    IPWRITE(UNIT_stdOut,*) 'zmax', BGMkmax,BGMCellYmax
    CALL abort(&
__STAMP__&
,' Element not located in FIBGM! iElem, myRank',iElem,REAL(PartMPI%MyRank))
  END IF
END DO ! iElem


DEALLOCATE(Elementfound)

! and get max number of bgm-elems
maxnBGMElems=0
DO kBGM = GEO%TFIBGMkmin,GEO%TFIBGMkmax
  DO jBGM = GEO%TFIBGMjmin,GEO%TFIBGMjmax
    DO iBGM = GEO%TFIBGMimin,GEO%TFIBGMimax
      !maxnBGMElems=MAX(maxnBGMElems,GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem)
      maxnBGMElems=MAX(maxnBGMElems,GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem)
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM
ALLOCATE(Distance    (1:maxnBGMElems) &
        ,ListDistance(1:maxnBGMElems) )


END SUBROUTINE AddHALOCellsToFIBGM


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
USE MOD_Mesh_Vars,          ONLY : nElems
USE MOD_Particle_Mesh_Vars, ONLY : GEO, WeirdElems
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
IF(WeirdElems.GT.0) THEN
  IPWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
  IPWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
  DO iElem = 1,WeirdElems
    IPWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
  !  DO iLocSide = 1,6
  !    DO iNode = 1,4
  !      SWRITE(UNIT_stdOut,*) GEO%NodeCoords(:,GEO%ElemSideNodeID(iNode,iLocSide,WeirdElemNbrs(iElem)))
  !    END DO
  !  END DO
  !  STOP
  END DO
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WeirdElementCheck



SUBROUTINE ReShapeBezierSides()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalBCSides,PartBCSideList,nTotalSides,nPartPeriodicSides
USE MOD_Mesh_Vars,              ONLY:nSides,nBCSides,NGeo,BC
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,SideType,SideDistance,SideNormVec
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
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
INTEGER           :: iSide,nOldBCSides,newBCSideID,BCInc,nPeriodicSidesTmp
REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummySideSlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummySideSlabIntervals               ! intervalls beta1, beta2, beta3
LOGICAL,ALLOCATABLE,DIMENSION(:)   :: DummyBoundingBoxIsEmpty
REAL,ALLOCATABLE                   :: DummyBezierControlPoints3D(:,:,:,:)                                
INTEGER,ALLOCATABLE                :: DummySideType(:)
REAL,ALLOCATABLE                   :: DummySideDistance(:)
REAL,ALLOCATABLE                   :: DummySideNormVec(:,:)
!===================================================================================================================================

nPeriodicSidesTmp=0
DO iSide=nBCSides+1,nSides+nPartPeriodicSides
  IF(BC(iSide).NE.0)THEN
    ! different list, contains ALL periodic sides (inner and duplicated)
    nPeriodicSidesTmp=nPeriodicSidesTmp+1
  END IF
END DO

! now, shrink partbcsidelist
nOldBCSides  =nTotalBCSides
nTotalBCSides=nTotalSides-nPartPeriodicSides-nSides+nBCSides+nPeriodicSidesTmp

IF(nTotalBCSides.EQ.0) RETURN

! allocate & fill dummy
! BezierControlPoints3D
ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides))
IF (.NOT.ALLOCATED(DummyBezierControlPoints3d)) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate DummyBezierControlPoints3D in ReshapeBezierSides')
IF (SIZE(DummyBezierControlPoints3D).NE.SIZE(BezierControlPoints3D)) CALL abort(&
__STAMP__&
,'size of DummyBezierControlPoionts3D and BezierControlPoints3D not equal!')
DummyBezierControlPoints3d=BezierControlPoints3d
DEALLOCATE(BezierControlPoints3D)
ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalBCSides),STAT=ALLOCSTAT)
BezierControlPoints3d=0.
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate BezierControlPoints3D in ReshapeBezierSides')
! SideSlabNormals
ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nTotalSides))
IF (.NOT.ALLOCATED(DummySideSlabNormals)) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate DummySideSlabNormals in ReshapeBezierSides')
DummySideSlabNormals=SideSlabNormals
DEALLOCATE(SideSlabNormals)
ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate SideSlabNormals in ReshapeBezierSides')
SideSlabNormals=0.
! SideSlabIntervals
ALLOCATE(DummySideSlabIntervals(1:6,1:nTotalSides))
IF (.NOT.ALLOCATED(DummySideSlabIntervals)) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate DummySideSlabIntervals in ReshapeBezierSides')
DummySideSlabIntervals=SideSlabIntervals
DEALLOCATE(SideSlabIntervals)
ALLOCATE(SideSlabIntervals(1:6,1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate ElemIndex in ReshapeBezierSides')
SideSlabIntervals=0.
! BoundingBoxIsEmpty
ALLOCATE(DummyBoundingBoxIsEmpty(1:nTotalSides))
IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate DummyBoundingBoxIsEmpty in ReshapeBezierSides')
DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
DEALLOCATE(BoundingBoxIsEmpty)
ALLOCATE(BoundingBoxIsEmpty(1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__& !wunderschoen!!!
,'Could not allocate BoundingBoxIsEmpty in ReshapeBezierSides')
BoundingBoxIsEmpty=.FALSE.
! side type
ALLOCATE(DummySideType(1:nOldBCSides))
IF (.NOT.ALLOCATED(DummySideType)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummySideType in ReshapeBezierSides')
DummySideType=SideType
DEALLOCATE(SideType)
ALLOCATE(SideType(1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate SideType in ReshapeBezierSides')
SideType=-1

ALLOCATE(DummySideDistance(1:nOldBCSides))
IF (.NOT.ALLOCATED(DummySideDistance)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummySideDistance in ReshapeBezierSides')
DummySideDistance=SideDistance
DEALLOCATE(SideDistance)
ALLOCATE(SideDistance(1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate SideDistance in ReshapeBezierSides')
SideDistance=0.

ALLOCATE(DummySideNormVec(1:3,1:nOldBCSides))
IF (.NOT.ALLOCATED(DummySideNormVec)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummySideNormVec in ReshapeBezierSides')
DummySideNormVec=SideNormVec
DEALLOCATE(SideNormVec)
ALLOCATE(SideNormVec(1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate SideNormVec in ReshapeBezierSides')
SideNormVec=0.


BCInc=0
!DO iSide=1,nSides
newBCSideID=0
DO iSide=1,nBCSides
  newBCSideID=newBCSideID+1
  BezierControlPoints3d(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
  SideSlabNormals          (1:3,1:3,          newBCSideID) =DummySideSlabNormals         (1:3,1:3,           iSide)
  SideSlabIntervals       (1:6,              newBCSideID) =DummySideSlabIntervals      (1:6,               iSide)
  BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
  SideType(newBCSideID) = DummySideType(newBCSideID)
  SideDistance(newBCSideID) = DummySideDistance(newBCSideID)
  SideNormVec(1:3,newBCSideID) = DummySideNormVec(1:3,newBCSideID)
END DO ! iSide

DO iSide=nBCSides+1,nSides+nPartPeriodicSides
  IF(BC(iSide).EQ.0) CYCLE
  newBCSideID=newBCSideID+1
  BezierControlPoints3d(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
  SideSlabNormals          (1:3,1:3,          newBCSideID) =DummySideSlabNormals         (1:3,1:3,           iSide)
  SideSlabIntervals       (1:6,              newBCSideID) =DummySideSlabIntervals      (1:6,               iSide)
  BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
  SideType(newBCSideID) = DummySideType(newBCSideID)
  SideDistance(newBCSideID) = DummySideDistance(newBCSideID)
  SideNormVec(1:3,newBCSideID) = DummySideNormVec(1:3,newBCSideID)
END DO ! iSide

DO iSide=nSides+nPartPeriodicSides+1,nTotalSides
  newBCSideID=newBCSideID+1
  BezierControlPoints3d(1:3,0:NGeo,0:NGeo,newBCSideID) =DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
  SideSlabNormals          (1:3,1:3,          newBCSideID) =DummySideSlabNormals         (1:3,1:3,           iSide)
  SideSlabIntervals       (1:6,              newBCSideID) =DummySideSlabIntervals      (1:6,               iSide)
  BoundingBoxIsEmpty   (                  newBCSideID) =DummyBoundingBoxIsEmpty  (                   iSide)
  SideType(newBCSideID) = DummySideType(newBCSideID)
  SideDistance(newBCSideID) = DummySideDistance(newBCSideID)
  SideNormVec(1:3,newBCSideID) = DummySideNormVec(1:3,newBCSideID)
END DO ! iSide

! create new mapping
SDEALLOCATE(PartBCSideList)
ALLOCATE(PartBCSideList(1:nTotalSides))
PartBCSideList=-1

newBCSideID=0
DO iSide=1,nBCSides
  newBCSideID=newBCSideID+1
  PartBCSideList(iSide)=newBCSideID
END DO

DO iSide=nBCSides+1,nSides+nPartPeriodicSides
  IF(BC(iSide).EQ.0) CYCLE
  newBCSideID=newBCSideID+1
  PartBCSideList(iSide)=newBCSideID
END DO ! iSide

DO iSide=nSides+nPartPeriodicSides+1,nTotalSides
  newBCSideID=newBCSideID+1
  PartBCSideList(iSide)=newBCSideID
END DO

! deallocate dummy buffer
DEALLOCATE(DummyBezierControlPoints3D)
DEALLOCATE(DummySideSlabNormals)
DEALLOCATE(DummySideSlabIntervals)
DEALLOCATE(DummyBoundingBoxIsEmpty)
DEALLOCATE(DummySideType)
DEALLOCATE(DummySideDistance)
DEALLOCATE(DummySideNormVec)


END SUBROUTINE ReShapeBezierSides


SUBROUTINE MapRegionToElem() 
!----------------------------------------------------------------------------------------------------------------------------------!
! map a particle region to element
! check only element barycenter, nothing else
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,     ONLY:NbrOfRegions, RegionBounds,GEO
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
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
      CALL abort(&
__STAMP__&
,'Defined regions are overlapping')
    END IF
  END DO ! iRegions=1,NbrOfRegions
END DO ! iElem=1,PP_nElems


END SUBROUTINE MapRegionToElem


SUBROUTINE PointToExactElement(X_In,Element,isInSide,doHalo)                                                         
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,     ONLY:Geo
USE MOD_Particle_Mesh_Vars,     ONLY:epsOneCell
USE MOD_Particle_Tracking_Vars, ONLY:ListDistance,Distance
USE MOD_Eval_xyz,               ONLY:GetPositionInRefElem
USE MOD_Utils,                  ONLY:InsertionSort !BubbleSortID
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)                   :: X_in(3)
LOGICAL,INTENT(IN)                :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL,INTENT(OUT)                :: isInside
INTEGER,INTENT(OUT)                :: Element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                              :: xi(1:3)
!REAL,PARAMETER                    :: eps=1e-8 ! same value as in eval_xyz_elem
!REAL,PARAMETER                    :: eps2=1e-3
!REAL                              :: epsOne,OneMeps
!===================================================================================================================================

!epsOne=1.0+epsInCell
!OneMeps=1.0-eps
isInside = .FALSE.
IF ( (X_in(1).LT.GEO%xmin).OR.(X_in(1).GT.GEO%xmax).OR. &
     (X_in(2).LT.GEO%ymin).OR.(X_in(2).GT.GEO%ymax).OR. &
     (X_in(3).LT.GEO%zmin).OR.(X_in(3).GT.GEO%zmax)) THEN
   RETURN
END IF

! --- get background mesh cell of particle
CellX = CEILING((X_in(1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
CellX = MIN(GEO%FIBGMimax,CellX)                             
CellY = CEILING((X_in(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MIN(GEO%FIBGMjmax,CellY) 
CellZ = CEILING((X_in(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MIN(GEO%FIBGMkmax,CellZ)

!--- check all cells associated with this beckground mesh cell
nBGMElems=GEO%FIBGM(CellX,CellY,CellZ)%nElem

! get closest element barycenter
Distance=-1.
ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  Distance(iBGMElem)=(X_in(1)-ElemBaryNGeo(1,ElemID))*(X_in(1)-ElemBaryNGeo(1,ElemID)) &
                    +(X_in(2)-ElemBaryNGeo(2,ElemID))*(X_in(2)-ElemBaryNGeo(2,ElemID)) &
                    +(X_in(3)-ElemBaryNGeo(3,ElemID))*(X_in(3)-ElemBaryNGeo(3,ElemID)) 
  Distance(iBGMElem)=SQRT(Distance(iBGMElem))
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element  
Element=-1
DO iBGMElem=1,nBGMElems
  ElemID=ListDistance(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  CALL GetPositionInRefElem(X_in(1:3),xi,ElemID)
  IF(ALL(ABS(Xi).LE.epsOneCell(ElemID))) THEN ! particle inside
    isInSide=.TRUE.
    Element=ElemID
    EXIT
  END IF
END DO ! iBGMElem

END SUBROUTINE PointToExactElement


SUBROUTINE BuildElementBasis()
!================================================================================================================================
! build the element local basis system 
! origin is located at xi=(0,0,0)^T
! each local coord system is pointing to an element side
!================================================================================================================================
USE MOD_Globals!,                  ONLY:CROSS
USE MOD_Preproc
USE MOD_Mesh_Vars,                ONLY:NGeo,XCL_NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D
USE MOD_Basis,                    ONLY:DeCasteljauInterpolation
USE MOD_Particle_Mesh_Vars,       ONLY:XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
USE MOD_Mesh_Vars,                ONLY:ElemBaryNgeo
USE MOD_Particle_Tracking_Vars,   ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,       ONLY:nTotalElems,PartElemToSide
USE MOD_Basis,                    ONLY:LagrangeInterpolationPolys
USE MOD_PICDepo_Vars,             ONLY:DepositionType,r_sf,ElemRadius2_sf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem,SideID,i,j,k,ilocSide, ALLOCSTAT
REAL                    :: Xi(3),XPos(3),Radius
REAL                    :: Lag(1:3,0:NGeo)
!================================================================================================================================

ElemRadiusNGeo=0.
ElemRadius2NGeo=0.
DO iElem=1,nTotalElems
  ! get point on each side 
  IF((iElem.LE.PP_nElems).OR.(DoRefMapping))THEN
    ! xi plus
    Xi=(/1.0,0.0,0.0/)
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,1,iElem)=xPos
    ! eta plus
    Xi=(/0.0,1.0,0.0/)
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,2,iElem)=xPos
    ! zeta plus
    Xi=(/0.0,0.0,1.0/)
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,3,iElem)=xPos
    ! xi minus
    Xi=(/-1.0,0.0,0.0/)
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,4,iElem)=xPos
    ! eta minus
    Xi=(/0.0,-1.0,0.0/)
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,5,iElem)=xPos
    ! zeta minus
    Xi=(/0.0,0.0,-1.0/)
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,6,iElem)=xPos
  ELSE ! compute particle position in physical space
    Xi=(/0.0,0.0,0.0/)
    SideID = PartElemToSide(1,XI_PLUS,iElem)
    CALL DeCasteljauInterpolation(NGeo,Xi(1:2),SideID,XiEtaZetaBasis(1:3,1,iElem))
    SideID = PartElemToSide(1,ETA_PLUS,iElem)
    CALL DeCasteljauInterpolation(NGeo,Xi(1:2),SideID,XiEtaZetaBasis(1:3,2,iElem))
    SideID = PartElemToSide(1,ZETA_PLUS,iElem)
    CALL DeCasteljauInterpolation(NGeo,Xi(1:2),SideID,XiEtaZetaBasis(1:3,3,iElem))
    SideID = PartElemToSide(1,XI_MINUS,iElem)
    CALL DeCasteljauInterpolation(NGeo,Xi(1:2),SideID,XiEtaZetaBasis(1:3,4,iElem))
    SideID = PartElemToSide(1,ETA_MINUS,iElem)
    CALL DeCasteljauInterpolation(NGeo,Xi(1:2),SideID,XiEtaZetaBasis(1:3,5,iElem))
    SideID = PartElemToSide(1,ZETA_MINUS,iElem)
    CALL DeCasteljauInterpolation(NGeo,Xi(1:2),SideID,XiEtaZetaBasis(1:3,6,iElem))
  END IF ! no ref mapping
  ! compute vector from each barycenter to sidecenter
  XiEtaZetaBasis(:,1,iElem)=XiEtaZetaBasis(:,1,iElem)-ElemBaryNGeo(:,iElem)
  XiEtaZetaBasis(:,2,iElem)=XiEtaZetaBasis(:,2,iElem)-ElemBaryNGeo(:,iElem)
  XiEtaZetaBasis(:,3,iElem)=XiEtaZetaBasis(:,3,iElem)-ElemBaryNGeo(:,iElem)
  XiEtaZetaBasis(:,4,iElem)=XiEtaZetaBasis(:,4,iElem)-ElemBaryNGeo(:,iElem)
  XiEtaZetaBasis(:,5,iElem)=XiEtaZetaBasis(:,5,iElem)-ElemBaryNGeo(:,iElem)
  XiEtaZetaBasis(:,6,iElem)=XiEtaZetaBasis(:,6,iElem)-ElemBaryNGeo(:,iElem)
  ! compute length
  slenXiEtaZetaBasis(1,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,1,iElem),XiEtaZetaBasis(:,1,iElem))
  slenXiEtaZetaBasis(2,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,2,iElem),XiEtaZetaBasis(:,2,iElem))
  slenXiEtaZetaBasis(3,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,3,iElem),XiEtaZetaBasis(:,3,iElem))
  slenXiEtaZetaBasis(4,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,4,iElem),XiEtaZetaBasis(:,4,iElem))
  slenXiEtaZetaBasis(5,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,5,iElem),XiEtaZetaBasis(:,5,iElem))
  slenXiEtaZetaBasis(6,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,6,iElem),XiEtaZetaBasis(:,6,iElem))

  Radius=0.
  IF(DoRefMapping)THEN ! thats not the bounding box, caution, this box is to small!
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=XCL_NGeo(:,i,j,k,iElem)-ElemBaryNGeo(:,iElem)
          Radius=MAX(Radius,SQRT(DOT_PRODUCT(xPos,xPos)))      
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
  ELSE
    ! required :(
   ! IF(iElem.GT.PP_nElems) CYCLE
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(SideID.EQ.-1) CYCLE
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=BezierControlPoints3D(:,i,j,SideID)-ElemBaryNGeo(:,iElem)
          Radius=MAX(Radius,SQRT(DOT_PRODUCT(xPos,xPos)))      
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO ! ilocSide
  END IF
  !ElemRadiusNGeo(iElem)=Radius
  ! elem radius containts 10% tolerance because we are not using the beziercontrolpoints
  ElemRadiusNGeo(iElem)=Radius
  IF(DoRefMapping)THEN
    !ElemRadius2NGeo(iElem)=(Radius*1.10)*(Radius*1.10)
    ElemRadius2NGeo(iElem)=(Radius*1.02)*(Radius*1.02)
  ELSE
    ElemRadius2NGeo(iElem)=Radius*Radius
  END IF
END DO ! iElem

IF (TRIM(DepositionType).EQ.'shape_function_simple')THEN
  ALLOCATE(ElemRadius2_sf(1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__ &
,' Cannot allocate ElemRadius2_sf!')
  DO iElem=1,PP_nElems
    ElemRadius2_sf(iElem)=(ElemRadiusNGeo(iElem)+r_sf)*(ElemRadiusNGeo(iElem)+r_sf)
  END DO ! iElem=1,PP_nElems
END IF


END SUBROUTINE BuildElementBasis


SUBROUTINE MapElemToFIBGM() 
!----------------------------------------------------------------------------------------------------------------------------------!
! here, the FIBGM range for each element is stored
! short list for intersection tracking, longer list for ref mapping tracking
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,nTotalElems
USE MOD_Mesh_Vars,              ONLY:XCL_NGeo
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ALLOCSTAT,iElem,lastElem
REAL              :: xmin,ymin,zmin,xmax,ymax,zmax
INTEGER           :: BGMimax,BGMimin,BGMjmax,BGMjmin,BGMkmax,BGMkmin
INTEGER           :: BGMCellXmax,BGMCellXmin,BGMCellYmax,BGMCellYmin,BGMCellZmax,BGMCellZmin
!===================================================================================================================================

!IF(.NOT.DoRefMapping) RETURN

IF(DoRefMapping) THEN
  LastElem=nTotalElems
ELSE
  LastElem=PP_nElems
END IF

ALLOCATE(GEO%ElemToFIBGM(1:6,1:LastElem),STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate GEO%ElemToFIBGM!')

! because I copy and past
BGMimax=GEO%FIBGMimax
BGMimin=GEO%FIBGMimin
BGMjmax=GEO%FIBGMjmax
BGMjmin=GEO%FIBGMjmin
BGMkmax=GEO%FIBGMkmax
BGMkmin=GEO%FIBGMkmin

DO iElem=1,LastElem
  xmin=HUGE(1.)
  ymin=HUGE(1.)
  zmin=HUGE(1.)
  xmax=-HUGE(1.)
  ymax=-HUGE(1.)
  zmax=-HUGE(1.)
  xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
  xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
  ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
  ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
  zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
  zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  !--- find minimum and maximum BGM cell for current element
  IF(GEO%nPeriodicVectors.EQ.0)THEN
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
  ELSE
    ! here fancy stuff, because element could be wide out of element range
    BGMCellXmax = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
    BGMCellXmax = MAX(MIN(BGMCellXmax,BGMimax),BGMimin)
    BGMCellXmin = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
    BGMCellXmin = MIN(MAX(BGMCellXmin,BGMimin),BGMimax)
    BGMCellYmax = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
    BGMCellYmax = MAX(MIN(BGMCellYmax,BGMjmax),BGMjmin)
    BGMCellYmin = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
    BGMCellYmin = MIN(MAX(BGMCellYmin,BGMjmin),BGMjmax)
    BGMCellZmax = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
    BGMCellZmax = MAX(MIN(BGMCellZmax,BGMkmax),BGMkmin)
    BGMCellZmin = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
    BGMCellZmin = MIN(MAX(BGMCellZmin,BGMkmin),BGMkmax)
  END IF
  GEO%ElemToFIBGM(1,iElem)=BGMCellXmin  
  GEO%ElemToFIBGM(3,iElem)=BGMCellYmin  
  GEO%ElemToFIBGM(5,iElem)=BGMCellZmin  

  GEO%ElemToFIBGM(2,iElem)=BGMCellXmax  
  GEO%ElemToFIBGM(4,iElem)=BGMCellYmax  
  GEO%ElemToFIBGM(6,iElem)=BGMCellZmax  
END DO ! iElem=1,nTotalElems

END SUBROUTINE MapElemToFIBGM


SUBROUTINE CountPartsPerElem(ResetNumberOfParticles)
!===================================================================================================================================
! count number of particles in element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_LoadBalance_Vars,        ONLY: nPartsPerElem
USE MOD_Particle_Vars,           ONLY: PDM,PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: ResetNumberOfParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPart, ElemID
!===================================================================================================================================
! DO NOT NULL this here, if e.g. this routine is called in between RK-stages in which particles are created
IF(ResetNumberOfParticles)THEN
  nPartsPerElem=0
END IF
! loop over all particles and add them up
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    ElemID = PEM%Element(iPart)
    IF(ElemID.LE.PP_nElems)THEN
      nPartsPerElem(ElemID)=nPartsPerElem(ElemID)+1
    END IF
  END IF
END DO ! iPart=1,PDM%ParticleVecLength

END SUBROUTINE CountPartsPerElem


SUBROUTINE CheckIfCurvedElem(IsCurved,XCL_NGeo)
!===================================================================================================================================
! check if element is curved
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Mesh_Vars,             ONLY:NGeo,Vdm_CLNGeo1_CLNGeo
USE MOD_ChangeBasis,           ONLY:changeBasis3D
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)      :: XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)  :: IsCurved
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                 :: XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0:NGeo)
INTEGER              :: NGeo3
!===================================================================================================================================

IsCurved=.FALSE.

! fill dummy
XCL_NGeo1(1:3,0,0,0) = XCL_NGeo(1:3, 0  , 0  , 0  )
XCL_NGeo1(1:3,1,0,0) = XCL_NGeo(1:3,NGeo, 0  , 0  )
XCL_NGeo1(1:3,0,1,0) = XCL_NGeo(1:3, 0  ,NGeo, 0  )
XCL_NGeo1(1:3,1,1,0) = XCL_NGeo(1:3,NGeo,NGeo, 0  )
XCL_NGeo1(1:3,0,0,1) = XCL_NGeo(1:3, 0  , 0  ,NGeo)
XCL_NGeo1(1:3,1,0,1) = XCL_NGeo(1:3,NGeo, 0  ,NGeo)
XCL_NGeo1(1:3,0,1,1) = XCL_NGeo(1:3, 0  ,NGeo,NGeo)
XCL_NGeo1(1:3,1,1,1) = XCL_NGeo(1:3,NGeo,NGeo,NGeo)

CALL ChangeBasis3D(3,1,NGeo,Vdm_CLNGeo1_CLNGeo,XCL_NGeo1,XCL_NGeoNew)
NGeo3=(NGeo+1)*(NGeo+1)*(NGeo+1)

! check 3D points
CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeo,IsCurved)

IF(.NOT.IsCurved)THEN
  ! set all elem sides to blabla
END IF

END SUBROUTINE CheckIfCurvedElem


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

SUBROUTINE InitElemBoundingBox() 
!===================================================================================================================================
! init of tight elem bounding box, constructed via beziercontrolpoints
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
!USE MOD_Particle_Surfaces,       ONLY:GetElemSlabNormalsAndIntervals
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:ExchangeBezierControlPoints3d
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#ifdef PARTICLES
#ifdef MPI
! first communicate the bezierControlPoints (slave information is missing)
CALL ExchangeBezierControlPoints3D()
#endif /*MPI*/
#endif /*PARTICLES*/

END SUBROUTINE InitElemBoundingBox


SUBROUTINE InsideElemBoundingBox(ParticlePosition,ElemID,InSide)
!================================================================================================================================
! check if the particles is inside the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Particle_Surfaces_Vars,  ONLY:ElemSlabNormals,ElemSlabIntervals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: ParticlePosition
INTEGER,INTENT(IN)                   :: ElemID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                  :: Inside
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: x,y,z,P(3)
!================================================================================================================================
P=ParticlePosition-ElemSlabNormals(1:3,0,ElemID)
! y is perpendicular to xi & eta directions --> check first, smallest intervall
y=DOT_PRODUCT(P,ElemSlabNormals(:,2,ElemID))
!IF((y.LT.ElemSlabIntervals(3,ElemID)-epsilontol).OR.(y.GT.ElemSlabIntervals(4,ElemID)+epsilontol))THEN
IF((y.LT.ElemSlabIntervals(3,ElemID)).OR.(y.GT.ElemSlabIntervals(4,ElemID)))THEN
  Inside=.FALSE.
  RETURN
END IF
! than xi
x=DOT_PRODUCT(P,ElemSlabNormals(:,1,ElemID))
!IF((x.LT.ElemSlabIntervals(1,ElemID)-epsilontol).OR.(x.GT.ElemSlabIntervals(2,ElemID)+epsilontol))THEN
IF((x.LT.ElemSlabIntervals(1,ElemID)).OR.(x.GT.ElemSlabIntervals(2,ElemID)))THEN
  Inside=.FALSE.
  RETURN
END IF
! than eta
z=DOT_PRODUCT(P,ElemSlabNormals(:,3,ElemID))
!IF((z.LT.ElemSlabIntervals(5,ElemID)-epsilontol).OR.(z.GT.ElemSlabIntervals(6,ElemID)+epsilontol))THEN
IF((z.LT.ElemSlabIntervals(5,ElemID)).OR.(z.GT.ElemSlabIntervals(6,ElemID)))THEN
  Inside=.FALSE.
  RETURN
END IF
Inside=.TRUE.
END SUBROUTINE InsideElemBoundingBox


SUBROUTINE GetElemAndSideType()
!===================================================================================================================================
! get the element and side type of each element (no mpi communication, this is done in exchange_halo_geometry)
! 1) Get Elem Type
! 2) Get Side Type
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_Vars,             ONLY:DoRefMapping
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D,BoundingBoxIsEmpty,SideType,SideNormVec,SideDistance
USE MOD_Particle_Mesh_Vars,                 ONLY:nTotalSides,nTotalElems,SidePeriodicType
USE MOD_Particle_Mesh_Vars,                 ONLY:ElemType,nPartSides
USE MOD_Mesh_Vars,                          ONLY:CurvedElem,XCL_NGeo,Vdm_CLNGeo1_CLNGeo,NGeo,Vdm_CLNGeo1_CLNGeo,ElemBaryNGeo
USE MOD_Particle_Mesh_Vars,                 ONLY:PartElemToSide,PartBCSideList,nTotalBCSides,GEO
USE MOD_ChangeBasis,                        ONLY:changeBasis3D
USE MOD_Particle_Vars,                      ONLY:PartMeshHasPeriodicBCs
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: iElem
INTEGER                                  :: iSide,SideID,TrueSideID,ilocSide,BCSideID
INTEGER                                  :: flip
REAL,DIMENSION(1:3)                      :: v1,v2,v3
LOGICAL,ALLOCATABLE                      :: SideIsDone(:)
REAL                                     :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                                     :: XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0:NGeo),Vec1(1:3)
INTEGER                                  :: NGeo3,NGeo2,PVID
REAL                                     :: XCL_NGeoSideNew(1:3,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoSideOld(1:3,0:NGeo,0:NGeo)
LOGICAL                                  :: isCurvedSide,isRectangular
REAL                                     :: ScalarProduct
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Get Element and Side Type ...'

! elements
ALLOCATE(CurvedElem(1:nTotalElems))
CurvedElem=.FALSE.
IF (.NOT.DoRefMapping) THEN
  ALLOCATE(ElemType(1:nTotalElems))
  ElemType=-1
END IF

! sides
IF(DoRefMapping)THEN
  ALLOCATE( SideType(nTotalBCSides)        &
          , SideDistance(nTotalBCSides)    &
          , SideIsDone(nTotalSides)        &
          , SideNormVec(1:3,nTotalBCSides) )
ELSE
  ALLOCATE( SideType(nTotalSides)        &
          , SideDistance(nTotalSides)    &
          , SideIsDone(nTotalSides)      &
          , SideNormVec(1:3,nTotalSides) )
END IF
SideIsDone=.FALSE.
SideType=-1

SideDistance=-0.
SideNormVec=0.

NGeo2=(NGeo+1)*(NGeo+1)
NGeo3=NGeo2*(NGeo+1)

! decide if element is (bi-)linear or curved
! decide if sides are planar-rect, planar-nonrect, planar-curved, bilinear or curved 
DO iElem=1,nTotalElems
  ! 1) check if elem is curved
  !   a) get the coordinates of the eight nodes of the hexahedral
  XCL_NGeo1(1:3,0,0,0) = XCL_NGeo(1:3, 0  , 0  , 0  ,iElem)
  XCL_NGeo1(1:3,1,0,0) = XCL_NGeo(1:3,NGeo, 0  , 0  ,iElem)
  XCL_NGeo1(1:3,0,1,0) = XCL_NGeo(1:3, 0  ,NGeo, 0  ,iElem)
  XCL_NGeo1(1:3,1,1,0) = XCL_NGeo(1:3,NGeo,NGeo, 0  ,iElem)
  XCL_NGeo1(1:3,0,0,1) = XCL_NGeo(1:3, 0  , 0  ,NGeo,iElem)
  XCL_NGeo1(1:3,1,0,1) = XCL_NGeo(1:3,NGeo, 0  ,NGeo,iElem)
  XCL_NGeo1(1:3,0,1,1) = XCL_NGeo(1:3, 0  ,NGeo,NGeo,iElem)
  XCL_NGeo1(1:3,1,1,1) = XCL_NGeo(1:3,NGeo,NGeo,NGeo,iElem)

  !  b) interpolate from the nodes to NGeo
  !     Compare the bi-linear mapping with the used mapping
  !     For NGeo=1, this should always be true, because the mappings are identical
  CALL ChangeBasis3D(3,1,NGeo,Vdm_CLNGeo1_CLNGeo,XCL_NGeo1,XCL_NGeoNew)
  ! check the coordinates of all Chebychev-Lobatto geometry points between the bi-linear and used
  ! mapping
  CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),CurvedElem(iElem))

  ! 2) check sides
  ! loop over all 6 sides of element
  ! a) check if the sides are straight
  ! b) use curved information to decide side type
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    flip  =PartElemToSide(E2S_FLIP,ilocSide,iElem)
    IF (SideID.LE.0) CYCLE
    IF (SideIsDone(SideID)) CYCLE
    IF(DoRefMapping)THEN
      TrueSideID=PartBCSideList(SideID)
      IF(TrueSideID.EQ.-1)CYCLE
    ELSE
      TrueSideID=SideID
    END IF
    IF(.NOT.CurvedElem(iElem))THEN
      ! linear element
      IF(BoundingBoxIsEmpty(SideID))THEN
        v1=(-BezierControlPoints3D(:,0,0   ,    SideID)+BezierControlPoints3D(:,NGeo,0   ,    SideID)   &
            -BezierControlPoints3D(:,0,NGeo,    SideID)+BezierControlPoints3D(:,NGeo,NGeo,    SideID) )
        
        v2=(-BezierControlPoints3D(:,0,0   ,    SideID)-BezierControlPoints3D(:,NGeo,0   ,    SideID)   &
            +BezierControlPoints3D(:,0,NGeo,    SideID)+BezierControlPoints3D(:,NGeo,NGeo,    SideID) )
        SideNormVec(:,TrueSideID) = CROSSNORM(v1,v2)
        v1=0.25*(BezierControlPoints3D(:,0,0      ,SideID)     &
                +BezierControlPoints3D(:,NGeo,0   ,SideID)  &
                +BezierControlPoints3D(:,0,NGeo   ,SideID)  &
                +BezierControlPoints3D(:,NGeo,NGeo,SideID))
        ! check if normal vector points outwards
        v2=v1-ElemBaryNGeo(:,iElem)
        IF(flip.EQ.0)THEN
          IF(DOT_PRODUCT(v2,SideNormVec(:,TrueSideID)).LT.0) SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID) 
        ELSE IF(flip.EQ.-1)THEN
          SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID) 
          PartElemToSide(E2S_FLIP,ilocSide,iElem) = 0
        ELSE
          IF(DOT_PRODUCT(v2,SideNormVec(:,TrueSideID)).GT.0) SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID)
        END IF
        SideDistance(TrueSideID)=DOT_PRODUCT(v1,SideNormVec(:,TrueSideID))
        ! check if it is rectangular
        isRectangular=.TRUE.
        v1=UNITVECTOR(BezierControlPoints3D(:,0   ,NGeo,SideID)-BezierControlPoints3D(:,0   ,0   ,SideID))
        v2=UNITVECTOR(BezierControlPoints3D(:,NGeo,0   ,SideID)-BezierControlPoints3D(:,0   ,0   ,SideID))
        v3=UNITVECTOR(BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID))
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        IF(isRectangular)THEN
          v1=UNITVECTOR(BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID))
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        END IF
        IF(isRectangular)THEN
          SideType(TrueSideID)=PLANAR_RECT
        ELSE
          SideType(TrueSideID)=PLANAR_NONRECT
        END IF
      ELSE
        v1=(-BezierControlPoints3D(:,0,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
            -BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
        v2=(-BezierControlPoints3D(:,0,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
            +BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
        SideNormVec(:,TrueSideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
        SideType(TrueSideID)=BILINEAR
      END IF
    ELSE
      ! possible curved face
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        XCL_NGeoSideOld=XCL_NGeo   (1:3,0,0:NGeo,0:NGeo,iElem)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0,0:NGeo,0:NGeo)
      CASE(XI_PLUS)
        XCL_NGeoSideOld=XCL_NGeo   (1:3,NGeo,0:NGeo,0:NGeo,iElem)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,NGeo,0:NGeo,0:NGeo)
      CASE(ETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeo   (1:3,0:NGeo,0,0:NGeo,iElem)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0,0:NGeo)
      CASE(ETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeo   (1:3,0:NGeo,NGeo,0:NGeo,iElem)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,NGeo,0:NGeo)
      CASE(ZETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeo   (1:3,0:NGeo,0:NGeo,0,iElem)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0)
      CASE(ZETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeo   (1:3,0:NGeo,0:NGeo,NGeo,iElem)
        XCL_NGeoSideNew=XCL_NGeoNEw(1:3,0:NGeo,0:NGeo,NGeo)
      END SELECT
      CALL PointsEqual(NGeo2,XCL_NGeoSideNew,XCL_NGeoSideOld,isCurvedSide)
      IF(isCurvedSide)THEn
        IF(BoundingBoxIsEmpty(SideID))THEN
          SideType(TrueSideID)=PLANAR_CURVED
          v1=(-BezierControlPoints3D(:,0,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
              -BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
          
          v2=(-BezierControlPoints3D(:,0,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
              +BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
          SideNormVec(:,TrueSideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints3D(:,0,0,SideID)     &
                  +BezierControlPoints3D(:,NGeo,0,SideID)  &
                  +BezierControlPoints3D(:,0,NGeo,SideID)  &
                  +BezierControlPoints3D(:,NGeo,NGeo,SideID))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,TrueSideID)).LT.0) SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID) 
          ELSE IF(flip.EQ.-1)THEN
            SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID) 
            PartElemToSide(E2S_FLIP,ilocSide,iElem) = 0
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,TrueSideID)).GT.0) SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID)
          END IF
          SideDistance(TrueSideID)=DOT_PRODUCT(v1,SideNormVec(:,TrueSideID))
        ELSE
          SideType(TrueSideID)=CURVED
        END IF
      ELSE
        IF(BoundingBoxIsEmpty(SideID))THEN
          v1=(-BezierControlPoints3D(:,0,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
              -BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
          
          v2=(-BezierControlPoints3D(:,0,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
              +BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
          SideNormVec(:,TrueSideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints3D(:,0,0,SideID)     &
                  +BezierControlPoints3D(:,NGeo,0,SideID)  &
                  +BezierControlPoints3D(:,0,NGeo,SideID)  &
                  +BezierControlPoints3D(:,NGeo,NGeo,SideID))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,TrueSideID)).LT.0) SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID) 
          ELSE IF(flip.EQ.-1)THEN
            SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID) 
            PartElemToSide(E2S_FLIP,ilocSide,iElem) = 0
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,TrueSideID)).GT.0) SideNormVec(:,TrueSideID)=-SideNormVec(:,TrueSideID)
          END IF
          SideDistance(TrueSideID)=DOT_PRODUCT(v1,SideNormVec(:,TrueSideID))
          ! check if it is rectangular
          isRectangular=.TRUE.
          v1=UNITVECTOR(BezierControlPoints3D(:,0   ,NGeo,SideID)-BezierControlPoints3D(:,0   ,0   ,SideID))
          v2=UNITVECTOR(BezierControlPoints3D(:,NGeo,0   ,SideID)-BezierControlPoints3D(:,0   ,0   ,SideID))
          v3=UNITVECTOR(BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID))
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
          IF(isRectangular)THEN
            v1=UNITVECTOR(BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID))
            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
          END IF
          IF(isRectangular)THEN
            SideType(TrueSideID)=PLANAR_RECT
          ELSE
            SideType(TrueSideID)=PLANAR_NONRECT
          END IF
        ELSE
          v1=(-BezierControlPoints3D(:,0,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
              -BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
          v2=(-BezierControlPoints3D(:,0,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
              +BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )
          SideNormVec(:,TrueSideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
          SideType(TrueSideID)=BILINEAR
        END IF
      END IF
    END IF
    SideIsDone(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iElem=1,nTotalElems

! sanity check for side periodic type
DO iSide=1,nPartSides
  IF(DoRefmapping)THEN
    BCSideID  =PartBCSideList(iSide)
    IF(BCSideID.LE.0) CYCLE
  ELSE
    BCSideID  =iSide
  END IF
  PVID=SidePeriodicType(iSide)
  IF(PVID.EQ.0) CYCLE
  IF(.NOT.PartMeshHasPeriodicBCs) CYCLE
  Vec1=SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
  ScalarProduct=DOT_PRODUCT(SideNormVec(1:3,BCSideID),Vec1)
  IF(ALMOSTEQUAL(ScalarProduct,GEO%PeriodicVectorsLength(ABS(PVID))))THEN
    SidePeriodicType(iSide)=-SidePeriodicType(iSide)
  ELSEIF(.NOT.ALMOSTEQUAL(ScalarProduct,-GEO%PeriodicVectorsLength(ABS(PVID))))THEN
    CALL abort(&
__STAMP__&
        , ' Missalignment between SideNormVec and PeriodicVector!',ABS(PVID),ScalarProduct)
  END IF
END DO ! iSide=1,nPartSides

! fill Element type checking sides
IF (.NOT.DoRefMapping) THEN
  DO iElem=1,nTotalElems
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      SELECT CASE(SideType(SideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT)
        IF (ElemType(iElem).GE.1) THEN
          CYCLE
        ELSE
          ElemType(iElem) = 1
        END IF
      CASE(BILINEAR)
        IF (ElemType(iElem).GE.2) THEN
          CYCLE
        ELSE
          ElemType(iElem) = 2
        END IF
      CASE(PLANAR_CURVED,CURVED)
        ElemType(iElem) = 3
        EXIT
      END SELECT
    END DO ! ilocSide=1,6
  END DO ! iElem=1,nTotalElems
END IF

END SUBROUTINE GetElemAndSideType


SUBROUTINE GetBCElemMap()
!===================================================================================================================================
! 1) Add BC and Halo-BC sides in halo_eps distance to a certain element
! 2) build epsOneCell for each element
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_IO_HDF5,                            ONLY:AddToElemData,ElementOut
USE MOD_Particle_Tracking_Vars,             ONLY:DoRefMapping
USE MOD_Mesh_Vars,                          ONLY:XCL_NGeo,nSides,NGeo,nBCSides,sJ,BC,nElems
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D
USE MOD_Particle_Mesh_Vars,                 ONLY:nTotalSides,IsTracingBCElem,nTotalElems,nTotalBCElems
USE MOD_Particle_Mesh_Vars,                 ONLY:TracingBCInnerSides,TracingBCTotalSides
USE MOD_Particle_Mesh_Vars,                 ONLY:PartElemToSide,BCElem,PartSideToElem,PartBCSideList,GEO
USE MOD_Particle_Mesh_Vars,                 ONLY:RefMappingEps,epsOneCell
USE MOD_Particle_Surfaces_Vars,             ONLY:sVdm_Bezier
USE MOD_Particle_MPI_Vars,                  ONLY:halo_eps,halo_eps2
USE MOD_ChangeBasis,                        ONLY:ChangeBasis2D
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: iElem,firstBezierPoint,lastBezierPoint
INTEGER                                  :: iSide,p,q,SideID,ilocSide,BCSideID2,BCSideID
INTEGER                                  :: nSideCount, s,r
INTEGER,ALLOCATABLE                      :: SideIndex(:)
REAL,DIMENSION(1:3)                      :: v1,NodeX,Vec1
REAL,DIMENSION(1:3,0:NGeo,0:NGeo)        :: xNodes
INTEGER                                  :: nLoop,iTest,nTest
REAL                                     :: scaleJ, Distance ,maxScaleJ,dx,dy,dz
LOGICAL                                  :: fullMesh, leave
!===================================================================================================================================
ALLOCATE(IsTracingBCElem(nTotalElems))
IsTracingBCElem=.FALSE.
IF(DoRefMapping)THEN
  ALLOCATE(TracingBCInnerSides(nTotalElems))
  TracingBCInnerSides=0
  ALLOCATE(TracingBCTotalSides(nTotalElems))
  TracingBCTotalSides=0
END IF

! decide if element:  
! DoRefMapping=T
! a) HAS own bc faces
! b) HAS bc-face in halo_eps distance
! DoRefMapping=F
! a) HAS own bc faces
IF(DoRefMapping)THEN
  ! mark elements as bc element if they have a local-BC side
  nTotalBCElems=0
  DO iElem=1,nTotalElems
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF (SideID.LE.0) CYCLE
      IF((SideID.LE.nBCSides).OR.(SideID.GT.nSides))THEN
        IF(.NOT.IsTracingBCElem(iElem))THEN
          IsTracingBCElem(iElem)=.TRUE.
          nTotalBCElems=nTotalBCElems+1
        END IF ! count only single
      END IF
    END DO ! ilocSide
  END DO ! iElem

  ! for simplifications
  ! get distance of diagonal of mesh
  V1(1) = GEO%xmaxglob-GEO%xminglob
  V1(2) = GEO%ymaxglob-GEO%yminglob
  V1(3) = GEO%zmaxglob-GEO%zminglob
  Distance=DOT_PRODUCT(V1,V1)
  fullMesh=.FALSE.
  ! build list with elements in halo-eps vicinity around bc-elements
  IF(Distance.LE.halo_eps2) fullMesh=.TRUE.  
  ! allocate the types for the element to bc-side mapping
  ALLOCATE( BCElem(1:nTotalElems) )
  ALLOCATE( SideIndex(1:nTotalSides) )
  ! for fullMesh, each element requires ALL BC faces
  IF(fullMesh)THEN
    DO iElem=1,nTotalElems
      ! mark my sides
      BCElem(iElem)%nInnerSides=0
      DO ilocSide=1,6
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        IF(SideID.LE.0) CYCLE
        IF(PartBCSideList(SideID).EQ.-1) CYCLE
        BCElem(iElem)%nInnerSides = BCElem(iElem)%nInnerSides+1
      END DO ! ilocSide=1,6
      BCElem(iElem)%lastSide=BCElem(iElem)%nInnerSides
      ! loop over all sides, exclusive of own sides
      SideIndex=0
      DO iSide=1,nTotalSides
        ! only bc sides
        BCSideID  =PartBCSideList(iSide)
        IF(BCSideID.EQ.-1) CYCLE
        ! ignore sides of the same element
        IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE
        IF(SideIndex(iSide).EQ.0)THEN
          BCElem(iElem)%lastSide=BCElem(iElem)%lastSide+1
          SideIndex(iSide)=BCElem(iElem)%lastSide
        END IF
      END DO ! iSide=1,nTotalSides
      IF(BCElem(iElem)%lastSide.EQ.0) CYCLE
      ! set true, only required for elements without an own bc side
      IsTracingBCElem(iElem)=.TRUE.
      ! allocate complete side list
      ALLOCATE( BCElem(iElem)%BCSideID(BCElem(iElem)%lastSide) )
      ! 1) inner sides
      nSideCount=0
      IF(BCElem(iElem)%nInnerSides.GT.0)THEN
        DO ilocSide=1,6
          SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
          IF(SideID.LE.0) CYCLE
          BCSideID=PartBCSideList(SideID)
          IF(BCSideID.LE.0) CYCLE
          nSideCount=nSideCount+1
          BCElem(iElem)%BCSideID(nSideCount)=SideID
        END DO ! ilocSide
      END IF ! nInnerSides.GT.0
      ! 2) outer sides
      DO iSide=1,nTotalSides
        IF(SideIndex(iSide).GT.0)THEN
          nSideCount=nSideCount+1
          BCElem(iElem)%BCSideID(nSideCount)=iSide !iSide
        END IF
      END DO  ! iSide
    END DO ! iElem=1,nTotalElems
  ELSE ! .NOT. fullMesh
    ! each element requires only the sides in its halo region
    DO iElem=1,nTotalElems
      ! mark my sides
      BCElem(iElem)%nInnerSides=0
      DO ilocSide=1,6
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        IF(SideID.LE.0) CYCLE
        IF(PartBCSideList(SideID).EQ.-1) CYCLE
        BCElem(iElem)%nInnerSides = BCElem(iElem)%nInnerSides+1
      END DO ! ilocSide=1,6
      BCElem(iElem)%lastSide=BCElem(iElem)%nInnerSides
      ! loop over all sides, to reduce required storage, if a side is marked once,
      ! it has not be checked for further sides
      SideIndex=0
      DO ilocSide=1,6
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        BCSideID2=SideID
        IF(SideID.GT.0) BCSideID2=PartBCSideList(SideID)
        IF (BCSideID2.GT.0) THEN
          xNodes(:,:,:)=BezierControlPoints3D(:,:,:,PartBCSideList(SideID))
          SELECT CASE(ilocSide)
          CASE(XI_MINUS,XI_PLUS)
            firstBezierPoint=0
            lastBezierPoint=NGeo
          CASE DEFAULT
            firstBezierPoint=1
            lastBezierPoint=NGeo-1
          END SELECT
        ELSE
          SELECT CASE(ilocSide)
          CASE(XI_MINUS)
            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,iElem),xNodes(:,:,:))
            firstBezierPoint=0
            lastBezierPoint=NGeo
          CASE(XI_PLUS)
            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,iElem),xNodes(:,:,:))
            firstBezierPoint=0
            lastBezierPoint=NGeo
          CASE(ETA_MINUS)
            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,iElem),xNodes(:,:,:))
            firstBezierPoint=1
            lastBezierPoint=NGeo-1
          CASE(ETA_PLUS)
            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,iElem),xNodes(:,:,:))
            firstBezierPoint=1
            lastBezierPoint=NGeo-1
          CASE(ZETA_MINUS)
            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,iElem),xNodes(:,:,:))
            firstBezierPoint=1
            lastBezierPoint=NGeo-1
          CASE(ZETA_PLUS)
            CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,iElem),xNodes(:,:,:))
            firstBezierPoint=1
            lastBezierPoint=NGeo-1
          END SELECT
        END IF
        DO iSide=1,nTotalSides
          ! only bc sides
          BCSideID  =PartBCSideList(iSide)
          IF(BCSideID.EQ.-1) CYCLE
          ! ignore sides of the same element
          IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE
          IF(SideIndex(iSide).EQ.0)THEN
            leave=.FALSE.
            nTest=1
            DO iTest=1,nTest
              Vec1=0.
              ! all points of bc side
              DO q=firstBezierPoint,lastBezierPoint
                DO p=firstBezierPoint,lastBezierPoint
                  NodeX(:) = BezierControlPoints3D(:,p,q,BCSideID)+Vec1
                  !all nodes of current side
                  DO s=firstBezierPoint,lastBezierPoint
                    DO r=firstBezierPoint,lastBezierPoint
                      dX=ABS(xNodes(1,r,s)-NodeX(1))
                      IF(dX.GT.halo_eps) CYCLE
                      dY=ABS(xNodes(2,r,s)-NodeX(2))
                      IF(dY.GT.halo_eps) CYCLE
                      dZ=ABS(xNodes(3,r,s)-NodeX(3))
                      IF(dZ.GT.halo_eps) CYCLE
                      IF(SQRT(dX*dX+dY*dY+dZ*dZ).LE.halo_eps)THEN
                        IF(SideIndex(iSide).EQ.0)THEN
                          BCElem(iElem)%lastSide=BCElem(iElem)%lastSide+1
                          SideIndex(iSide)=BCElem(iElem)%lastSide
                          leave=.TRUE.
                          EXIT
                        END IF
                      END IF
                    END DO ! r
                    IF(leave) EXIT
                  END DO ! s
                  IF(leave) EXIT
                END DO ! p
                IF(leave) EXIT
              END DO ! q
              IF(leave) EXIT
            END DO ! iTest=1,nTest
          END IF ! SideIndex(iSide).EQ.0
        END DO ! iSide=1,nTotalSides
      END DO ! ilocSide=1,6
      IF(BCElem(iElem)%lastSide.EQ.0) CYCLE
      ! set true, only required for elements without an own bc side
      IsTracingBCElem(iElem)=.TRUE.
      ! allocate complete side list
      ALLOCATE( BCElem(iElem)%BCSideID(BCElem(iElem)%lastSide) )
      ! 1) inner sides
      nSideCount=0
      IF(BCElem(iElem)%nInnerSides.GT.0)THEN
        DO ilocSide=1,6
          SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
          IF(SideID.LE.0) CYCLE
          BCSideID=PartBCSideList(SideID)
          IF(BCSideID.LE.0) CYCLE
          nSideCount=nSideCount+1
          BCElem(iElem)%BCSideID(nSideCount)=SideID
        END DO ! ilocSide
      END IF ! nInnerSides.GT.0
      ! 2) outer sides
      DO iSide=1,nTotalSides
        IF(SideIndex(iSide).GT.0)THEN
          nSideCount=nSideCount+1
          BCElem(iElem)%BCSideID(nSideCount)=iSide !iSide
        END IF
      END DO  ! iSide
    END DO ! iElem=1,nTotalElems
  END IF ! fullMesh
ELSE ! .NOT.DoRefMapping
  ! tracing
  ! mark only elements with bc-side
  nTotalBCElems=0
  DO iElem=1,nTotalElems
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF (SideID.LE.0) CYCLE
      IF(SideID.LE.nBCSides)THEN ! non-halo elements
        IF(.NOT.IsTracingBCElem(iElem))THEN
          IsTracingBCElem(iElem)=.TRUE.
          nTotalBCElems=nTotalBCElems+1
        END IF ! count only single
      END IF
#ifdef MPI
      IF(SideID.GT.nSides)THEN ! halo elements
        IF(BC(SideID).NE.0)THEN
          IF(.NOT.IsTracingBCElem(iElem))THEN
            IsTracingBCElem(iElem)=.TRUE.
            nTotalBCElems=nTotalBCElems+1
          END IF ! count only single
        END IF
      END IF ! SideID.GT.nSides
#endif
    END DO ! ilocSide
  END DO ! iElem
END IF

IF(DoRefMapping)THEN
  DO iElem=1,nTotalElems
    TracingBCInnerSides(iElem) = BCElem(iElem)%nInnerSides
    TracingBCTotalSides(iElem) = BCElem(iElem)%lastSide
  END DO ! iElem
  CALL AddToElemData(ElementOut,'TracingBCInnerSides',IntArray=TracingBCInnerSides(1:nElems))
  CALL AddToElemData(ElementOut,'TracingBCTotalSides',IntArray=TracingBCTotalSides(1:nElems))
END IF

CALL AddToElemData(ElementOut,'IsTracingBCElem'    ,LogArray=IsTracingBCElem(    1:nElems))

! finally, build epsonecell per element
IF(DoRefMapping)THEN
  ALLOCATE(epsOneCell(1:nTotalElems))
ELSE
  ALLOCATE(epsOneCell(1:PP_nElems))
END IF
epsOneCell=0.

nLoop=nTotalElems
IF(.NOT.DoRefMapping) nLoop=PP_nElems
maxScaleJ=0.
DO iElem=1,PP_nElems
  scaleJ=MAXVAL(sJ(:,:,:,iElem))/MINVAL(sJ(:,:,:,iElem))
  epsOneCell(iElem)=1.0+SQRT(3.0*scaleJ*RefMappingEps)
  maxScaleJ=MAX(scaleJ,maxScaleJ)
END DO ! iElem=1,nLoop
DO iElem=PP_nElems+1,nLoop
  epsOneCell(iElem)=1.0+SQRT(maxScaleJ*RefMappingEps)
END DO ! iElem=1,nLoop
CALL AddToElemData(ElementOut,'epsOneCell',RealArray=epsOneCell(1:nElems))

END SUBROUTINE GetBCElemMap


SUBROUTINE CalcElemAndSideNum()
!===================================================================================================================================
! calculate number of different elements types and side types for local and halo region
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_Vars,             ONLY:DoRefMapping
USE MOD_Mesh_Vars,                          ONLY:CurvedElem,nGlobalElems
USE MOD_Particle_Surfaces_Vars,             ONLY:SideType
USE MOD_Particle_Mesh_Vars,                 ONLY:nTotalSides,IsTracingBCElem,nTotalElems
USE MOD_Particle_Mesh_Vars,                 ONLY:nPartSides
USE MOD_Particle_Mesh_Vars,                 ONLY:nTotalBCSides
USE MOD_Particle_MPI_Vars,                  ONLY:PartMPI
#ifdef MPI
USE MOD_Particle_MPI_HALO,                  ONLY:WriteParticlePartitionInformation
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: iElem
INTEGER                                  :: iSide, nDummy
INTEGER                                  :: nBCElems,nBCelemsTot
INTEGER                                  :: nPlanarRectangular, nPlanarNonRectangular,nPlanarCurved,nBilinear,nCurved
INTEGER                                  :: nPlanarRectangularTot, nPlanarNonRectangularTot,nPlanarCurvedTot,nBilinearTot,nCurvedTot
INTEGER                                  :: nLinearElems, nCurvedElems, nCurvedElemsTot
#ifdef MPI
INTEGER                                  :: nPlanarRectangularHalo, nPlanarNonRectangularHalo,nPlanarCurvedHalo, &
                                            nBilinearHalo,nCurvedHalo,nCurvedElemsHalo,nLinearElemsHalo,nBCElemsHalo
#endif /*MPI*/
INTEGER                                  :: nLoop
!===================================================================================================================================

! zero counter for side and elem types
nPlanarRectangular         = 0
nPlanarNonRectangular      = 0
nPlanarCurved              = 0
nBilinear                  = 0
nCurved                    = 0
nBCElems                   = 0
nCurvedElems               = 0
nLinearElems               = 0
#ifdef MPI
nPlanarRectangularHalo     = 0
nPlanarNonRectangularHalo  = 0
nPlanarCurvedHalo          = 0
nBilinearHalo              = 0
nCurvedHalo                = 0
nCurvedElemsHalo           = 0
nLinearElemsHalo           = 0
nBCElemsHalo               = 0
#endif /*MPI*/

DO iElem=1,nTotalElems
  ! count elements by type and in own and halo region
  IF(iElem.LE.PP_nElems)THEN
    IF(CurvedElem(iElem))THEN
      nCurvedElems=nCurvedElems+1
    ELSE
      nLinearElems=nLinearElems+1
    END IF
    IF(IsTracingBCElem(iElem))THEN
      nBCElems=nBCElems+1
    END IF ! count only single
#ifdef MPI
  ELSE
    IF(CurvedElem(iElem)) THEN
      nCurvedElemsHalo=nCurvedElemsHalo+1
    ELSE
      nLinearElemsHalo=nLinearElemsHalo+1
    END IF
    IF(IsTracingBCElem(iElem))THEN
      nBCElemsHalo=nBCElemsHalo+1
    END IF ! count only single
#endif /*MPI*/
  END IF
END DO
nLoop = nTotalSides
IF (DoRefMapping) nLoop = nTotalBCSides
DO iSide=1,nLoop
  IF (iSide.LE.nPartSides) THEN
    SELECT CASE(SideType(iSide))
    CASE (PLANAR_RECT)
      nPlanarRectangular=nPlanarRectangular+1
    CASE (PLANAR_NONRECT)
      nPlanarNonRectangular=nPlanarNonRectangular+1
    CASE (BILINEAR)
      nBilinear = nBilinear+1
    CASE (PLANAR_CURVED)
      nPlanarCurved = nPlanarCurved+1
    CASE (CURVED)
      nCurved = nCurved+1
    END SELECT
#ifdef MPI
  ELSE IF (iSide.GT.nPartSides) THEN
    SELECT CASE(SideType(iSide))
    CASE (PLANAR_RECT)
      nPlanarRectangularHalo=nPlanarRectangularHalo+1
    CASE (PLANAR_NONRECT)
      nPlanarNonRectangularHalo=nPlanarNonRectangularHalo+1
    CASE (BILINEAR)
      nBilinearHalo = nBilinearHalo+1
    CASE (PLANAR_CURVED)
      nPlanarCurvedHalo = nPlanarCurvedHalo+1
    CASE (CURVED)
      nCurvedHalo = nCurvedHalo+1
    END SELECT
#endif /*MPI*/
  END IF
END DO

#ifdef MPI
IF(MPIRoot) THEN
  CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nCurvedElems,nCurvedElemsTot,1,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(DoRefMapping) CALL MPI_REDUCE(nBCElems,nBCElemsTot ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
ELSE ! no Root
  CALL MPI_REDUCE(nPlanarRectangular     ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nPlanarNonRectangular  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nBilinear              ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nPlanarCurved          ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nCurved                ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nCurvedElems           ,nDummy,1,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(DoRefMapping) CALL MPI_REDUCE(nBCElems  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
END IF
#else
nPlanarRectangularTot   =nPlanarRectangular
nPlanarNonRectangularTot=nPlanarNonRectangular
nBilinearTot            =nBilinear
nPlanarCurvedTot        =nPlanarCurved
nCurvedTot              =nCurved
nCurvedElemsTot         =nCurvedElems
IF(DorefMapping) nBCElemstot=nBCElems
#endif /*MPI*/

SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar-rectangular     faces: ', nPlanarRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar-non-rectangular faces: ', nPlanarNonRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of bi-linear              faces: ', nBilineartot
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar-curved          faces: ', nPlanarCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of curved                 faces: ', nCurvedtot
! and add number of curved elems
IF(DoRefMapping)THEN
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of BC-adjoined            elems: ', nBCElemstot
END IF
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of (bi-)linear            elems: ', nGlobalElems-nCurvedElemsTot
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of curved                 elems: ', nCurvedElemsTot
SWRITE(UNIT_StdOut,'(132("-"))')
#ifdef MPI
CALL WriteParticlePartitionInformation(nPlanarRectangular+nPlanarNonRectangular,nBilinear,nCurved+nPlanarCurved,                    &
                                       nPlanarRectangularHalo+nPlanarNonRectangularHalo,nBilinearHalo,nCurvedHalo+nPlanarCurvedHalo &
                                      ,nBCElems,nLinearElems,nCurvedElems,nBCElemsHalo,nLinearElemsHalo,nCurvedElemsHalo)
#endif

END SUBROUTINE CalcElemAndSideNum


SUBROUTINE GetLinearSideBaseVectors()
!===================================================================================================================================
! computes the face base vector for linear (planar or bilinear) face intersection calculation
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_Vars,        ONLY:DoRefMapping
USE MOD_Mesh_Vars,                     ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,        ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,        ONLY:BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,BaseVectorsScale
!USE MOD_Particle_Surfaces_Vars,        ONLY:BaseVectors0flip,BaseVectors1flip,BaseVectors2flip,BaseVectors3flip
! USE MOD_Particle_Surfaces_Vars,        ONLY:SideID2PlanarSideID
! USE MOD_Particle_Surfaces_Vars,        ONLY:SideType
USE MOD_Particle_Mesh_Vars,            ONLY:nTotalSides,nTotalBCSides
USE MOD_Particle_Mesh_Vars,            ONLY:PartBCSideList
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: iSide, BCSide
REAL                                  :: crossVec(3)
! INTEGER                               :: iSide_temp
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' GET LINEAR SIDE BASEVECTORS...'
IF(.NOT.DoRefMapping)THEN
!   ALLOCATE(SideID2PlanarSideID(1:nTotalSides))
!   SideID2PlanarSideID(:) = 0
!   iSide_temp = 0
!   DO iSide=1,nTotalSides
!     IF (SideType(iSide).EQ.PLANAR_RECT) THEN
!       iSide_temp = iSide_temp + 1
!       SideID2PlanarSideID(iSide) = iSide_temp
!     END IF
!   END DO
  
  ALLOCATE( BaseVectors0(1:3,1:nTotalSides),&
            BaseVectors1(1:3,1:nTotalSides),&
            BaseVectors2(1:3,1:nTotalSides),&
            BaseVectors3(1:3,1:nTotalSides),&
            BaseVectorsScale(1:nTotalSides))
  !IF (GEO%nPeriodicVectors.GT.0) THEN
  !  ALLOCATE( BaseVectors0flip(1:3,1:nTotalSides),&
  !            BaseVectors1flip(1:3,1:nTotalSides),&
  !            BaseVectors2flip(1:3,1:nTotalSides),&
  !            BaseVectors3flip(1:3,1:nTotalSides))
  !END IF
   
  DO iSide=1,nTotalSides
    ! extension for periodic sides
!     IF ((SideType(iSide).EQ.PLANAR_RECT) &
!        .OR. (SideType(iSide).EQ.PLANAR_NONRECT) .OR. (SideType(iSide).EQ.BILINEAR)))THEN
!       iSide_temp = SideID2PlanarSideID(iSide)
    BaseVectors0(:,iSide) = (+BezierControlPoints3D(:,0,0,iSide)+BezierControlPoints3D(:,NGeo,0,iSide)   &
                              +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3D(:,0,0,iSide)+BezierControlPoints3D(:,NGeo,0,iSide)   &
                              -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3D(:,0,0,iSide)-BezierControlPoints3D(:,NGeo,0,iSide)   &
                              +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3D(:,0,0,iSide)-BezierControlPoints3D(:,NGeo,0,iSide)   &
                              -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
ELSE
  ALLOCATE( BaseVectors0(1:3,1:nTotalBCSides),&
            BaseVectors1(1:3,1:nTotalBCSides),&
            BaseVectors2(1:3,1:nTotalBCSides),&
            BaseVectors3(1:3,1:nTotalBCSides),&
            BaseVectorsScale(1:nTotalBCSides))
  DO iSide=1,nTotalSides
    BCSide = PartBCSideList(iSide)
    ! extension for periodic sides
    IF(BCSide.EQ.-1) CYCLE
    BaseVectors0(:,BCSide) = (+BezierControlPoints3D(:,0,0,BCSide)+BezierControlPoints3D(:,NGeo,0,BCSide)   &
                              +BezierControlPoints3D(:,0,NGeo,BCSide)+BezierControlPoints3D(:,NGeo,NGeo,BCSide) )
    BaseVectors1(:,BCSide) = (-BezierControlPoints3D(:,0,0,BCSide)+BezierControlPoints3D(:,NGeo,0,BCSide)   &
                              -BezierControlPoints3D(:,0,NGeo,BCSide)+BezierControlPoints3D(:,NGeo,NGeo,BCSide) )
    BaseVectors2(:,BCSide) = (-BezierControlPoints3D(:,0,0,BCSide)-BezierControlPoints3D(:,NGeo,0,BCSide)   &
                              +BezierControlPoints3D(:,0,NGeo,BCSide)+BezierControlPoints3D(:,NGeo,NGeo,BCSide) )
    BaseVectors3(:,BCSide) = (+BezierControlPoints3D(:,0,0,BCSide)-BezierControlPoints3D(:,NGeo,0,BCSide)   &
                              -BezierControlPoints3D(:,0,NGeo,BCSide)+BezierControlPoints3D(:,NGeo,NGeo,BCSide) )
    crossVec = CROSS(BaseVectors1(:,BCSide),BaseVectors2(:,BCSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(BCSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
END IF

SWRITE(UNIT_stdOut,'(A)')' GET LINEAR SIDE BASEVECTORS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE GetLinearSideBaseVectors


SUBROUTINE ElemConnectivity() 
!===================================================================================================================================
! computes the element connectivity between different elements, inclusive the halo region
! and mortar interfaces
! CAUTION: the assumption is, that one element is only linked once or twice with another element
!          one link: normal inner connection or periodic connection
!          two links: one normal connection PLUS one periodic connection
!          more than 2 links: funny.
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,  ONLY:PartElemToElemGlob, PartElemToElemAndSide,nTotalElems,PartElemToSide,PartBCSideList &
                                 ,SidePeriodicType,ElemToGlobalElemID
USE MOD_Particle_MPI_Vars,   ONLY:PartHaloElemToProc
USE MOD_Mesh_Vars,           ONLY:OffSetElem,BC,BoundaryType,MortarType
USE MOD_Particle_Surfaces_Vars, ONLY:SideNormVec
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
#ifdef MPI
USE MOD_MPI_Vars,            ONLY:OffSetElemMPI
#endif /*MPI*/
USE MOD_Mesh_vars
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,ilocSide,iMortar,ProcID,ilocSide2,iMortar2,NbElemID,ElemID,BCID,SideID,BCSideID
INTEGER(KIND=8)               :: GlobalElemID
LOGICAL                       :: found
REAL                          :: Vec1(1:3)
#ifdef MPI
INTEGER                       :: iHaloElem
INTEGER(KIND=8)               :: HaloGlobalElemID
#endif /*MPI*/
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' BUILD MESH-CONNECTIVITY ... '

SDEALLOCATE(PartElemToElemAndSide)
ALLOCATE(PartElemToElemAndSide(1:8,1:6,1:nTotalElems))
                      ! [1]1:4 - MortarNeighborElemID
                      ! [1]5:8 -       Neighbor locSideID
                      ! [2]1:6 - locSideID
                      ! [3]    - nTotalElems 
                      ! if the connections points to an element which is not in MY region (MY elems + halo elems)
                      ! then this connection points to -1
ALLOCATE(ElemToGlobalElemID(1:nTotalElems))
! nullify
PartElemToElemAndSide=-1
ElemToGlobalElemID=-1

! now, map the PartElemToElemGlob to local elemids
! loop over all Elems and map the neighbor element to local coordinates
DO iElem=1,nTotalElems
  IF(iElem.LE.nElems)THEN
    ElemToGlobalElemID(iElem)=offSetElem+iElem
#ifdef MPI
  ELSE
    ProcID=PartHaloElemToProc(NATIVE_PROC_ID,iElem)
    ElemToGlobalElemID(iElem)=offSetElemMPI(ProcID) + PartHaloElemToProc(NATIVE_ELEM_ID,iElem)
#endif /*MPI*/
  END IF
  DO ilocSide=1,6
    DO iMortar=1,4
      GlobalElemID=PartElemToElemGlob(iMortar,ilocSide,iElem)
      IF(GlobalElemID.LE.0) CYCLE
      ! check if the element is in MY range of elements
      IF((GlobalElemID.GE.OffSetElem+1).AND.(GlobalElemID.LE.(OffSetElem+PP_nElems)))THEN
        PartElemToElemAndSide(iMortar,ilocSide,iElem)=INT(GlobalElemID-OffSetElem,4)
        CYCLE
      END IF
#ifdef MPI
      ! neighbor element not found, hence, it can be a halo element
      DO iHaloElem=PP_nElems+1,nTotalElems
        ProcID=PartHaloElemToProc(NATIVE_PROC_ID,iHaloElem)
        HaloGlobalElemID=offSetElemMPI(ProcID) + PartHaloElemToProc(NATIVE_ELEM_ID,iHaloElem)
        CHECKSAFEINT(HaloGlobalElemID,4)
        IF(HaloGlobalElemID.EQ.GlobalElemID)THEN
          PartElemToElemAndSide(iMortar,ilocSide,iElem)=iHaloElem
          EXIT
        END IF
      END DO ! iHaloElem=1,nTotalElems
#endif /*MPI*/
    END DO ! iMortar=1,4
  END DO ! ilocSide=1,6
END DO ! iElem=1,PP_nElems

! which local side of neighbor element is connected to MY element
DO iElem=1,nTotalElems
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)    
    ! check for ref-mapping or tracing
    IF(DoRefMapping)THEN
      IF(SideID.GT.0)THEN
        BCSideID=PartBCSideList(SideID)
      ELSE
        BCSideID=-1
      END IF
    ELSE
      BCSideID=SideID
    END IF
    ! disable BCSideID, if it is NOT a periodic side
    IF(BCSideID.GT.0)THEN ! only BC faces 
      IF(SidePeriodicType(SideID).NE.0)THEN ! only periodic sides
        Vec1=SideNormVec(1:3,BCSideID)
        IF(ALMOSTZERO(DOT_PRODUCT(Vec1,Vec1))) CALL abort(&
__STAMP__&
        , ' Error in ElemConnectivity. No SideNormVec!',iElem,REAL(ilocSide))
      ELSE ! disable non-periodic  sides
        Vec1=0.
        BCSideID=-1
      END IF
    END IF
    IF(BCSideID.GT.0)THEN ! periodic sides
      DO iMortar=1,4
        NBElemID=PartElemToElemAndSide(iMortar,ilocSide,iElem)
        IF(NBElemID.EQ.-1) CYCLE
        found=.FALSE.
        ! loop  over all local sides of neighbor element to find the right face
        DO ilocSide2=1,6
          DO iMortar2=1,4
            ElemID=PartElemToElemAndSide(iMortar2,ilocSide2,NBElemID)
            IF(ElemID.LE.0) CYCLE
            IF(ElemID.EQ.iElem) THEN
              ! check if periodic side
              SideID=PartElemToSide(E2S_SIDE_ID,ilocSide2,NBElemID)    
              ! check for ref-mapping or tracing
              IF(DoRefMapping)THEN
                IF(SideID.GT.0)THEN
                  BCSideID=PartBCSideList(SideID)
                ELSE
                  BCSideID=-1
                END IF
              ELSE
                BCSideID=SideID
              END IF
              IF(BCSideID.GT.0)THEN ! only BC faces 
                IF(SidePeriodicType(SideID).NE.0)THEN ! only periodic sides
                  IF(ALMOSTEQUAL(ABS(DOT_PRODUCT(Vec1,SideNormVec(1:3,BCSideID))),1.0))THEN
                    ! finally, found matching local sides
                    PartElemToElemAndSide(iMortar+4,ilocSide,iElem)=ilocSide2
                    Found=.TRUE.
                    EXIT
                  ELSE
                    CYCLE
                  END IF
                ELSE ! disable non-periodic  sides
                  CYCLE
                END IF
              ELSE
                CYCLE
              END IF
            END IF
          END DO ! iMortar=1,4
          IF(Found) EXIT
        END DO ! ilocSide=1,6
      END DO ! iMortar=1,4
    ELSE ! non-periodic sides
      DO iMortar=1,4
        NBElemID=PartElemToElemAndSide(iMortar,ilocSide,iElem)
        IF(NBElemID.EQ.-1) CYCLE
        found=.FALSE.
        ! loop  over all local sides of neighbor element to find the right face
        DO ilocSide2=1,6
          DO iMortar2=1,4
            ElemID=PartElemToElemAndSide(iMortar2,ilocSide2,NBElemID)
            IF(ElemID.LE.0) CYCLE
            IF(ElemID.EQ.iElem) THEN
              ! check if periodic side
              SideID=PartElemToSide(E2S_SIDE_ID,ilocSide2,NBElemID)    
              ! check for ref-mapping or tracing
              IF(DoRefMapping)THEN
                IF(SideID.GT.0)THEN
                  BCSideID=PartBCSideList(SideID)
                ELSE
                  BCSideID=-1
                END IF
              ELSE
                BCSideID=SideID
              END IF
              IF(BCSideID.GT.0)THEN ! BC face?
                IF(SidePeriodicType(SideID).NE.0)THEN ! only non-periodic sides
                  CYCLE
                ELSE ! enable non-periodic  sides
                  ! finally, found matching local sides
                  PartElemToElemAndSide(iMortar+4,ilocSide,iElem)=ilocSide2
                  Found=.TRUE.
                  EXIT
                END IF
              ELSE
                ! finally, found matching local sides
                PartElemToElemAndSide(iMortar+4,ilocSide,iElem)=ilocSide2
                Found=.TRUE.
                EXIT
              END IF
            END IF
          END DO ! iMortar=1,4
          IF(Found) EXIT
        END DO ! ilocSide=1,6
      END DO ! iMortar=1,4
    END IF ! periodic sides
  END DO ! ilocSide=1,6
END DO ! iElem=1,PP_nElems

! sanity check
DO iElem=1,nTotalElems
  DO ilocSide=1,6
    DO iMortar=1,4
      IF((PartElemToElemAndSide(iMortar,ilocSide,iElem).GT.0).AND.(PartElemToElemAndSide(iMortar+4,ilocSide,iElem).EQ.-1))THEN
        IPWRITE(UNIT_StdOut,*) ' iElem:     ', iElem
        IPWRITE(UNIT_StdOut,*) ' ilocSide:  ', ilocSide
        IPWRITE(UNIT_StdOut,*) ' NBElem-ID: ', PartElemToElemAndSide(iMortar,ilocSide,iElem)
        IPWRITE(UNIT_StdOut,*) ' NBElem-ID: ', PartElemToElemAndSide(iMortar+4,ilocSide,iElem)
        CALL abort(&
__STAMP__&
        , ' Error in ElemConnectivity. Found no neighbor locSideID. iElem,ilocSide',iElem,REAL(ilocSide))
      END IF
      IF((PartElemToElemAndSide(iMortar,ilocSide,iElem).EQ.-1).AND.(PartElemToElemAndSide(iMortar+4,ilocSide,iElem).GT.-1))THEN
        IPWRITE(UNIT_StdOut,*) ' iElem:     ', iElem
        IPWRITE(UNIT_StdOut,*) ' ilocSide:  ', ilocSide
        IPWRITE(UNIT_StdOut,*) ' NBElem-ID: ', PartElemToElemAndSide(iMortar,ilocSide,iElem)
        CALL abort(&
__STAMP__&
        , ' Error in ElemConnectivity. Found no neighbor ElemID. iElem,ilocSide',iElem,REAL(ilocSide))
      END IF
    END DO ! iMortar=1,4
  END DO ! ilocSide=1,6
END DO

! check is working on CONFORM mesh!!!
DO iElem=1,nTotalElems
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)    
    IF(DoRefMapping)THEN
      IF(SideID.LT.1) CYCLE
    ELSE
      IF(SideID.LE.0) CALL abort(&
__STAMP__&
       , ' Error in PartElemToSide! No SideID for side!. iElem,ilocSide',iElem,REAL(ilocSide))
    END IF
    IF(MortarType(1,SideID).NE.0) CYCLE
    BCID=BC(SideID)
    IF(BCID.NE.0)THEN
      IF(BoundaryType(BCID,BC_TYPE).GT.1) CYCLE
    END IF
    IF(PartElemToElemAndSide(1,ilocSide,iElem).LT.1)THEN
       CALL abort(&
__STAMP__&
      , ' Error in ElemConnectivity. Found no neighbor ElemID. iElem,ilocSide',iElem,REAL(ilocSide))
      END IF
  END DO ! ilocSide=1,6
END DO

#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iERROR)
#endif
SWRITE(UNIT_stdOut,'(A)')' BUILD MESH-CONNECTIVITY SUCCESSFUL '
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE ElemConnectivity


SUBROUTINE NodeNeighbourhood()
!===================================================================================================================================
!> Subroutine for initialization of neighbourhood with nodes using GEO container
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: nElems, nNodes
USE MOD_Particle_Mesh_Vars ,ONLY: GEO, PartElemToElemAndSide
#ifdef CODE_ANALYZE
#ifdef MPI
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_Particle_MPI_Vars  ,ONLY: PartHaloElemToProc
USE MOD_MPI_Vars           ,ONLY: offsetElemMPI
#endif /*MPI*/
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE tNodeToElem
  INTEGER :: ElemID(50)
END TYPE
TYPE(tNodeToElem)      :: TempNodeToElem(1:nNodes)
INTEGER                :: TempElemsOnNode(1:nNodes)
INTEGER                :: Element, iLocSide, k, l
LOGICAL                :: ElemExists
INTEGER                :: iElem, jNode
INTEGER                :: iNode
INTEGER                :: TempHaloElems(1:500)
INTEGER                :: TempHaloNumElems
#ifdef MPI
LOGICAL                :: HaloNeighNode(1:nNodes)
#endif /*MPI*/
LOGICAL                :: ElemDone
REAL                   :: MPINodeCoord(3), ElemCoord(3)
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' BUILD NODE-NEIGHBOURHOOD ... '

#ifdef MPI
! set nodes of sides with halo element connected to it as HaloNeighNodes
HaloNeighNode(:) = .FALSE.
DO iElem=1,nElems
  DO iLocSide = 1,6
    IF (PartElemToElemAndSide(1,iLocSide,iElem).GT.PP_nElems) THEN
      DO iNode = 1,4
        HaloNeighNode(GEO%ElemSideNodeID(iNode,iLocSide,iElem)) = .TRUE.
      END DO
    END IF
  END DO
END DO
#endif /*MPI*/

#ifdef CODE_ANALYZE
DO iNode=1,nNodes
  IF (HaloNeighNode(iNode)) THEN
    print*,'Rank: ',MyRank,'---- local Node: ',iNode,' is halo node'
  END IF
END DO
#endif /*CODE_ANALYZE*/

ALLOCATE(GEO%NumNeighborElems(1:PP_nElems))
ALLOCATE(GEO%ElemToNeighElems(1:PP_nElems))
GEO%NumNeighborElems(:)=0
TempElemsOnNode(:)=0
DO iNode=1,nNodes
  TempNodeToElem(iNode)%ElemID=-1
END DO

! find all real neighbour elements for elements with halo neighbours 
! recursively checking connected halo area
DO iElem =1, PP_nElems
  TempHaloElems(1:500) = 0
  TempHaloNumElems = 0
  ! now check every side for neighbours, add valid neighbour to corresponding array and proceed recursively until neighbourhood
  ! is finished
  DO iLocSide = 1, 6
    ElemExists = .FALSE.
    Element = PartElemToElemAndSide(1,iLocSide,iElem)
    IF (Element.GT.0) THEN !side has neighbour element
      DO l=1, TempHaloNumElems
        IF(Element.EQ.TempHaloElems(l)) THEN
          ElemExists=.TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT.ElemExists) THEN
        TempHaloNumElems = TempHaloNumElems + 1
        TempHaloElems(TempHaloNumElems) = Element
      END IF
      CALL RecurseCheckNeighElems(iElem,Element,TempHaloNumElems,TempHaloElems)
    END IF
  END DO
  IF (TempHaloNumElems.LE.0) CALL abort(&
__STAMP__&
,'ERROR in FindNeighbourElems! no neighbour elements found for Element',iElem)
  ! write local variables into global array
  GEO%NumNeighborElems(iElem) = TempHaloNumElems
  ALLOCATE(GEO%ElemToNeighElems(iElem)%ElemID(1:GEO%NumNeighborElems(iElem)))
  GEO%ElemToNeighElems(iElem)%ElemID(1:GEO%NumNeighborElems(iElem)) = TempHaloElems(1:TempHaloNumElems)
  ! add neighbour elements to respective nodes of iElem if nodes of neighbour are the same to those of iElem
  DO l=1,GEO%NumNeighborElems(iElem)
    Element = GEO%ElemToNeighElems(iElem)%ElemID(l)
    DO iNode = 1, 8
      ElemExists = .FALSE.
      DO k = 1,TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem))
        IF (TempNodeToElem(GEO%ElemToNodeID(iNode,iElem))%ElemID(k).EQ.Element) THEN
          ElemExists = .TRUE.
          EXIT
        END IF
      END DO
      IF (ElemExists) CYCLE
      ElemCoord(1:3) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
      DO jNode = 1, 8
        ElemDone = .FALSE.
        MPINodeCoord(1:3) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(jNode,Element))
        IF(ALMOSTEQUAL(MPINodeCoord(1),ElemCoord(1)).AND.ALMOSTEQUAL(MPINodeCoord(2),ElemCoord(2)) &
            .AND.ALMOSTEQUAL(MPINodeCoord(3),ElemCoord(3))) THEN
          TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem)) = TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem)) + 1
          TempNodeToElem(GEO%ElemToNodeID(iNode,iElem))%ElemID(TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem))) = Element
          ElemDone = .TRUE.
        END IF
        IF (ElemDone) EXIT
      END DO ! jNode=1,8
      IF (ElemDone) CYCLE
    END DO ! iNode=1,8
  END DO ! l=1,NumNeihborElems
END DO ! iElem=1,PP_nElems

! check if current element already added to every node of the element and add to missing (elements for corner nodes not added yet)
! this part needs its own loop over all elements
DO iElem=1,PP_nElems
  DO iNode = 1,8
    ElemExists = .FALSE.
    DO l = 1,TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem))
      IF (TempNodeToElem(GEO%ElemToNodeID(iNode,iElem))%ElemID(l).EQ.iElem) THEN
        ElemExists = .TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT.ElemExists) THEN
      TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem)) = TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem)) + 1
      TempNodeToElem(GEO%ElemToNodeID(iNode,iElem))%ElemID(TempElemsOnNode(GEO%ElemToNodeID(iNode,iElem))) = Element
    END IF
  END DO
END DO

! write number of elements for corresponding proc global nodes into GEO structure
ALLOCATE(GEO%ElemsOnNode(1:nNodes))
ALLOCATE(GEO%NodeToElem(1:nNodes))
DO iNode=1,nNodes
  GEO%ElemsOnNode(iNode) = TempElemsOnNode(iNode)
  ALLOCATE(GEO%NodeToElem(iNode)%ElemID(1:GEO%ElemsOnNode(iNode)))
  GEO%NodeToElem(iNode)%ElemID(1:GEO%ElemsOnNode(iNode)) = TempNodeToElem(iNode)%ElemID(1:TempElemsOnNode(iNode))
END DO

! fill array of neighbour nodes to proc global nodes
ALLOCATE(GEO%NeighNodesOnNode(1:nNodes))
ALLOCATE(GEO%NodeToNeighNode(1:nNodes))
DO iNode=1,nNodes
  TempHaloElems(:) = 0
  TempHaloNumElems = 0
  DO iElem=1,GEO%ElemsOnNode(iNode)
    DO jNode=1,8
      IF (GEO%ElemToNodeID(jNode,GEO%NodeToElem(iNode)%ElemID(iElem)).EQ.iNode) CYCLE
      ElemExists=.false.
      DO l=1, TempHaloNumElems
        IF(GEO%ElemToNodeID(jNode,GEO%NodeToElem(iNode)%ElemID(iElem)).EQ.TempHaloElems(l)) THEN
          ElemExists=.true.
          CYCLE
        END IF
      END DO
      IF(.NOT.ElemExists) THEN
        TempHaloNumElems = TempHaloNumElems + 1
        TempHaloElems(TempHaloNumElems) = GEO%ElemToNodeID(jNode,GEO%NodeToElem(iNode)%ElemID(iElem))
      END IF
    END DO
  END DO
  ! write local variables into global array
  ALLOCATE(GEO%NodeToNeighNode(iNode)%ElemID(1:TempHaloNumElems))
  GEO%NeighNodesOnNode(iNode)=TempHaloNumElems
  GEO%NodeToNeighNode(iNode)%ElemID(1:TempHaloNumElems) = TempHaloElems(1:TempHaloNumElems)
END DO

#ifdef CODE_ANALYZE
! write some code analyze output of connectivity
DO iElem=1,PP_nElems
#ifdef MPI
  print*,'Rank: ',MyRank,'------ Element: ',iElem+offsetElem,' has ',GEO%NumNeighborElems(iElem),' Neighbours'
  print*,'Rank: ',MyRank,'------ Neighbours are:'
  DO l=1,GEO%NumNeighborElems(iElem)
    IF (GEO%ElemToNeighElems(iElem)%ElemID(l).GT.PP_nElems) THEN
      print*,'Rank: ',MyRank,offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,GEO%ElemToNeighElems(iElem)%ElemID(l))) &
          + PartHaloElemToProc(NATIVE_ELEM_ID,GEO%ElemToNeighElems(iElem)%ElemID(l))
    ELSE
      print*,'Rank: ',MyRank,GEO%ElemToNeighElems(iElem)%ElemID(l) + offsetElem
    END IF
  END DO
#else
  print*,'Rank: ',MyRank,'------ Element: ',iElem,' has ',GEO%NumNeighborElems(iElem),' Neighbours'
  print*,'Rank: ',MyRank,'------ Neighbours are:',GEO%ElemToNeighElems(iElem)%ElemID(:)
#endif /*MPI*/
END DO

DO iNode=1,nNodes
  print*,'Rank: ',MyRank,'------ Node: ',iNode,' has: ',GEO%ElemsOnNode(iNode),' Elements'
END DO
#endif /*CODE_ANALYZE*/

SWRITE(UNIT_stdOut,'(A)')' BUILD NODE-NEIGHBOURHOOD SUCCESSFUL '
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE NodeNeighbourhood


RECURSIVE SUBROUTINE RecurseCheckNeighElems(StartElem,HaloElem,TempHaloNumElems,TempHaloElems)
!===================================================================================================================================
!> Subroutine for recursively checking halo neighbourhood for connectivity to current elem
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars ,ONLY: GEO, PartElemToElemAndSide
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)  :: StartElem,HaloElem,TempHaloElems(1:500), TempHaloNumElems
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: iNode, jNode
INTEGER                :: iLocSide, l
INTEGER                :: currentElem
LOGICAL                :: ElemExists, ElemDone
REAL                   :: MPINodeCoord(3), ElemCoord(3)
!===================================================================================================================================
DO iLocSide = 1,6
  ElemExists = .FALSE.
  currentElem = PartElemToElemAndSide(1,iLocSide,HaloElem)
  IF (currentElem.GT.0 .AND. currentElem.NE.StartElem) THEN
  !IF (currentElem.GT.PP_nElems) THEN
    DO l=1, TempHaloNumElems
      IF(currentElem.EQ.TempHaloElems(l)) THEN
        ElemExists=.TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT.ElemExists) THEN
      ElemDone = .FALSE.
      DO iNode = 1, 8
        DO jNode = 1, 8
          MPINodeCoord(1:3) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(jNode,currentElem))
          ElemCoord(1:3) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,StartElem))
          IF(ALMOSTEQUAL(MPINodeCoord(1),ElemCoord(1)).AND.ALMOSTEQUAL(MPINodeCoord(2),ElemCoord(2)) &
              .AND.ALMOSTEQUAL(MPINodeCoord(3),ElemCoord(3))) THEN
            TempHaloNumElems = TempHaloNumElems + 1
            TempHaloElems(TempHaloNumElems) = currentElem
            ElemDone = .TRUE.
            CALL RecurseCheckNeighElems(StartElem,currentElem,TempHaloNumElems,TempHaloElems)
          END IF
          IF (ElemDone) EXIT
        END DO
        IF (ElemDone) EXIT
      END DO
    END IF
  END IF
END DO

END SUBROUTINE RecurseCheckNeighElems


SUBROUTINE DuplicateSlavePeriodicSides() 
!===================================================================================================================================
! duplicate periodic sides from old nPartSides=nSides to nPartSidesNew=nSides+nDuplicatePeriodicSides
! without MPI 
! 1) loop over all sides and detect periodic sides
! 2) increase BezierControlPoints and SideXXX from old nSides to nSides+nDuplicatePeriodicSides
! 3) loop over the OLD sides and copy the corresponding SideXXX. The missing BezierControlPoints (periodic shifted values) 
!    are build from the other element. Now two BezierControlPoints existes which are shifted by the sideperiodicvector
! 4) shift and map sideperiodicvector and displacement to match new sides
! with MPI
! 1) loop over all sides and detect periodic sides
! 2) increase BezierControlPoints and SideXXX from old nSides to nSides+nDuplicatePeriodicSides
! 3) loop over the OLD sides and copy the corresponding SideXXX. The missing BezierControlPoints (periodic shifted values) 
!    are build  the other element. Now two BezierControlPoints existes which are shifted by the sideperiodicvector
! 3b) newSideId depends on localSideID and yourMPISide
!     a) both periodic sides are on proc:
!        * duplicate side and two separate sideids with changes in partsidetoelem
!        * build missing side with own data
!     b) periodic side is MPI Side
!        I) mySide (Master)-Side
!           *  nothing to due, old side can be reused
!        II) yourSide (Slave)-Side
!           *  build new Side with own data
! 4) shift and map sideperiodicvector and displacement to match new sides
! Note:
! periodic sides are unique for the DG operator and duplicated for the particle tracking
! CAUTION:
! Routine has to be called before MarkAllBCSides
!===================================================================================================================================
! MODULES                                                                                                                          !
USE  MOD_GLobals
USE MOD_Mesh_Vars              ,ONLY: MortarType,BC,NGeo,nBCs,nSides,BoundaryType,MortarSlave2MasterInfo,nElems,XCL_NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide,PartSideToElem,nTotalSides,SidePeriodicType,nPartPeriodicSides,GEO &
                                     ,nTotalBCSides,nPartSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Mesh_Vars              ,ONLY: NGeoElevated
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals,RotateMasterToSlave,GetBezierControlPoints3D
USE MOD_Particle_Surfaces_vars ,ONLY: BezierControlPoints3D,SideSlabIntervals,BezierControlPoints3DElevated &
                                        ,SideSlabIntervals,SideSlabNormals,BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars ,ONLY: CartesianPeriodic
USE MOD_Particle_MPI_Vars      ,ONLY: printBezierControlPointsWarnings
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iSide,NBElemID,tmpnSides,NBlocSideID,ElemID,newSideID,locSideID,PVID
INTEGER                              :: BCID,iBC,flip,ilocSide,iElem,SideID,idir
REAL,ALLOCATABLE                     :: DummyBezierControlPoints3D(:,:,:,:)                                
REAL,ALLOCATABLE                     :: DummyBezierControlPoints3DElevated(:,:,:,:)                                
REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: DummySideSlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)      :: DummySideSlabIntervals               ! intervalls beta1, beta2, beta3
!REAL,ALLOCATABLE,DIMENSION(:,:)      :: DummySidePeriodicDisplacement        ! intervalls beta1, beta2, beta3
LOGICAL,ALLOCATABLE,DIMENSION(:)     :: DummyBoundingBoxIsEmpty
INTEGER,ALLOCATABLE,DIMENSION(:)     :: DummyBC
INTEGER,ALLOCATABLE,DIMENSION(:)     :: DummyMortarSlave2MasterInfo
INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: DummyMortarType
INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: DummyPartSideToElem
INTEGER,ALLOCATABLE,DIMENSION(:)     :: DummySidePeriodicType
LOGICAL                              :: MapPeriodicSides
REAL                                 :: MinMax(1:2),MinMaxGlob(1:6),xTest(1:3)
!===================================================================================================================================

! 1) loop over all sides and detect periodic sides
nPartPeriodicSides=0
MapPeriodicSides=.FALSE.
IF(.NOT.CartesianPeriodic .AND. GEO%nPeriodicVectors.GT.0)THEN
  DO iSide=1,nSides
    IF(SidePeriodicType(iSide).NE.0)THEN
      ! abort if particles are traced over mortar sides
      IF(MortarSlave2MasterInfo(iSide).NE.-1.OR.MortarType(1,iSide).GE.0) CALL abort(&
__STAMP__&
      , ' Periodic tracing over mortar sides is not implemented!')
      ! ignore MPI sides, these have NOT to be mirrored
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      IF(ElemID.EQ.-1) THEN
        ! master side is NOT on proc, hence, the side must NOT BE DUPLICATED
        MapPeriodicSides=.TRUE.
        CYCLE
      END IF
      NBElemID=PartSideToElem(S2E_NB_ELEM_ID,iSide)
      ! only master side is on proc, nothing to do
      IF(NBElemID.LT.1) CYCLE
      IF(NBElemID.GT.nElems) CYCLE
      ! if master and slave side are on proc, duplicate
      nPartPeriodicSides=nPartPeriodicSides+1
      MapPeriodicSides=.TRUE.
    END IF
  END DO
END IF

!IF(nPartPeriodicSides.GT.0)THEN
IF(MapPeriodicSides)THEN
  ! map min-max glob to local array
  MinMaxGlob(1)=GEO%xminglob
  MinMaxGlob(2)=GEO%yminglob
  MinMaxGlob(3)=GEO%zminglob
  MinMaxGlob(4)=GEO%xmaxglob
  MinMaxGlob(5)=GEO%ymaxglob
  MinMaxGlob(6)=GEO%zmaxglob

  ! 2) increase BezierControlPoints and SideXXX from old nSides to nSides+nDuplicatePeriodicSides
  ALLOCATE(DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides))
  ALLOCATE(DummyBezierControlPoints3dElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides))
  ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nTotalSides))
  ALLOCATE(DummySideSlabIntervals(1:6,1:nTotalSides))
  !ALLOCATE(DummySidePeriodicDisplacement(1:3,1:nTotalSides))
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nTotalSides))
  ALLOCATE(DummyBC(1:nTotalSides))
  ALLOCATE(DummyMortarType(1:2,1:nTotalSides))
  ALLOCATE(DummyPartSideToElem(1:5,1:nTotalSides))
  ALLOCATE(DummySidePeriodicType(1:nTotalSides))
  ALLOCATE(DummyMortarSlave2MasterInfo(1:nTotalSides))
  
  ! copy data to backup
  DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides) = BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides)
  DummyBezierControlPoints3dElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides) &
     = BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides)
  DummySideSlabNormals(1:3,1:3,1:nTotalSides)                 = SideSlabNormals(1:3,1:3,1:nTotalSides)
  DummySideSlabIntervals(1:6,1:nTotalSides)                   = SideSlabIntervals(1:6,1:nTotalSides)
  !DummySidePeriodicDisplacement(1:3,1:nTotalSides)            = SidePeriodicDisplacement(1:3,1:nTotalSides)
  DummyBoundingBoxIsEmpty(1:nTotalSides)                      = BoundingBoxIsEmpty(1:nTotalSides)
  DummyBC(1:nTotalSides)                                      = BC(1:nTotalSides)
  DummyMortarSlave2MasterInfo(1:nTotalSides)                  = MortarSlave2MasterInfo(1:nTotalSides)
  DummyMortarType(1:2,1:nTotalSides)                          = MortarType(1:2,1:nTotalSides)
  DummyPartSideTOElem(1:5,1:nTotalSides)                      = PartSideTOElem(1:5,1:nTotalSides)
  DummySidePeriodicType(1:nTotalSides)                        = SidePeriodicType(1:nTotalSides)

  ! deallocate old values and reallocate
  DEALLOCATE(BezierControlPoints3D)
  DEALLOCATE(BezierControlPoints3DElevated)
  DEALLOCATE(SideSlabNormals)
  DEALLOCATE(SideSlabIntervals)
  DEALLOCATE(BoundingBoxIsEmpty)
  DEALLOCATE(MortarSlave2MasterInfo)
  DEALLOCATE(BC)
  DEALLOCATE(MortarType)
  DEALLOCATE(PartSideToElem)
  DEALLOCATE(SidePeriodicType)

  tmpnSides  =nTotalSides 
  nTotalSides=nTotalSides+nPartPeriodicSides
  ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides))
  BezierControlPoints3d=-1.
  ALLOCATE(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nTotalSides))
  BezierControlPoints3DElevated=-1.
  ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalSides))
  SideSlabNormals=-1.
  ALLOCATE(SideSlabIntervals(1:6,1:nTotalSides))
  SideSlabIntervals=-1.
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalSides))
  BoundingBoxIsEmpty=.FALSE.
  ALLOCATE(BC(1:nTotalSides))
  BC=-3
  ALLOCATE(MortarSlave2MasterInfo(1:nTotalSides))
  MortarSlave2MasterInfo=-1
  ALLOCATE(MortarType(1:2,1:nTotalSides))
  MortarType=-1
  ALLOCATE(PartSideToElem(1:5,1:nTotalSides))
  PartSideToElem=-1
  ALLOCATE(SidePeriodicType(1:nTotalSides))
  SidePeriodicType=0
  !ALLOCATE(SidePeriodicDisplacement(1:3,1:nTotalSides))

  BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:tmpnSides) = DummyBezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:tmpnSides)
  BezierControlPoints3dElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:tmpnSides) &
     = DummyBezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:tmpnSides)
  SideSlabNormals(1:3,1:3,1:tmpnSides)                 = DummySideSlabNormals(1:3,1:3,1:tmpnSides)
  SideSlabIntervals(1:6,1:tmpnSides)                   = DummySideSlabIntervals(1:6,1:tmpnSides)
  BoundingBoxIsEmpty(1:tmpnSides)                      = DummyBoundingBoxIsEmpty(1:tmpnSides)
  !SidePeriodicDisplacement(1:3,1:nTotalSides)          = DummySidePeriodicDisplacement(1:3,1:nTotalSides)
  BC(1:tmpnSides)                                      = DummyBC(1:tmpnSides)
  MortarSlave2MasterInfo(1:tmpnSides)                  = DummyMortarSlave2MasterInfo(1:tmpnSides)
  MortarType(1:2,1:tmpnSides)                          = DummyMortarType(1:2,1:tmpnSides)
  PartSideToElem(1:5,1:tmpnSides)                      = DummyPartSideTOElem(1:5,1:tmpnSides)
  SidePeriodicType(1:tmpnSides)                        = DummySidePeriodicType(1:tmpnSides)

  ! 3) loop over the OLD sides and copy the corresponding SideXXX. The missing BezierControlPoints (periodic shifted values) 
  !    are build from the other element. Now two BezierControlPoints exists which are shifted by the SidePeriodicVector
  nPartPeriodicSides=0
  DO iSide=1,tmpnSides
    IF(SidePeriodicType(iSide).NE.0)THEN
      NBElemID=PartSideToElem(S2E_NB_ELEM_ID,iSide)
      IF(NBElemID.LT.1) CYCLE
      IF(NBElemID.GT.nElems) CYCLE
      NBlocSideID=PartSideToElem(S2E_NB_LOC_SIDE_ID,iSide)
      flip=PartSideToElem(S2E_FLIP,iSide)
      locSideID=PartSideToElem(S2E_LOC_SIDE_ID,iSide)
      ElemID   =PartSideToElem(S2E_ELEM_ID,iSide)
      ! 3b) set newSideID and sidedata
      IF(ElemID.EQ.-1) THEN
        ! MPI side
        newSideID=iSide
        PVID=SidePeriodicType(iSide)
        SidePeriodicType(newSideID)=-SidePeriodicType(iSide) ! stored the inital alpha value
      ELSE
        nPartPeriodicSides=nPartPeriodicSides+1
        newSideID=tmpnSides+nPartPeriodicSides
        ! bc
        BCID = BC(iSide)
        PVID = BoundaryType(BCID,BC_ALPHA) 
        ! loop over bc to get the NEW BC type
        DO iBC = 1,nBCs
          IF(BoundaryType(iBC,BC_ALPHA).EQ.-PVID) THEN
            BC(newSideID)=iBC 
            EXIT
          END IF
        END DO
        MortarSlave2MasterInfo(newSideID) = DummyMortarSlave2MasterInfo(iSide)
        PVID=SidePeriodicType(iSide)
        SidePeriodicType(newSideID)=-SidePeriodicType(iSide) ! stored the inital alpha value
        ! rebuild sides for sanity
        CALL GetBezierControlPoints3D(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID),ElemID,ilocSide_In=locSideID,SideID_In=iSide)
        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)                             &
                                           ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide)     &
                                           ,SideSlabNormals(1:3,1:3,iSide)                                             &
                                           ,SideSlabInterVals(1:6,iSide)                                               &
                                           ,BoundingBoxIsEmpty(iSide)                                                  )
        ! sanity check
        ! check flip of master element, has to be zero, because of master side
        IF(PartElemToSide(E2S_FLIP   ,locSideID,ElemID).NE.0)THEN
          CALL abort(&
                __STAMP__&
                , ' Periodic Side is no master side!')
        END IF
        xTest(1:3) = BezierControlPoints3D(1:3,0,0,iSide)
        xTest      = xTest + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
        IF(xTest(1)+1e-8.LT.MinMaxGlob(1)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
        IF(xTest(2)+1e-8.LT.MinMaxGlob(2)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
        IF(xTest(3)+1e-8.LT.MinMaxGlob(3)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
        IF(xTest(1)-1e-8.GT.MinMaxGlob(4)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
        IF(xTest(2)-1e-8.GT.MinMaxGlob(5)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
        IF(xTest(3)-1e-8.GT.MinMaxGlob(6)) SidePeriodicType(iSide)=-SidePeriodicType(iSide)
      END IF
      ! the flip has to be set to -1, artificial master side
      PartElemToSide(E2S_FLIP   ,NBlocSideID,NBElemID) = 0
      PartElemToSide(E2S_SIDE_ID,NBlocSideID,NBElemID) = newSideID
      ! rebuild BezierControlPoints3D (simplified version, all BezierControlPoints3D are rebuild)
      CALL GetBezierControlPoints3D(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,NBElemID),NBElemID,ilocSide_In=NBlocSideID,SideID_In=NewSideID)
      ! remains equal because of MOVEMENT and MIRRORING of periodic side
      ! periodic displacement 
      !DO q=0,NGeo
      !  DO p=0,NGeo
      !    BezierControlPoints3d(1:3,p,q,newSideID)  = DummyBezierControlPoints3d(1:3,p,q,iSide) &
      !                                              + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
      !  END DO ! p=0,NGeo
      !END DO ! q=0,NGeo
      !! recompute quark
      !CALL RotateMasterToSlave(flip,NBlocSideID,BezierControlPoints3d(1:3,0:NGeo,0:NGeo,newSideID))
      DO idir=1,3
        MinMax(1)=MINVAL(BezierControlPoints3d(iDir,:,:,newSideID))
        MinMax(2)=MAXVAL(BezierControlPoints3d(iDir,:,:,newSideID))
        ! this may be required a tolerance due to periodic displacement
        IF(.NOT.ALMOSTEQUALRELATIVE(MinMax(1),MinMaxGlob(iDir),1e-10))THEN
          IF(MinMax(1).LT.MinMaxGlob(iDir)) THEN
            IPWRITE(UNIT_stdOut,*) ' Min-comparison. MinValue, GlobalMin ', MinMax(1),MinMaxGlob(iDir)
            CALL abort(&
                __STAMP__&
                , ' BezierControlPoints3d is moved outside of minvalue of GEO%glob! Direction', iDir)
          END IF
        ELSE
          IF(printBezierControlPointsWarnings)THEN
            IPWRITE(UNIT_stdOut,*) ' WARNING: Min-comparison. MinValue, GlobalMin ', MinMax(1),MinMaxGlob(iDir)
          END IF
        END IF
        IF(.NOT.ALMOSTEQUALRELATIVE(MinMax(2),MinMaxGlob(iDir+3),1e-10))THEN
          IF(MinMax(2).GT.MinMaxGlob(iDir+3)) THEN
            IPWRITE(UNIT_stdOut,*) ' Max-comparison MaxValue, GlobalMax ', MinMax(2),MinMaxGlob(iDir+3)
            CALL abort(&
                __STAMP__&
                , ' BezierControlPoints3d is moved outside of maxvalue of GEO%glob! Direction', iDir)
          END IF
        ELSE
          IF(printBezierControlPointsWarnings)THEN
            IPWRITE(UNIT_stdOut,*) ' WARNING: Max-comparison MaxValue, GlobalMax ', MinMax(2),MinMaxGlob(iDir+3)
          END IF
        END IF
      END DO

      ! fill partsidetoelem
      PartSideToElem(S2E_ELEM_ID,newSideID)=NBElemID
      PartSideToElem(S2E_NB_ELEM_ID,newSideID)=ElemID
      PartSideToElem(S2E_FLIP,newSideID)=-1
      PartSideToElem(S2E_LOC_SIDE_ID,newSideID)=NBlocSideID
      PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID)=locSideID
      ! mortar type
      MortarType(1:2,newSideID) = DummyMortarType(1:2,iSide)
      ! bounding box, etc...
      CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)                         &
                                         ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,newSideID) &
                                         ,SideSlabNormals(1:3,1:3,newSideID)                                         &
                                         ,SideSlabInterVals(1:6,newSideID)                                           &
                                         ,BoundingBoxIsEmpty(newSideID)                                              )

      ! sanity check
      xTest(1:3) = BezierControlPoints3D(1:3,0,0,newSideID)
      PVID=SidePeriodicType(newSideID)
      xTest      = xTest + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
      IF(xTest(1)+1e-8.LT.MinMaxGlob(1)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
      IF(xTest(2)+1e-8.LT.MinMaxGlob(2)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
      IF(xTest(3)+1e-8.LT.MinMaxGlob(3)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
      IF(xTest(1)-1e-8.GT.MinMaxGlob(4)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
      IF(xTest(2)-1e-8.GT.MinMaxGlob(5)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
      IF(xTest(3)-1e-8.GT.MinMaxGlob(6)) SidePeriodicType(newSideID)=-SidePeriodicType(newSideID)
    END IF
  END DO ! iSide=1,tmpnSides
  ! deallocate dummy  
  DEALLOCATE(DummyBezierControlPoints3D)
  DEALLOCATE(DummySideSlabNormals)
  DEALLOCATE(DummySideSlabIntervals)
  DEALLOCATE(DummyBoundingBoxIsEmpty)
  DEALLOCATE(DummyBC)
  DEALLOCATE(DummyMortarType)
  DEALLOCATE(DummyPartSideToElem)
  DEALLOCATE(DummySidePeriodicType)

END IF ! nPartPeriodicSides .GT.0
! reset side-counter
nPartSides     =nPartPeriodicSides+nSides
nTotalBCSides =nPartPeriodicSides+nSides

! sanity check for PartElemToSide
DO iElem=1,nElems
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    IF(MortarType(1,SideID).EQ.0)THEN
      IF(SideID.LE.0)THEN
        CALL abort(&
__STAMP__&
      , ' No Side ID set. critical error!',iElem,REAL(ilocSide))
      END IF
    END IF
  END DO
END DO ! iElem=1,PP_nElems

#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iERROR)
#endif
SWRITE(UNIT_StdOut,'(A)') ' Sanity check of duplication successful!'

END SUBROUTINE DuplicateSlavePeriodicSides


SUBROUTINE MarkAllBCSides() 
!===================================================================================================================================
! CAUTION: nTotalBCSides is reset from old value to new current,process-local BCSides
! The PartBCSideList contains a mapping from the global side list to a local, pure BC side list
! DG-SideList
! 1:nBCSides - nInnerSides - nMortarSides - nMPISides
! DG: periodic sides are no BC sides
! ParticleTracking treats periodic sides as BC sides and the process needs the actual side,
! hence it may be required to be duplicated (side at correct position)
! Particle-Tracking-List before MarkAllBCSides
! 1:nBCSides - nInnerSides - nSomePeriodicSides - nMortarSides - nMPISides - nMissingPeriodicSides
! As RefMapping requires only the BC sides, a shorter list is generated over all
! nTotalBCSides which is NOW smaller than nPartSides or nTotalSides
! CAUTION and BRAIN-FUCK: 
! This smaller list is used to build: from 1:nTotalBCSides < nTotalSides and is used for
! SideNormVec,SideTypes,SideDistance
! BUT: 1:nTotalSides is STILL used for 
! BezierControlPoints3D, SideSlabInterVals,SideSlabNormals,BoundingBoxIsEmpty
! and are NOT reshaped yet, hence, the length of the array remains nTotalSides
! BRAIN-FUCK CONTINUOUS: 
! During building of the HALO region, the BezierControlPoints variables are further increased with nTotalSides while the 
! already small arrays increases with nTotalBCSides
! After building the HALO region, the actual arrays are reshaped and a stored in shorter arrays
! 
! AND no rule without a break:
! SidePeriodicType is still on nTotalSides and NOT reshaped
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Mesh_Vars,               ONLY:nSides
USE MOD_Particle_Mesh_Vars,      ONLY:PartBCSideList,nTotalSides,nPartPeriodicSides,nTotalBCSides,nPartSides
USE MOD_Mesh_Vars,               ONLY:BC,nBCSides
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
USE MOD_Globals
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iSide
!===================================================================================================================================

! PartBCSideList is increased, due to the periodic sides
IF(.NOT.DoRefMapping) RETURN

DEALLOCATE(PartBCSideList)
ALLOCATE(PartBCSideList(1:nTotalSides))
! BC Sides 
PartBCSideList=-1
nTotalBCSides=0
DO iSide=1,nBCSides
  nTotalBCSides=nTotalBCSides+1
  PartBCSideList(iSide)=nTotalBCSides
END DO

DO iSide=nBCSides+1,nSides+nPartPeriodicSides
  IF(BC(iSide).EQ.0) CYCLE
  nTotalBCSides=nTotalBCSides+1
  PartBCSideList(iSide)=nTotalBCSides
END DO ! iSide

! nPartsides 
nPartSides   =nPartPeriodicSides+nSides

END SUBROUTINE MarkAllBCSides


SUBROUTINE BGMIndexOfElement(ElemID,ElemToBGM) 
!===================================================================================================================================
! computes the element indices of an given element in the BGM-mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ChangeBasis,                        ONLY:ChangeBasis2D
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Surfaces_Vars,             ONLY:sVdm_Bezier
USE MOD_Mesh_Vars,                          ONLY:XCL_NGeo
USE MOD_Mesh_Vars,                          ONLY:NGeo
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO
USE MOD_Particle_Tracking_Vars,             ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,                 ONLY:PartElemToSide
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)        :: ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: ElemToBGM(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ilocSide, SideID
REAL                      :: xmin,xmax,ymin,ymax,zmin,zmax
REAL                      :: BezierControlPoints3D_tmp(1:3,0:NGeo,0:NGeo)
!===================================================================================================================================

xmin = HUGE(1.0)
xmax =-HUGE(1.0)
ymin = HUGE(1.0)
ymax =-HUGE(1.0)
zmin = HUGE(1.0)
zmax =-HUGE(1.0)

! get min,max of BezierControlPoints of Element
DO iLocSide = 1,6
  SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, ElemID)
  IF(DoRefMapping)THEN
    IF(SideID.GT.0)THEN
      IF(PartElemToSide(E2S_FLIP,ilocSide,ElemID).EQ.0)THEN
        BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
      ELSE
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
        CASE(XI_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ETA_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ETA_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ZETA_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
        CASE(ZETA_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
        END SELECT
      END IF
    ELSE
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
      CASE(XI_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ZETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
      CASE(ZETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
      END SELECT
    END IF
  ELSE ! pure tracing
    BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
  END IF
  xmin=MIN(xmin,MINVAL(BezierControlPoints3D_tmp(1,:,:)))
  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D_tmp(1,:,:)))
  ymin=MIN(ymin,MINVAL(BezierControlPoints3D_tmp(2,:,:)))
  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D_tmp(2,:,:)))
  zmin=MIN(zmin,MINVAL(BezierControlPoints3D_tmp(3,:,:)))
  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D_tmp(3,:,:)))
END DO ! ilocSide

ElemToBGM(1) = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
ElemToBGM(2) = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
ElemToBGM(3) = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
ElemToBGM(4) = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
ElemToBGM(5) = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
ElemToBGM(6) = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))


END SUBROUTINE BGMIndexOfElement


SUBROUTINE GetFIBGMMinMax() 
!===================================================================================================================================
! computes the minimum and maximum value of the FIBGM mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO
USE MOD_Mesh_Vars,                          ONLY:MortarSlave2MasterInfo
USE MOD_Particle_Mesh_Vars,                 ONLY:GEO,nTotalSides
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D
#ifdef MPI
USE MOD_Particle_MPI_Vars,                  ONLY:PartMPI
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iSide
REAL            :: xmin, xmax, ymin, ymax, zmin, zmax
!===================================================================================================================================

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
!DO iElem=1,nTotalElems
!  xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
!  xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
!  ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
!  ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
!  zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
!  zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
!END DO ! iElem

! bounding box!!
DO iSide=1,nTotalSides
  IF(MortarSlave2MasterInfo(iSide).NE.-1) CYCLE
  xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,iSide)))
  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,iSide)))
  ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,iSide)))
  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,iSide)))
  zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,iSide)))
  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,iSide)))
END DO ! iSide

GEO%xmin=xmin
GEO%xmax=xmax
GEO%ymin=ymin
GEO%ymax=ymax
GEO%zmin=zmin
GEO%zmax=zmax

#ifdef MPI
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

END SUBROUTINE GetFIBGMMinMax


SUBROUTINE GetSideOriginAndRadius(nTotalBCSides,SideOrigin,SideRadius)
!===================================================================================================================================
! ONLY RefMapping
! Computes the side origin and radius for each BC Side
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Mesh_Vars,              ONLY:NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:PartBCSideList,nTotalSides
USE MOD_Basis,                  ONLY:DeCasteljauInterpolation
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3d
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)       :: nTotalBCSides
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)         :: SideOrigin(1:3,1:nTotalBCSides),SideRadius(1:nTotalBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iSide, BCSideID,p,q
REAL                     :: Xi(1:2), Origin(1:3), Radius, RadiusMax, Vec(1:3)
!===================================================================================================================================

SideOrigin=0.
SideRadius=0.

DO iSide=1,nTotalSides
  BCSideID=PartBCSideList(iSide)
  IF(BCSideID.LT.1) CYCLE
  Xi=0.
  CALL DeCasteljauInterpolation(NGeo,Xi(1:2),BCSideID,Origin)
  SideOrigin(1:3,BCSideID) = Origin
  Radius=0.
  RadiusMax=0.
  DO q=0,NGeo
    DO p=0,NGeo
      Vec(1:3) = BezierControlPoints3d(:,p,q,BCSideID)-Origin
      Radius=DOT_PRODUCT(Vec,Vec)
      RadiusMax=MAX(RadiusMax,Radius)
    END DO ! p=0,NGeo
  END DO ! q=0,NGeo
  SideRadius(BCSideID)=SQRT(RadiusMax)
END DO ! iSide=1,nTotalSides

END SUBROUTINE GetSideOriginAndRadius


SUBROUTINE GetElemToSideDistance(nTotalBCSides,SideOrigin,SideRadius)
!===================================================================================================================================
! computes the distance between each element and it associated sides for DoRefMapping=T
! only sides for which ElemToSideDistance<lengthPartTrajectory have to be checked during the current tracing step
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
USE MOD_Particle_Mesh_Vars,     ONLY:IsTracingBCElem,ElemRadiusNGeo,BCElem,PartBCSideList,nTotalElems
USE MOD_Utils,                  ONLY:InsertionSort
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)       :: nTotalBCSides
REAL,INTENT(IN)          :: SideOrigin(1:3,1:nTotalBCSides),SideRadius(1:nTotalBCSides)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,ilocSide,SideID,BCSideID
REAL                     :: Vec(1:3)
REAL                     :: Origin(1:3)
!===================================================================================================================================

! loop over all  elements
DO iElem=1,nTotalElems
  IF(.NOT.IsTracingBCElem(iElem)) CYCLE
  ALLOCATE( BCElem(iElem)%ElemToSideDistance(BCElem(iElem)%lastSide) )
  BCElem(iElem)%ElemToSideDistance(BCElem(iElem)%lastSide)=0.
  Origin(1:3) = ElemBaryNGeo(1:3,iElem)
  ! loop over all associated sides
  DO iLocSide=1,BCElem(iElem)%lastSide
    SideID=BCElem(iElem)%BCSideID(ilocSide)
    BCSideID=PartBCSideList(SideID)
    Vec=Origin - SideOrigin(1:3,BCSideID)
    BCElem(iElem)%ElemToSideDistance(ilocSide) = SQRT(DOT_PRODUCT(Vec,Vec))-ElemRadiusNGeo(iElem)-SideRadius(BCSideID)
  END DO ! iLocSide=1,BCElem(iElem)%lastSide
  ! sort each side distance for each element according to it's distance
  CALL InsertionSort(BCElem(iElem)%ElemToSideDistance(:),BCElem(iElem)%BCSideID(:),BCElem(iElem)%lastSide)
END DO ! iElem=1,PP_nElems

END SUBROUTINE GetElemToSideDistance


SUBROUTINE MarkAuxBCElems()
!===================================================================================================================================
! check if auxBCs are inside BoundingBox of Elems
! -- plane: use plane equation f=a1*x+a2*y+a3*z+a4=0 and insert corresponding intervals of box -> fmin and fmax
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,                 ONLY:ElemHasAuxBCs
USE MOD_Particle_Boundary_Vars,             ONLY:nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,iAuxBC,icoord,dir(3),positiontype,positiontype_tmp
REAL                     :: r_vec(3),n_vec(3),fmin,fmax,Bounds(1:2,1:3),radius,BoundsBC(1:2,1:3)
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
      CALL BoundsOfElement(iElem,Bounds)
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
        CALL BoundsOfElement(iElem,Bounds)
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

SUBROUTINE BoundsOfElement(ElemID,Bounds)
!===================================================================================================================================
! computes the min/max of element in xyz (Based on BGMIndexOfElement)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ChangeBasis,                        ONLY:ChangeBasis2D
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D,sVdm_Bezier
USE MOD_Mesh_Vars,                          ONLY:XCL_NGeo
USE MOD_Mesh_Vars,                          ONLY:NGeo
USE MOD_Particle_Tracking_Vars,             ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,                 ONLY:PartElemToSide
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: Bounds(1:2,1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ilocSide, SideID
REAL                      :: xmin,xmax,ymin,ymax,zmin,zmax
REAL                      :: BezierControlPoints3D_tmp(1:3,0:NGeo,0:NGeo)
!===================================================================================================================================

xmin = HUGE(1.0)
xmax =-HUGE(1.0)
ymin = HUGE(1.0)
ymax =-HUGE(1.0)
zmin = HUGE(1.0)
zmax =-HUGE(1.0)

! get min,max of BezierControlPoints of Element
DO iLocSide = 1,6
  SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, ElemID)
  IF(DoRefMapping)THEN
    IF(SideID.GT.0)THEN
      IF(PartElemToSide(E2S_FLIP,ilocSide,ElemID).EQ.0)THEN
        BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
      ELSE
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
        CASE(XI_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ETA_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ETA_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ZETA_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
        CASE(ZETA_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
        END SELECT
      END IF
    ELSE
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
      CASE(XI_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ZETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
      CASE(ZETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
      END SELECT
    END IF
  ELSE ! pure tracing
    BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
  END IF
  xmin=MIN(xmin,MINVAL(BezierControlPoints3D_tmp(1,:,:)))
  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D_tmp(1,:,:)))
  ymin=MIN(ymin,MINVAL(BezierControlPoints3D_tmp(2,:,:)))
  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D_tmp(2,:,:)))
  zmin=MIN(zmin,MINVAL(BezierControlPoints3D_tmp(3,:,:)))
  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D_tmp(3,:,:)))
END DO ! ilocSide
Bounds(:,1)=(/xmin,xmax/)
Bounds(:,2)=(/ymin,ymax/)
Bounds(:,3)=(/zmin,zmax/)

END SUBROUTINE BoundsOfElement


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

END MODULE MOD_Particle_Mesh
