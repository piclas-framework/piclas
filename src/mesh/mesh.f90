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

MODULE MOD_Mesh
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

INTERFACE GetMeshMinMaxBoundaries
  MODULE PROCEDURE GetMeshMinMaxBoundaries
END INTERFACE

PUBLIC::InitMesh
PUBLIC::FinalizeMesh
PUBLIC::GetMeshMinMaxBoundaries
PUBLIC::DefineParametersMesh

INTEGER,PARAMETER :: PRM_P_ADAPTION_ZERO = 0  ! deactivate
INTEGER,PARAMETER :: PRM_P_ADAPTION_RDN  = 1  ! random
INTEGER,PARAMETER :: PRM_P_ADAPTION_NPB  = 2  ! Elements with non-periodic boundary sides get NMax
INTEGER,PARAMETER :: PRM_P_ADAPTION_HH   = 3  ! Elements in the lower half domain in x-direction are set to NMin and the upper half are set to NMax. Origin must be at x=0.

INTEGER,PARAMETER :: PRM_P_ADAPTION_LVL_MINTWO  = -2 ! Directly connected elements are set to NMax and 2nd layer elements are set to NMin+1
INTEGER,PARAMETER :: PRM_P_ADAPTION_LVL_MINONE  = -1 ! Directly connected elements are set to NMin+1
INTEGER,PARAMETER :: PRM_P_ADAPTION_LVL_DEFAULT = 1 ! Directly connected elements are set to NMax
INTEGER,PARAMETER :: PRM_P_ADAPTION_LVL_TWO     = 2 ! Directly connected elements and elements connected with these are set to NMax
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for Mesh
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile'            , "(relative) path to meshfile (mandatory)\n(HALOWIKI:) usually located in directory of project.ini")
CALL prms%CreateLogicalOption( 'useCurveds'          , "Controls usage of high-order information in mesh. Turn off to discard high-order data and treat curved meshes as linear meshes."    , '.FALSE.')
CALL prms%CreateRealOption(    'meshScale'           , "Scale the mesh by this factor (shrink/enlarge)."                                                                                    , '1.0')
CALL prms%CreateLogicalOption( 'meshdeform'          , "Apply simple sine-shaped deformation on cartesion mesh (for testing)."                                                              , '.FALSE.')
CALL prms%CreateLogicalOption( 'meshCheckRef'        , "Flag if the mesh Jacobians should be checked in the reference system in addition to the computational system."                      , '.TRUE.')
CALL prms%CreateLogicalOption( 'CalcMeshInfo'        , 'Calculate and output elem data for myrank, ElemID and tracking info to ElemData'                                                    , '.FALSE.')
CALL prms%CreateLogicalOption( 'readFEMconnectivity' , 'Activate reading the FEM connectivity arrays EdgeInfo, EdgeConnectInfo, VertexInfo and VertexConnectInfo from the mesh file'        , '.FALSE.')
CALL prms%CreateLogicalOption( 'crossProductMetrics' , "Compute mesh metrics using cross product form. Caution: in this case free-stream preservation is only guaranteed for N=3*NGeo."     , '.FALSE.')
CALL prms%CreateStringOption(  'BoundaryName'        , "Names of boundary conditions to be set (must be present in the mesh!). For each BoundaryName a BoundaryType needs to be specified." , multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType'        , "Type of boundary conditions to be set. Format: (BC_TYPE, BC_STATE)"                                                                 , multiple=.TRUE. , no=2)
CALL prms%CreateLogicalOption( 'writePartitionInfo'  , "Write information about MPI partitions into a file."                                                                                , '.FALSE.')

! p-adaption
CALL prms%CreateIntFromStringOption('pAdaptionType', "Method for initial polynomial degree distribution among the elements: \n"//&
                                    '           none ('//TRIM(int2strf(PRM_P_ADAPTION_ZERO))//'): default for setting all elements to N\n'//&
                                    '         random ('//TRIM(int2strf(PRM_P_ADAPTION_RDN))//'): elements get random polynomial degree between NMin and NMax\n'//&
                                    'non-periodic-BC ('//TRIM(int2strf(PRM_P_ADAPTION_NPB))//'): elements with non-periodic boundary conditions receive NMax\n'//&
                                    '      half-half ('//TRIM(int2strf(PRM_P_ADAPTION_HH))//'): elements in the lower half domain in x-direction are set to NMin and the upper half are set to NMax. The domain must be centered around x=0.\n'&
                                   ,'none')

CALL addStrListEntry('pAdaptionType' , 'none'            , PRM_P_ADAPTION_ZERO)
CALL addStrListEntry('pAdaptionType' , 'random'          , PRM_P_ADAPTION_RDN)
CALL addStrListEntry('pAdaptionType' , 'non-periodic-BC' , PRM_P_ADAPTION_NPB)
CALL addStrListEntry('pAdaptionType' , 'half-half'       , PRM_P_ADAPTION_HH)

CALL prms%CreateIntFromStringOption('pAdaptionBCLevel', "Only for pAdaptionType=non-periodic-BC: Number/Depth of elements connected to a boundary that are set to NMax.\n"//&
                                    '1st-and-2nd-NMin+1 ('//TRIM(int2strf(PRM_P_ADAPTION_LVL_MINTWO))//'): elements with non-periodic boundary conditions receive NMax, 2nd layer receive NMin+1\n'//&
                                    '    directly-connected-NMin+1 ('//TRIM(int2strf(PRM_P_ADAPTION_LVL_MINONE))//'): elements with non-periodic boundary conditions receive NMin+1\n'//&
                                    '      directly-connected-NMax ('//TRIM(int2strf(PRM_P_ADAPTION_LVL_DEFAULT))//'): elements with non-periodic boundary conditions receive NMax\n'//&
                                    '  1st-and-2nd-NMax ('//TRIM(int2strf(PRM_P_ADAPTION_LVL_TWO))//'): first two elements with non-periodic boundary conditions receive NMax\n'&
                                   ,'directly-connected')

CALL addStrListEntry('pAdaptionBCLevel' , '1st-and-2nd-NMin+1'        , PRM_P_ADAPTION_LVL_MINTWO)
CALL addStrListEntry('pAdaptionBCLevel' , 'directly-connected-NMin+1' , PRM_P_ADAPTION_LVL_MINONE)
CALL addStrListEntry('pAdaptionBCLevel' , 'directly-connected-NMax'   , PRM_P_ADAPTION_LVL_DEFAULT)
CALL addStrListEntry('pAdaptionBCLevel' , '1st-and-2nd-NMax'          , PRM_P_ADAPTION_LVL_TWO)

END SUBROUTINE DefineParametersMesh

!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!==================================================================================================================================
SUBROUTINE InitMesh(meshMode,MeshFile_IN)
!===================================================================================================================================
! Read Parameter from inputfile
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_PreProc
USE MOD_Analyze_Vars           ,ONLY: CalcMeshInfo
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_HDF5_Input
USE MOD_Interpolation_Vars     ,ONLY: InterpolationInitIsDone,Nmin,Nmax
USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
USE MOD_Mappings               ,ONLY: InitMappings
USE MOD_Mesh_Vars
USE MOD_Mesh_ReadIn            ,ONLY: ReadMesh
USE MOD_Metrics                ,ONLY: BuildElem_xGP,CalcMetrics,CalcSurfMetrics
USE MOD_Prepare_Mesh           ,ONLY: setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools            ,ONLY: PrintOption
USE MOD_ReadInTools            ,ONLY: GETLOGICAL,GETSTR,GETREAL,GETINT,GETREALARRAY
#if USE_HDG
USE MOD_Symmetry_Vars          ,ONLY: Symmetry
#endif
#if USE_MPI
USE MOD_Prepare_Mesh           ,ONLY: exchangeFlip
!USE MOD_DG_Vars                ,ONLY: N_DG_Mapping_Shared_Win
!USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_LEADERS_SHARED, MPI_COMM_SHARED, myComputeNodeRank, myleadergrouprank, nComputeNodeProcessors
!USE MOD_MPI_Shared_Vars        ,ONLY: nLeaderGroupProcs
!USE MOD_Particle_Mesh_Vars     ,ONLY: offsetComputeNodeElem
USE MOD_MPI_Shared
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Metrics    ,ONLY: ExchangeVolMesh,ExchangeMetrics
USE MOD_LoadBalance_Vars       ,ONLY: DoLoadBalance,PerformLoadBalance,UseH5IOLoadBalance
USE MOD_Output_Vars            ,ONLY: DoWriteStateToHDF5
USE MOD_Restart_Vars           ,ONLY: DoInitialAutoRestart
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_Particle_Vars          ,ONLY: usevMPF
#endif
#if USE_HDG && USE_LOADBALANCE
USE MOD_Mesh_Tools             ,ONLY: BuildSideToNonUniqueGlobalSide
#endif /*USE_HDG && USE_LOADBALANCE*/
USE MOD_DG_Vars                ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave
USE MOD_Particle_Mesh_Vars     ,ONLY: meshScale
USE MOD_Mesh_Vars              ,ONLY: firstMortarInnerSide,lastMortarInnerSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !<  0: only read and build Elem_xGP,
                               !< -2: as 0 + build connectivity and always set ReadNodes=T: read node info (automatically read for PARTICLES=ON) + build FacexGP and keep NodeCoords
                               !< -1: as 0 + build connectivity and always set ReadNodes=T: read node info (automatically read for PARTICLES=ON)
                               !<  1: as 0 + build connectivity
                               !<  2: as 1 + calc metrics
                               !<  3: as 2 but skip InitParticleMesh
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: MeshFile_IN !< file name of mesh to be read
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: x(3)
REAL,POINTER        :: Coords(:,:,:,:,:)
INTEGER             :: iElem,i,j,k,nElemsLoc
INTEGER             :: Nloc,iSide,NSideMin
LOGICAL             :: validMesh,ReadNodes
#if USE_HDG
INTEGER             :: iMortar,nMortars,MortarSideID
INTEGER             :: SideID,locSide
#endif /*USE_HDG*/
!===================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL abort(__STAMP__,'InitMesh not ready to be called or already called.')
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'
#if defined(PARTICLES)
ReadNodes  =.TRUE.
#else
ReadNodes  =.FALSE.
IF(meshMode.LT.0) ReadNodes  =.TRUE.
#endif /*defined(PARTICLES)*/

! Output of myrank, ElemID and tracking info
CalcMeshInfo = GETLOGICAL('CalcMeshInfo')

! prepare pointer structure (get nElems, etc.)
IF (PRESENT(MeshFile_IN)) THEN
  MeshFile = MeshFile_IN
ELSE
  MeshFile = GETSTR('MeshFile')
END IF

#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  validMesh = ISVALIDMESHFILE(MeshFile)
  IF(.NOT.validMesh) CALL CollectiveStop(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/


useCurveds=GETLOGICAL('useCurveds')
#if USE_LOADBALANCE
IF ( (DoLoadBalance.OR.DoInitialAutoRestart) .AND. (.NOT.DoWriteStateToHDF5) .AND. UseH5IOLoadBalance) THEN
  DoWriteStateToHDF5=.TRUE.
  CALL PrintOption('Loadbalancing or InitialAutoRestart enabled: DoWriteStateToHDF5','INFO',LogOpt=DoWriteStateToHDF5)
END IF
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)
  CALL PrintOption('NGeo','INFO',IntOpt=NGeo)
  CALL CloseDataFile()
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

IF(useCurveds)THEN
  IF(PP_N.LT.NGeo)THEN
    LBWRITE(UNIT_stdOut,'(A)') 'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                             & can cause problems on periodic boundaries! Set N>NGeo'
  END IF
ELSE
  IF(NGeo.GT.1) THEN
    LBWRITE(UNIT_stdOut,'(A)') ' WARNING: Using linear elements although NGeo>1! NGeo will be set to 1 after coordinates read-in!'
  END IF
END IF

meshScale = GETREAL('meshScale') ! default is 1.0
! Sanity check
IF(ABS(meshScale).LE.0.) CALL abort(__STAMP__,'meshScale is zero')

! Activate reading the FEM connectivity arrays EdgeInfo, EdgeConnectInfo, VertexInfo and VertexConnectInfo from the mesh file
readFEMconnectivity = GETLOGICAL('readFEMconnectivity') ! This flag is required in ReadMesh() to load the FEM info from the .h5 mesh file
#if USE_MPI
IF(readFEMconnectivity)THEN
  ELEM_RANK     = 11
  ELEM_HALOFLAG = 12
ELSE
  ELEM_RANK     = 7
  ELEM_HALOFLAG = 8
END IF ! readFEMconnectivity
#endif  /*USE_MPI*/

CALL ReadMesh(MeshFile,ReadNodes) !set nElems

!schmutz fink
PP_nElems=nElems

Coords=>NodeCoords
nElemsLoc=nElems

! scale and deform mesh only if not already done in ReadMeshNodes()
IF(.NOT.ReadNodes)THEN
  IF(ABS(meshScale-1.).GT.1e-14) Coords = Coords*meshScale
END IF ! .NOT.ReadNodes

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  DO iElem=1,nElems
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PI*x(1))*SIN(PI*x(2))*SIN(PI*x(3))
    END DO; END DO; END DO;
  END DO
END IF

! Return if no connectivity and metrics are required (e.g. for visualization mode)
IF (ABS(meshMode).GT.0) THEN
  LBWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
  CALL setLocalSideIDs()

#if USE_MPI
  ! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
  LBWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
  CALL exchangeFlip()
#endif

  ! Set the side ranges here (because by now nMortarInnerSides, nMPISides_MINE and nMPISides_YOUR have been determined)
  ! and calculate nGlobalUniqueSides (also required are nBCSides and nInnerSides, which have been determined in ReadMesh())
  ! Requires:
  !   nBCSides          is set in ReadMesh()
  !   nMortarInnerSides is set in setLocalSideIDs()
  !   nInnerSides       is set in ReadMesh()
  !   nMPISides_MINE    is set in setLocalSideIDs()
  !   nMPISides_YOUR    is set in setLocalSideIDs()
  CALL setSideRanges()

  ! fill ElemToSide, SideToElem,BC
  ALLOCATE(ElemToSide(2,6,nElems))
  ALLOCATE(SideToElem(5,nSides))
  ALLOCATE(BC(1:nSides))
  ALLOCATE(AnalyzeSide(1:nSides))
  ElemToSide  = 0
  SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
  BC          = 0
  AnalyzeSide = 0

  ! fill output definition for InnerBCs
#if defined(PARTICLES) || USE_HDG
  ALLOCATE(GlobalUniqueSideID(1:nSides))
  GlobalUniqueSideID(:)=-1
#endif /*defined(PARTICLES) || USE_HDG*/

  !NOTE: nMortarSides=nMortarInnerSides+nMortarMPISides
  ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
  ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
  ALLOCATE(MortarSlave2MasterInfo(1:nSides))
  MortarType=-1
  MortarInfo=-1

  LBWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
#if USE_HDG && USE_LOADBALANCE
  ! Call with meshMode to check whether, e.g., HDG load balance info need to be determined or not
  CALL fillMeshInfo(meshMode) ! Fills SideToElem and ElemToSide
#else
  CALL fillMeshInfo() ! Fills SideToElem and ElemToSide
#endif /*USE_HDG && USE_LOADBALANCE*/

  ! build necessary mappings
  ALLOCATE(N_Mesh(Nmin:Nmax))
  DO Nloc = Nmin, Nmax
    CALL InitMappings(Nloc, N_Mesh(Nloc)%VolToSideA, N_Mesh(Nloc)%VolToSideIJKA, N_Mesh(Nloc)%FS2M)
  END DO ! Nloc = Nmin, Nmax

END IF ! meshMode.GT.0

CALL InitpAdaption() ! Calls Build_N_DG_Mapping() which builds N_DG_Mapping()

#if USE_LOADBALANCE
IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
!IF (PerformLoadBalance) THEN
  CALL ExchangeVolMesh() !  Allocates N_VolMesh(1:nElems) after communication of Elem_xGP
ELSE
  ! Build the coordinates of the solution gauss points in the volume
#endif /*USE_LOADBALANCE*/
  ! Build Elem_xGP
  ALLOCATE(N_VolMesh(1:nElems))
  ALLOCATE(N_VolMesh2(1:nElems))
  CALL BuildElem_xGP(NodeCoords) ! Builds N_VolMesh(iElem)%Elem_xGP, requires N_DG_Mapping() for Nloc
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

CALL DG_ProlongDGElemsToFace() ! Builds DG_Elems_master and ,DG_Elems_slave, requires SideToElem()

! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----

IF (ABS(meshMode).GT.1) THEN
#if USE_LOADBALANCE
  IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  !IF (PerformLoadBalance) THEN
    ! Shift metric arrays during load balance
    CALL ExchangeMetrics()
  ELSE
#endif /*USE_LOADBALANCE*/

    NGeoRef=3*NGeo ! build jacobian at higher degree
    ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

    ! volume data
    DO iElem = 1, nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      ALLOCATE(N_VolMesh(iElem)%XCL_N(       3,  0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%Metrics_fTilde(3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%Metrics_gTilde(3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%Metrics_hTilde(3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%sJ            (  0:Nloc,0:Nloc,0:Nloc))
    END DO ! iElem = 1, nElems

    ! volume data
    DO iElem = 1, nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      ALLOCATE(N_VolMesh2(iElem)%dXCL_N(     3,3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh2(iElem)%JaCL_N(     3,3,0:Nloc,0:Nloc,0:Nloc))
    END DO ! iElem = 1, nElems

    ! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
    ! assign 1/detJ (sJ)
    ! assign normal and tangential vectors and surfElems on faces

    crossProductMetrics=GETLOGICAL('crossProductMetrics')
#if USE_HDG
    IF(Symmetry%Axisymmetric.AND..NOT.crossProductMetrics) THEN
      crossProductMetrics = .TRUE.
      CALL PrintOption('WARNING: axisymmetric HDG simulations require crossProductMetrics','INFO',LogOpt=crossProductMetrics)
    END IF
#endif /*USE_HDG*/
    LBWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."

    ! get XCL_NGeo
    ALLOCATE(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
    XCL_NGeo = 0.
#if defined(PARTICLES)
    ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
    dXCL_NGeo = 0.
#endif /*defined(PARTICLES)*/

#ifdef PARTICLES
    CALL CalcMetrics(XCL_NGeo_Out=XCL_NGeo,dXCL_NGeo_Out=dXCL_NGeo)
#else
    CALL CalcMetrics(XCL_NGeo_Out=XCL_NGeo)
#endif
#if USE_LOADBALANCE
  END IF
#endif /*USE_LOADBALANCE*/

  ! surface data
  ALLOCATE(N_SurfMesh(1:nSides))
  DO iSide = 1, nSides

    ! Allocate with max. polynomial degree of the two master-slave sides
    Nloc     = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
    NSideMin = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
    N_SurfMesh(iSide)%NSideMin = NSideMin
    ALLOCATE(N_SurfMesh(iSide)%Face_xGP (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%NormVec  (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%TangVec1 (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%TangVec2 (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%SurfElem (  0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%SurfElemMin(0:NSideMin,0:NSideMin))
    N_SurfMesh(iSide)%Face_xGP = 0.
    N_SurfMesh(iSide)%NormVec  = 0.
    N_SurfMesh(iSide)%TangVec1 = 0.
    N_SurfMesh(iSide)%TangVec2 = 0.
    N_SurfMesh(iSide)%SurfElem = 0.
    N_SurfMesh(iSide)%SurfElemMin= 0.
  END DO ! iSide = 1, nSides

#if USE_HDG
  ! Mortars: Get minimum of all sides the Mortar interface
  DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
    nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
    locSide  = MortarType(2,MortarSideID)
    NSideMin = N_SurfMesh(MortarSideID)%NSideMin
    DO iMortar = 1,nMortars
      SideID = MortarInfo(MI_SIDEID,iMortar,locSide) !small SideID
      NSideMin = MIN(NSideMin,N_SurfMesh(SideID)%NSideMin)
    END DO !iMortar

    N_SurfMesh(MortarSideID)%NSideMin = NSideMin
    DO iMortar = 1,nMortars
      SideID = MortarInfo(MI_SIDEID,iMortar,locSide) !small SideID
      N_SurfMesh(SideID)%NSideMin = NSideMin
    END DO !iMortar
  END DO !MortarSideID
#endif /*USE_HDG*/

  ! Due to possible load balance, this is done outside of CalcMetrics() now
  DO iElem = 1, nElems
    CALL CalcSurfMetrics(iElem)
  END DO ! iElem = 1, nElems

#if defined(PARTICLES) || USE_HDG
  ! Compute element bary and element radius for processor-local elements (without halo region)
  ALLOCATE(ElemBaryNGeo(1:3,1:nElems))
  CALL BuildElementOrigin()
#endif /*defined(PARTICLES) || USE_HDG*/

#ifndef PARTICLES
  ! dealloacte pointers
  LBWRITE(UNIT_stdOut,'(A)') "InitMesh: NOW CALLING deleteMeshPointer..."
  CALL deleteMeshPointer()
#endif

  ! Initialize element volumes and characteristic lengths
  CALL InitElemVolumes()

#ifndef PARTICLES
  IF(meshMode.GT.1) DEALLOCATE(NodeCoords)
#endif

  IF((ABS(meshMode).NE.3).AND.(meshMode.GT.1))THEN
#ifdef PARTICLES
    IF(RadialWeighting%DoRadialWeighting) THEN
      usevMPF = .TRUE.
    ELSE
      usevMPF = GETLOGICAL('Part-vMPF','.FALSE.')
    END IF
#endif /* PARTICLES */
  END IF ! meshMode.NE.3
END IF ! meshMode.GT.1

IF(CalcMeshInfo)THEN
  CALL AddToElemData(ElementOut,'myRank',IntScalar=myRank)
  !#ifdef PARTICLES
  ALLOCATE(ElemGlobalID(1:nElems))
  DO iElem=1,nElems
    ElemGlobalID(iElem)=OffsetElem+iElem
  END DO ! iElem=1,nElems
  CALL AddToElemData(ElementOut,'ElemID',LongIntArray=ElemGlobalID)
  !#endif /*PARTICLES*/
END IF

#if USE_HDG && USE_LOADBALANCE
IF (ABS(meshMode).GT.0) CALL BuildSideToNonUniqueGlobalSide() ! requires ElemInfo
#endif /*USE_HDG && USE_LOADBALANCE*/
!DEALLOCATE(ElemInfo,SideInfo)
DEALLOCATE(SideInfo)
IF(readFEMconnectivity)THEN
  SDEALLOCATE(EdgeInfo)
  SDEALLOCATE(VertexInfo)
  SDEALLOCATE(EdgeConnectInfo)
  SDEALLOCATE(VertexConnectInfo)
END IF

MeshInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


!===================================================================================================================================
!> Get local N for each element and side
!===================================================================================================================================
SUBROUTINE InitpAdaption()
! MODULES
USE MOD_Globals
USE MOD_PreProc
!USE MOD_DG_Vars            ,ONLY: DG_Elems_master,DG_Elems_slave
USE MOD_DG_Vars            ,ONLY: N_DG,pAdaptionType,pAdaptionBCLevel,NDGAllocationIsDone
USE MOD_IO_HDF5            ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars          ,ONLY: nElems,SideToElem,nBCSides,Boundarytype,BC,readFEMconnectivity,NodeCoords
USE MOD_ReadInTools        ,ONLY: GETINTFROMSTR
USE MOD_Interpolation_Vars ,ONLY: NMax,NMin
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,BCSideID,BCType
REAL    :: RandVal,x
LOGICAL :: SetBCElemsToNMax
!===================================================================================================================================
! Set defaults
SetBCElemsToNMax = .FALSE. ! Initialize
pAdaptionBCLevel = -1
NDGAllocationIsDone = .FALSE.

! Read p-adaption specific input data
pAdaptionType = GETINTFROMSTR('pAdaptionType')

! Allocate arrays and initialize local polynomial degree
! This happens here because nElems is determined here and N_DG is required below for the mesh initialisation
ALLOCATE(N_DG(1:nElems))
N_DG = -1
!CALL AddToElemData(ElementOut,'Nloc',IntArray=N_DG_Mapping(2,1+offSetElem:nElems+offSetElem)) ! TODO: Why does this not work?
! Add array containing the local polynomial degree to the hdf5 output
CALL AddToElemData(ElementOut,'Nloc',IntArray=N_DG)

SELECT CASE(pAdaptionType)
CASE(PRM_P_ADAPTION_ZERO)
  N_DG = PP_N ! By default, the initial degree is set to PP_N
CASE(PRM_P_ADAPTION_RDN) ! Random between NMin and NMax
  DO iElem=1,nElems
    CALL RANDOM_NUMBER(RandVal)
    N_DG(iElem) = NMin + INT(RandVal*(NMax-NMin+1))
  END DO
CASE(PRM_P_ADAPTION_NPB) ! Non-periodic BCs are set to NMax
  ! Get depth of increased polynomial degree
  pAdaptionBCLevel = GETINTFROMSTR('pAdaptionBCLevel')
  N_DG = Nmin ! By default, the initial degree is set to Nmin
  SetBCElemsToNMax = .TRUE.
CASE(PRM_P_ADAPTION_HH) ! Elements in the lower half domain in x-direction are set to NMin and the upper half are set to NMax
  DO iElem=1,nElems
    x = (MAXVAL(NodeCoords(1,:,:,:,iElem))+MINVAL(NodeCoords(1,:,:,:,iElem)))/2.0
    IF(x.GT.0.0)THEN
      N_DG(iElem) = NMax
    ELSE
      N_DG(iElem) = NMin
    END IF ! x.GT.0.0
  END DO
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown pAdaptionType!' ,IntInfo=pAdaptionType)
END SELECT

! Check if BC elements are to be set to a higher polynomial degree
IF(ABS(pAdaptionBCLevel).GT.1.OR.&
  (ABS(pAdaptionBCLevel).GT.0.AND.readFEMconnectivity))THEN
  CALL SetpAdaptionBCLevel()
ELSE
  ! Check if all BC elements are set to Nmax
  IF(SetBCElemsToNMax)THEN
    DO BCSideID=1,nBCSides
      BCType=Boundarytype(BC(BCSideID),BC_TYPE)
      IF(BCType.EQ.1) CYCLE ! Skip periodic sides
      iElem       = SideToElem(S2E_ELEM_ID,BCSideID)
      N_DG(iElem) = NMax
    END DO ! BCSideID=1,nBCSides
  END IF ! SetBCElemsToNMax
END IF

! Sanity check
DO iElem=1,nElems
  IF((N_DG(iElem).LT.NMin).OR.(N_DG(iElem).GT.NMax))THEN
    IPWRITE(*,*) "iElem       = ", iElem
    IPWRITE(*,*) "N_DG(iElem) = ", N_DG(iElem)
    IPWRITE(*,*) "NMin        = ", NMin
    IPWRITE(*,*) "NMax        = ", NMax
  END IF
  IF(N_DG(iElem).LT.NMin) CALL abort(__STAMP__,'N_DG(iElem)<NMin')
  IF(N_DG(iElem).GT.NMax) CALL abort(__STAMP__,'N_DG(iElem)>NMax')
END DO

! Initialize element containers
CALL Build_N_DG_Mapping()

END SUBROUTINE InitpAdaption


!===================================================================================================================================
!> Set Nloc of BC elements to a higher polynomial degree than in the volume
!>
!> 1. Loop over the processor-local elements
!> 2. Loop over the corner vertices of the element
!> 3. Use VertexConnectInfo to get the neighbour element index and vertex for all possible connection (also periodic)
!> 4. Check the sides of connected to the neighbour node and find out if the side is a BC side
!> 5. For further layers, loop over the elements again and check for already marked elements instead of the sides
!===================================================================================================================================
SUBROUTINE SetpAdaptionBCLevel()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: nFEMEdges, nFEMVertices, readFEMconnectivity, offsetElem, nElems
USE MOD_Mesh_Vars          ,ONLY: ElemInfo,SideInfo,EdgeInfo,EdgeConnectInfo,VertexInfo,VertexConnectInfo
USE MOD_Mesh_Vars          ,ONLY: SideToElem,BoundaryType,BC
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_DG_Vars            ,ONLY: N_DG,pAdaptionType,pAdaptionBCLevel,N_DG_Mapping
USE MOD_Interpolation_Vars ,ONLY: NMax,NMin
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping_Shared_Win
USE MOD_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,BCSideID,BCType,NonUniqueGlobalSideID,iGlobalElemID,BCIndex,ElemType,OffsetCounter
INTEGER :: iVertex,iVertexConnect,GlobalNbElemID,GlobalNbLocVertexID,LocSideList(3),iLocSideList,SideID,iLocSide,localSideID
INTEGER :: FirstSideInd,LastSideInd,FirstElemInd,LastElemInd
INTEGER :: FirstEdgeInd,LastEdgeInd,FirstEdgeConnectInd,LastEdgeConnectInd
INTEGER :: nVertexIDs,offsetVertexID,nVertexConnectIDs,offsetVertexConnectID
INTEGER :: FirstVertexInd,LastVertexInd,FirstVertexConnectInd,LastVertexConnectInd
!===================================================================================================================================
! Do not re-allocate during load balance here as it is communicated between the processors
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

! Sanity check: This routine requires FEM connectivity
IF(.NOT.readFEMconnectivity) CALL abort(__STAMP__,'Error in p-adaption init: readFEMconnectivity=T is required!')

! Element index
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems

! Loop over the process-local global elements indices
GlobalElemIDLoop: DO iGlobalElemID = FirstElemInd, LastElemInd
  iElem = iGlobalElemID - offsetElem
  ElemType = ElemInfo(ELEM_TYPE,iGlobalElemID)
  ! Sanity check: currently only hexahedral elements are implemented
  SELECT CASE(ElemType)
  CASE(108,118,208)
    ! Hexahedral elements
  CASE DEFAULT
    CALL abort(__STAMP__,'Element type not implemented: ElemType =',IntInfoOpt=ElemType)
  END SELECT
  ! Get local VertexInfo of current element
  FirstVertexInd = ElemInfo(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1
  LastVertexInd  = ElemInfo(ELEM_LASTVERTEXIND,iGlobalElemID)
  ! Get local vertex connectivity
  FirstVertexConnectInd = VertexInfo(VERTEX_FIRSTCONNECTIND,FirstVertexInd)+1
  LastVertexConnectInd  = VertexInfo(VERTEX_LASTCONNECTIND,LastVertexInd)
  VertexConnectLoop: DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
    ! Check if current element has already been flagged
    IF(N_DG(iElem).EQ.NMax) EXIT VertexConnectLoop
    ! Get neighbour infos
    GlobalNbElemID      = ABS(VertexConnectInfo(VERTEXCONNECT_NBELEMID   ,iVertexConnect))
    GlobalNbLocVertexID = VertexConnectInfo(VERTEXCONNECT_NBLOCNODEID,iVertexConnect)
    ! Set sides depending on the element type: Only implemented for Hexahedral elements
    CALL GetLocSideList(ElemType,GlobalNbLocVertexID,LocSideList)
    LocSideListLoop: DO iLocSideList = 1, 3
      ! Check if current element has already been flagged
      IF(N_DG(iElem).EQ.NMax) EXIT VertexConnectLoop
      iLocSide = LocSideList(iLocSideList)
      NonUniqueGlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalNbElemID) + iLocSide
      BCIndex = SideInfo_Shared(SIDE_BCID,NonUniqueGlobalSideID)
      IF(BCIndex.LE.0) CYCLE LocSideListLoop ! Skip inner sides
      BCType = BoundaryType(BCIndex,BC_TYPE)
      IF(BCType.LE.1) CYCLE LocSideListLoop ! Skip periodic sides
      IF(BCType.EQ.10) CYCLE LocSideListLoop ! Skip Neumann sides
      IF(pAdaptionBCLevel.EQ.-1)THEN
        N_DG(iElem) = NMin+1
      ELSE
        N_DG(iElem) = NMax
      END IF ! pAdaptionBCLevel.EQ.-1
    END DO LocSideListLoop ! iLocSideList = 1, 3
  END DO VertexConnectLoop ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
END DO GlobalElemIDLoop ! iGlobalElemID = FirstElemInd, LastElemInd


! For further layers, loop over the elements again and check for already marked elements instead of the sides
IF(ABS(pAdaptionBCLevel).GT.1)THEN
  ! Allocate the shared memory container and associate pointer: N_DG_Mapping
  CALL Allocate_N_DG_Mapping()

  ! Set Nloc in N_DG_Mapping
  CALL Set_N_DG_Mapping(OffsetCounter)

  ! Loop over the process-local global elements indices
  iGlobalElemID_loop: DO iGlobalElemID = FirstElemInd, LastElemInd
    iElem = iGlobalElemID - offsetElem
    IF(N_DG(iElem).GT.NMin) CYCLE iGlobalElemID_loop
    ! Get local VertexInfo of current element
    FirstVertexInd = ElemInfo(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1
    LastVertexInd  = ElemInfo(ELEM_LASTVERTEXIND,iGlobalElemID)
    ! Get local vertex connectivity
    FirstVertexConnectInd = VertexInfo(VERTEX_FIRSTCONNECTIND,FirstVertexInd)+1
    LastVertexConnectInd  = VertexInfo(VERTEX_LASTCONNECTIND,LastVertexInd)
    iVertexConnect_loop: DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
      ! Check if current element has already been flagged
      IF(N_DG(iElem).GT.NMin) EXIT iVertexConnect_loop
      ! Get neighbour infos
      GlobalNbElemID = ABS(VertexConnectInfo(VERTEXCONNECT_NBELEMID   ,iVertexConnect))
      IF(N_DG_Mapping(2,GlobalNbElemID).EQ.NMax)THEN
        IF(pAdaptionBCLevel.EQ.-2)THEN
          N_DG(iElem) = NMin+1
        ELSE
          N_DG(iElem) = NMax
        END IF ! pAdaptionBCLevel.EQ.-2
      END IF ! N_DG_Mapping(2,GlobalNbElemID).EQ.NMax
    END DO iVertexConnect_loop ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
  END DO iGlobalElemID_loop ! iGlobalElemID = FirstElemInd, LastElemInd

#if USE_MPI
  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/
END IF ! pAdaptionBCLevel.GT.1


END SUBROUTINE SetpAdaptionBCLevel

!===================================================================================================================================
!> Returns a list of sides (depending on the CGNS ordering) for a given local corner node index
!===================================================================================================================================
SUBROUTINE GetLocSideList(ElemType,iLocNode,LocSideList)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)  :: ElemType
INTEGER, INTENT(IN)  :: iLocNode
INTEGER, INTENT(OUT) :: LocSideList(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Set sides depending on the element type: Only implemented for Hexahedral elements
SELECT CASE(ElemType)
CASE(108,118,208)
  ! Hexahedral elements
  SELECT CASE(iLocNode)
  CASE(1)
    LocSideList=(/1,2,5/)
  CASE(2)
    LocSideList=(/1,2,3/)
  CASE(3)
    LocSideList=(/1,3,4/)
  CASE(4)
    LocSideList=(/1,4,5/)
  CASE(5)
    LocSideList=(/2,5,6/)
  CASE(6)
    LocSideList=(/2,3,6/)
  CASE(7)
    LocSideList=(/3,4,6/)
  CASE(8)
    LocSideList=(/4,5,6/)
  CASE DEFAULT
    CALL abort(__STAMP__,'Wrong iLocNode',IntInfoOpt=iLocNode)
  END SELECT
CASE DEFAULT
  CALL abort(__STAMP__,'Element type not implemented: ElemType =',IntInfoOpt=ElemType)
END SELECT
END SUBROUTINE GetLocSideList


!===================================================================================================================================
!> Allocate the shared memory container N_DG_Mapping_Shared
!===================================================================================================================================
SUBROUTINE Allocate_N_DG_Mapping()
! MODULES
USE MOD_Globals
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping,N_DG_Mapping_Shared,N_DG,NDGAllocationIsDone
USE MOD_Mesh_Vars       ,ONLY: nElems,offSetElem,nGlobalElems
#if USE_MPI
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping_Shared_Win
USE MOD_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED, myComputeNodeRank
USE MOD_MPI_Shared
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
CALL Allocate_Shared((/3,nGlobalElems/),N_DG_Mapping_Shared_Win,N_DG_Mapping_Shared)
CALL MPI_WIN_LOCK_ALL(0,N_DG_Mapping_Shared_Win,IERROR)
N_DG_Mapping => N_DG_Mapping_Shared
IF (myComputeNodeRank.EQ.0) N_DG_Mapping = 0
CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(N_DG_Mapping(3,nElems))
N_DG_Mapping = 0
#endif /*USE_MPI*/
NDGAllocationIsDone = .TRUE.
END SUBROUTINE Allocate_N_DG_Mapping

!===================================================================================================================================
!> Loop over all CN elements and set N_DG_Mapping(1:2,iElem+offSetElem). N_DG_Mapping(2,:) is only correct for single-process
!> execution as the communication between the processes is done later.
!===================================================================================================================================
SUBROUTINE Set_N_DG_Mapping(OffsetCounter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY: nElems,offSetElem
USE MOD_DG_Vars   ,ONLY: N_DG_Mapping,N_DG
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: OffsetCounter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: locDofs,iElem,iCNElem,Nloc
!===================================================================================================================================
OffsetCounter = 0
! Loop all CN elements
DO iCNElem = 1+offSetElem,nElems+offSetElem
  iElem = iCNElem-offSetElem
  Nloc = N_DG(iElem)
  locDofs = (Nloc+1)**3
  N_DG_Mapping(2,iCNElem) = Nloc
  N_DG_Mapping(1,iCNElem) = OffsetCounter
  OffsetCounter = OffsetCounter + locDofs
END DO ! iCNElem
END SUBROUTINE Set_N_DG_Mapping


!===================================================================================================================================
!> Create shared memory array N_DG_Mapping containing the global element information
!>   N_DG_Mapping(1,nElems+offSetElem): DOF offset
!>   N_DG_Mapping(2,nElems+offSetElem): element polynomial degree Nloc
!===================================================================================================================================
SUBROUTINE Build_N_DG_Mapping()
! MODULES
USE MOD_Globals
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping,displsDofs, recvcountDofs, N_DG_Mapping_Shared, nDofsMapping, N_DG
USE MOD_DG_Vars            ,ONLY: NDGAllocationIsDone
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem,nGlobalElems
#if USE_MPI
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping_Shared_Win
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_LEADERS_SHARED, MPI_COMM_SHARED, myComputeNodeRank, myleadergrouprank
USE MOD_MPI_Shared_Vars    ,ONLY: nLeaderGroupProcs,nComputeNodeProcessors
USE MOD_Particle_Mesh_Vars ,ONLY: offsetComputeNodeElem
USE MOD_MPI_Shared
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: OffsetCounter,OffsetN_DG_Mapping
#if USE_MPI
INTEGER             :: iProc
INTEGER             :: sendbuf,recvbuf
#endif
!===================================================================================================================================
! Do not re-allocate during load balance here as it is communicated between the processors
#if USE_LOADBALANCE
IF(PerformLoadBalance)THEN
#endif /*USE_LOADBALANCE*/

  ! N_DG_Mapping is already set
  !OffsetCounter = N_DG_Mapping(1,nElems+offSetElem) + (N_DG_Mapping(2,nElems+offSetElem)+1)**3

#if USE_LOADBALANCE
ELSE
#endif /*USE_LOADBALANCE*/
  IF(.NOT.NDGAllocationIsDone)THEN
    ! Allocate the shared memory container and associate pointer: N_DG_Mapping
    CALL Allocate_N_DG_Mapping()
  END IF ! .NOT.NDGAllocationIsDone

  ! Set Nloc in N_DG_Mapping
  CALL Set_N_DG_Mapping(OffsetCounter)

#if USE_MPI
  ! Set the correct offsets in N_DG_Mapping(1,:) for parallel simulation
  sendbuf = OffsetCounter
  recvbuf = 0
  CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
  OffsetN_DG_Mapping   = recvbuf
  ! The last process knows CN total number of connected CN elements
  sendbuf = OffsetN_DG_Mapping + OffsetCounter
  CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nProcessors-1,MPI_COMM_PICLAS,iError)
  nDofsMapping = sendbuf

  N_DG_Mapping(1,1+offSetElem:nElems+offSetElem) = N_DG_Mapping(1,1+offSetElem:nElems+offSetElem) + OffsetN_DG_Mapping
  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win,MPI_COMM_SHARED)

  ! Communication between nodes
  IF (nComputeNodeProcessors.NE.nProcessors.AND.myComputeNodeRank.EQ.0) THEN
    ! Arrays for the compute node to hold the elem offsets
    ALLOCATE(displsDofs(   0:nLeaderGroupProcs-1), recvcountDofs(0:nLeaderGroupProcs-1))
    displsDofs(myLeaderGroupRank) = offsetComputeNodeElem
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsDofs,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountDofs(iProc-1) = displsDofs(iProc)-displsDofs(iProc-1)
    END DO
    recvcountDofs(nLeaderGroupProcs-1) = nGlobalElems - displsDofs(nLeaderGroupProcs-1)

    CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
        , 0                             &
        , MPI_DATATYPE_NULL             &
        , N_DG_Mapping               &
        , 3*recvcountDofs   &
        , 3*displsDofs      &
        , MPI_INTEGER          &
        , MPI_COMM_LEADERS_SHARED       &
        , IERROR)

    displsDofs(myLeaderGroupRank) = N_DG_Mapping(1,1+offSetElem)
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsDofs,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountDofs(iProc-1) = displsDofs(iProc)-displsDofs(iProc-1)
    END DO
    recvcountDofs(nLeaderGroupProcs-1) = nDofsMapping - displsDofs(nLeaderGroupProcs-1)
  END IF

  CALL BARRIER_AND_SYNC(N_DG_Mapping_Shared_Win ,MPI_COMM_SHARED)

#else
  OffsetN_DG_Mapping = 0
  nDofsMapping = OffsetCounter
#endif /*USE_MPI*/

#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

END SUBROUTINE Build_N_DG_Mapping


!==================================================================================================================================
!> Set DG_Elems_slave and DG_Elems_master information
!==================================================================================================================================
SUBROUTINE DG_ProlongDGElemsToFace()
! MODULES
USE MOD_PreProc
USE MOD_GLobals
USE MOD_DG_Vars   ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave,N_DG_Mapping!,pAdaptionType
USE MOD_Mesh_Vars ,ONLY: SideToElem,nSides,nBCSides, offSetElem

USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide,MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
#if USE_MPI
USE MOD_MPI       ,ONLY: StartExchange_DG_Elems,FinishExchangeMPIData
USE MOD_MPI_Vars  ,ONLY: DataSizeSideSend,DataSizeSideRec,nNbProcs,nMPISides_rec,nMPISides_send,OffsetMPISides_rec
USE MOD_MPI_Vars  ,ONLY: OffsetMPISides_send
USE MOD_MPI_Vars  ,ONLY: DataSizeSurfSendMax,DataSizeSurfRecMax, DataSizeSurfSendMin,DataSizeSurfRecMin
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iSide,ElemID,nbElemID,nMortars, locSide,SideID, iMortar,flip!, iElem
#if USE_MPI
INTEGER                      :: iNbProc,Nloc
INTEGER, DIMENSION(nNbProcs) :: RecRequest_U,SendRequest_U,RecRequest_U2,SendRequest_U2
#endif /*USE_MPI*/
!==================================================================================================================================
! Side containers
ALLOCATE(DG_Elems_master(1:nSides))
ALLOCATE(DG_Elems_slave (1:nSides))
! Initialize with element-local N
!IF(pAdaptionType.EQ.0)THEN
  DG_Elems_master = PP_N
  DG_Elems_slave  = PP_N
!ELSE
!  DO iSide = 1, nSides
!    iElem = SideToElem(S2E_ELEM_ID,iSide)
!    DG_Elems_master(iSide) = N_DG_Mapping(2,iElem+offSetElem)
!    DG_Elems_slave(iSide)  = N_DG_Mapping(2,iElem+offSetElem)
!  END DO ! iSide = 1, nSides
!END IF ! pAdaptionType.EQ.0

! set information which polynomial degree elements adjacent to a side have
DO iSide = 1,nSides
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,iSide)
  ! Master sides
  IF(ElemID  .GT.0) DG_Elems_master(iSide) = N_DG_Mapping(2,ElemID+offSetElem)
  ! Slave side (ElemID,locSide and flip =-1 if not existing)
  IF(nbElemID.GT.0) DG_Elems_slave( iSide) = N_DG_Mapping(2,nbElemID+offSetElem)
  ! Boundaries
  IF(iSide.LE.nBCSides) DG_Elems_slave( iSide) = DG_Elems_master(iSide)
END DO

DO iSide = firstMortarInnerSide,lastMortarInnerSide
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nMortars=MERGE(4,2,MortarType(1,iSide).EQ.1)
  locSide=MortarType(2,iSide)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    SELECT CASE(flip)
     CASE(0) ! master side
       DG_Elems_master(SideID) = N_DG_Mapping(2,ElemID+offSetElem)
     CASE(1:4) ! slave side
       DG_Elems_slave(SideID) = N_DG_Mapping(2,ElemID+offSetElem)
    END SELECT !f
  END DO
END DO

#if USE_MPI
DO iSide = firstMortarMPISide,lastMortarMPISide
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nMortars=MERGE(4,2,MortarType(1,iSide).EQ.1)
  locSide=MortarType(2,iSide)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    SELECT CASE(flip)
     CASE(0) ! master side
       DG_Elems_master(SideID) = N_DG_Mapping(2,ElemID+offSetElem)
     CASE(1:4) ! slave side
       DG_Elems_slave(SideID) = N_DG_Mapping(2,ElemID+offSetElem)
    END SELECT !f
  END DO
END DO
! Exchange element local polynomial degree (N_LOC)
CALL StartExchange_DG_Elems(DG_Elems_slave ,1,nSides,SendRequest_U ,RecRequest_U ,SendID=2)  ! RECEIVE MINE, SEND YOUR / DG_Elems_slave:  slave  -> master
CALL StartExchange_DG_Elems(DG_Elems_master,1,nSides,SendRequest_U2,RecRequest_U2,SendID=1)  ! RECEIVE YOUR, SEND MINE / DG_Elems_master: master -> slave
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U ,RecRequest_U ,SendID=2) ! Send YOUR - receive MINE
CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=1) ! Send YOUR - receive MINE
! Initialize the send/rec face sizes for master/slave communication
ALLOCATE(DataSizeSideSend(nNbProcs,2)   , DataSizeSideRec(nNbProcs,2))
! Min/Max is only required for HDG
ALLOCATE(DataSizeSurfSendMax(nNbProcs,2), DataSizeSurfRecMax(nNbProcs,2))
ALLOCATE(DataSizeSurfSendMin(nNbProcs,2), DataSizeSurfRecMin(nNbProcs,2))
DataSizeSideSend=0
DataSizeSideRec=0
DataSizeSurfSendMax =0; DataSizeSurfRecMax = 0; DataSizeSurfSendMin =0; DataSizeSurfRecMin = 0
DO iNbProc = 1, nNbProcs

  ! 1: Set number of sides and offset for SEND MINE - RECEIVE YOUR case
  IF(nMPISides_rec(iNbProc,1).GT.0) THEN
    DO iSide = OffsetMPISides_rec(iNbProc-1,1)+1, OffsetMPISides_rec(iNbProc,1)
      Nloc = DG_Elems_slave(iSide) ! polynomial degree of the sending slave side
      DataSizeSideRec(iNbProc,1) = DataSizeSideRec(iNbProc,1)  + (Nloc+1)**2
      Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
      DataSizeSurfRecMax(iNbProc,1) = DataSizeSurfRecMax(iNbProc,1)  + (Nloc+1)**2
      Nloc = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
      DataSizeSurfRecMin(iNbProc,1) = DataSizeSurfRecMin(iNbProc,1)  + (Nloc+1)**2
    END DO ! iSide = 1, nMPISides_rec(iNbProc,1)
  END IF

  IF(nMPISides_send(iNbProc,1).GT.0) THEN
    DO iSide = OffsetMPISides_send(iNbProc-1,1)+1, OffsetMPISides_send(iNbProc,1)
      Nloc = DG_Elems_slave(iSide) ! polynomial degree of the sending slave side
      DataSizeSideSend(iNbProc,1) = DataSizeSideSend(iNbProc,1)  + (Nloc+1)**2
      Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
      DataSizeSurfSendMax(iNbProc,1) = DataSizeSurfSendMax(iNbProc,1)  + (Nloc+1)**2
      Nloc = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
      DataSizeSurfSendMin(iNbProc,1) = DataSizeSurfSendMin(iNbProc,1)  + (Nloc+1)**2
    END DO ! iSide = 1, nMPISides_rec(iNbProc,1)
  END IF

  ! Important: Keep symmetry
  DataSizeSideSend(iNbProc,2) = DataSizeSideRec( iNbProc,1)
  DataSizeSideRec( iNbProc,2) = DataSizeSideSend(iNbProc,1)
  DataSizeSurfSendMax(iNbProc,2) = DataSizeSurfRecMax( iNbProc,1)
  DataSizeSurfRecMax( iNbProc,2) = DataSizeSurfSendMax(iNbProc,1)
  DataSizeSurfSendMin(iNbProc,2) = DataSizeSurfRecMin( iNbProc,1)
  DataSizeSurfRecMin( iNbProc,2) = DataSizeSurfSendMin(iNbProc,1)
END DO ! iNbProc = 1, nNbProcs

#endif /*USE_MPI*/
END SUBROUTINE DG_ProlongDGElemsToFace


SUBROUTINE GetMeshMinMaxBoundaries()
!============================================================================================================================
! Get min and max coordinates of the face xGP and store them in "xyzMinMax" array
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Mesh_Vars ,ONLY: xyzMinMax,GetMeshMinMaxBoundariesIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
#if USE_MPI
REAL              :: xmin_loc(3),xmax_loc(3) ! for sending via MPI
REAL              :: xmin(3),xmax(3) ! for receiving via MPI
#endif /*USE_MPI*/
!============================================================================================================================
! check if already called
IF(GetMeshMinMaxBoundariesIsDone) RETURN

! Get global min/max for all directions x, y and z
#if USE_MPI

   ! Map global min/max values from xyzMinMax(1:6) that were determined in metric.f90 for each processor
   xmin_loc(1) = xyzMinMax(1)
   xmin_loc(2) = xyzMinMax(3)
   xmin_loc(3) = xyzMinMax(5)

   xmax_loc(1) = xyzMinMax(2)
   xmax_loc(2) = xyzMinMax(4)
   xmax_loc(3) = xyzMinMax(6)

   ! Find global min
   CALL MPI_ALLREDUCE(xmin_loc(1:3),xmin(1:3), 3, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_PICLAS, IERROR)
   ! Find global max
   CALL MPI_ALLREDUCE(xmax_loc(1:3),xmax(1:3), 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_PICLAS, IERROR)

   ! Map global min/max values to xyzMinMax(1:6)
   xyzMinMax(1) = xmin(1)
   xyzMinMax(3) = xmin(2)
   xyzMinMax(5) = xmin(3)

   xyzMinMax(2) = xmax(1)
   xyzMinMax(4) = xmax(2)
   xyzMinMax(6) = xmax(3)
#else
   ! Already doe in metric.f90
#endif /*USE_MPI*/

! don't call twice
GetMeshMinMaxBoundariesIsDone=.TRUE.

END SUBROUTINE GetMeshMinMaxBoundaries


#if defined(PARTICLES) || USE_HDG
SUBROUTINE BuildElementOrigin()
!================================================================================================================================
! compute the element origin at xi=(0,0,0)^T and set it as ElemBaryNGeo
!================================================================================================================================
USE MOD_Globals!,                  ONLY:CROSS
USE MOD_Preproc
USE MOD_Mesh_Vars,                ONLY:NGeo,XCL_NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Mesh_Vars,                ONLY:ElemBaryNGeo
USE MOD_Basis,                    ONLY:LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem,i,j,k
REAL                    :: XPos(3),buf
REAL                    :: Lag(1:3,0:NGeo)
!================================================================================================================================

ElemBaryNGeo=0.
DO iElem=1,PP_nElems
  ! evaluate the polynomial at origin: Xi=(/0.0,0.0,0.0/)
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  xPos=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buf=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*buf
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo
  ElemBaryNGeo(:,iElem)=xPos
END DO ! iElem

END SUBROUTINE BuildElementOrigin
#endif /*defined(PARTICLES) || USE_HDG*/


SUBROUTINE InitElemVolumes()
!===================================================================================================================================
! Calculate Element volumes for later use in particle routines
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals            ,ONLY: UNIT_StdOut,abort
USE MOD_Interpolation_Vars ,ONLY: N_Inter
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_Mesh_Vars          ,ONLY: nElems,N_VolMesh, offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: LocalVolume,MeshVolume
USE MOD_Particle_Mesh_Vars ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
USE MOD_ReadInTools
#if USE_MPI
USE MPI
USE MOD_MPI_Shared
USE MOD_Globals            ,ONLY: IERROR,MPIRoot
#ifdef PARTICLES
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars ,ONLY: nComputeNodeElems,offsetComputeNodeElem
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED,myComputeNodeRank,MPI_COMM_LEADERS_SHARED
USE MOD_Particle_Mesh_Vars ,ONLY: ElemVolume_Shared_Win,ElemCharLength_Shared_Win
#else
USE MOD_Globals            ,ONLY: MPI_COMM_PICLAS
#endif /*PARTICLES*/
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,CNElemID
INTEGER                         :: i,j,k,Nloc
INTEGER                         :: offsetElemCNProc
#if USE_MPI
#ifdef PARTICLES
REAL                            :: CNVolume                       ! Total CN volume
#endif /*PARTICLES*/
#endif /*USE_MPI*/
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT ELEMENT GEOMETRY INFORMATION ...'
! Sanity check
IF(.NOT.ALLOCATED(N_Inter)) CALL abort(__STAMP__,'Error in InitElemVolumes(): N_Inter is not allocated')

#if USE_MPI && defined(PARTICLES)
! J_N is only built for local DG elements. Therefore, array is only filled for elements on the same compute node
offsetElemCNProc = offsetElem - offsetComputeNodeElem
#else
offsetElemCNProc = 0
#endif  /*USE_MPI && defined(PARTICLES)*/

! In case of MPI=ON and PARTICLES=OFF, no shared array is created and all arrays are processor-local
#if USE_MPI && defined(PARTICLES)
CALL Allocate_Shared((/nComputeNodeElems/),ElemVolume_Shared_Win,ElemVolume_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemVolume_Shared_Win,IERROR)
CALL Allocate_Shared((/nComputeNodeElems/),ElemCharLength_Shared_Win,ElemCharLength_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCharLength_Shared_Win,IERROR)

! Only root nullifies
IF (myComputeNodeRank.EQ.0) THEN
  ElemVolume_Shared(:)     = 0.
  ElemCharLength_Shared(:) = 0.
END IF
CALL BARRIER_AND_SYNC(ElemVolume_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemCharLength_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(ElemVolume_Shared(nElems))
ALLOCATE(ElemCharLength_Shared(nElems))
ElemVolume_Shared(:)     = 0.
ElemCharLength_Shared(:) = 0.
#endif  /*USE_MPI && defined(PARTICLES)*/

! Calculate element volumes and characteristic lengths
DO iElem = 1,nElems
  CNElemID=iElem+offsetElemCNProc
  !--- Calculate and save volume of element iElem
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    ElemVolume_Shared(CNElemID) = ElemVolume_Shared(CNElemID) + &
                                  N_Inter(Nloc)%wGP(i)*N_Inter(Nloc)%wGP(j)*N_Inter(Nloc)%wGP(k)/N_VolMesh(iElem)%sJ(i,j,k)
  END DO; END DO; END DO
  !---- Calculate characteristic cell length: V^(1/3)
  ElemCharLength_Shared(CNElemID) = ElemVolume_Shared(CNElemID)**(1./3.)
END DO

#if USE_MPI && defined(PARTICLES)
CALL BARRIER_AND_SYNC(ElemVolume_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemCharLength_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI && defined(PARTICLES)*/

! Proc-local mesh volume
LocalVolume = SUM(ElemVolume_Shared(offsetElemCNProc+1:offsetElemCNProc+nElems))

#if USE_MPI
#ifdef PARTICLES
! Compute-node mesh volume
CALL MPI_ALLREDUCE(LocalVolume,CNVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SHARED,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  ! All-reduce between node leaders
  CALL MPI_ALLREDUCE(CNVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
! Broadcast from node leaders to other processors on the same node
CALL MPI_BCAST(MeshVolume,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iERROR)
#else
! In this case, no shared array is created and all arrays are processor-local
CALL MPI_ALLREDUCE(LocalVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,IERROR)
#endif /*PARTICLES*/
#else
MeshVolume = LocalVolume
#endif /*USE_MPI*/

LBWRITE(UNIT_StdOut,'(A,E18.8)') ' |              Total MESH Volume |                ', MeshVolume
LBWRITE(UNIT_stdOut,'(A)')' INIT ELEMENT GEOMETRY INFORMATION DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitElemVolumes


SUBROUTINE setSideRanges()
!----------------------------------------------------------------------------------------------------------------------------------!
! Set the ranges in the different side lists
!
!-----------------|-----------------|-------------------|
!    U_master     | U_slave         |    FLUX           |
!-----------------|-----------------|-------------------|
!  BCsides        |                 |    BCSides        |
!  InnerMortars   |                 |    InnerMortars   |
!  InnerSides     | InnerSides      |    InnerSides     |
!  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
!                 | MPI_YOUR sides  |    MPI_YOUR sides |
!  MPIMortars     |                 |    MPIMortars     |
!-----------------|-----------------|-------------------|
!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals   ,ONLY: abort
USE MOD_Mesh_Vars ,ONLY: firstBCSide,firstMortarInnerSide,firstInnerSide,firstMPISide_MINE,firstMPISide_YOUR
USE MOD_Mesh_Vars ,ONLY: nMPISides_MINE,nMPISides_YOUR,nInnerSides,nMortarInnerSides,nBCSides
USE MOD_Mesh_Vars ,ONLY: lastBCSide,lastMortarInnerSide,lastInnerSide,lastMPISide_MINE,lastMPISide_YOUR,lastMortarMPISide
USE MOD_Mesh_Vars ,ONLY: firstMortarMPISide,nSides,nSidesMaster,nSidesSlave
#if USE_HDG
USE MOD_Globals   ,ONLY: UNIT_StdOut
USE MOD_Mesh_Vars ,ONLY: nGlobalUniqueSidesFromMesh,nGlobalUniqueSides,nMortarMPISides,nUniqueSides
#if USE_MPI
USE MOD_Globals   ,ONLY: myrank,MPI_COMM_PICLAS,iError
USE mpi
#endif /*USE_MPI*/
#endif /*USE_HDG*/
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: firstMasterSide     ! lower side ID of array U_master/gradUx_master...
INTEGER             :: lastMasterSide      ! upper side ID of array U_master/gradUx_master...
INTEGER             :: firstSlaveSide      ! lower side ID of array U_slave/gradUx_slave...
INTEGER             :: lastSlaveSide       ! upper side ID of array U_slave/gradUx_slave...
!===================================================================================================================================

firstBCSide          = 1
firstMortarInnerSide = firstBCSide         +nBCSides          ! nBCSides is set in ReadMesh()
firstInnerSide       = firstMortarInnerSide+nMortarInnerSides ! nMortarInnerSides is set in setLocalSideIDs()
firstMPISide_MINE    = firstInnerSide      +nInnerSides       ! nInnerSides is set in ReadMesh()
firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE    ! nMPISides_MINE is set in setLocalSideIDs()
firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR    ! nMPISides_YOUR is set in setLocalSideIDs()

lastBCSide           = firstMortarInnerSide-1
lastMortarInnerSide  = firstInnerSide    -1
lastInnerSide        = firstMPISide_MINE -1
lastMPISide_MINE     = firstMPISide_YOUR -1
lastMPISide_YOUR     = firstMortarMPISide-1
lastMortarMPISide    = nSides

firstMasterSide = 1
lastMasterSide  = nSides
firstSlaveSide  = firstInnerSide
lastSlaveSide   = lastMPISide_YOUR
nSidesMaster    = lastMasterSide-firstMasterSide+1
nSidesSlave     = lastSlaveSide -firstSlaveSide+1

! Set nGlobalUniqueSides: Note that big mortar sides are appended to the end of the list
#if USE_HDG
nUniqueSides = lastMPISide_MINE + nMortarMPISides !big mortars are at the end of the side list!
#if USE_MPI
CALL MPI_ALLREDUCE(nUniqueSides,nGlobalUniqueSides,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#else
nGlobalUniqueSides=nSides
#endif /*USE_MPI*/
! Sanity check: Compare the number of global unique sides with the value that is read from the mesh file
IF(nGlobalUniqueSidesFromMesh.NE.nGlobalUniqueSides)THEN!.AND.ChangedPeriodicBC) THEN ! FUTURE: use this variable when
                                                          ! nGlobalUniqueSides is calculated from mesh info
  IPWRITE(UNIT_StdOut,*) "nUniqueSides              =",nUniqueSides
  IPWRITE(UNIT_StdOut,*) "nGlobalUniqueSidesFromMesh=",nGlobalUniqueSidesFromMesh
  IPWRITE(UNIT_StdOut,*) "nGlobalUniqueSides        =",nGlobalUniqueSides
  CALL abort( &
      __STAMP__, &
      "nGlobalUniqueSides for HDG not equal the one from meshfile.")
      !"nGlobalUniqueSides for HDG not equal the one from meshfile even though no periodic sides have been changed to non-periodic.")
END IF
#endif /*HDG*/

END SUBROUTINE setSideRanges


SUBROUTINE FinalizeMesh()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
#if defined(PARTICLES) && USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
USE MOD_DG_Vars          ,ONLY: DG_Elems_master,DG_Elems_slave, N_DG_Mapping_Shared, N_DG
#if USE_MPI
USE MOD_DG_Vars          ,ONLY: N_DG_Mapping_Shared_Win
USE MOD_MPI_Shared_Vars  ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_MPI_Vars         ,ONLY: DataSizeSideSend,DataSizeSurfSendMax,DataSizeSurfSendMin
USE MOD_MPI_Vars         ,ONLY: DataSizeSideRec,DataSizeSurfRecMax,DataSizeSurfRecMin
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance!,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(ElemInfo)
! mapping from elems to sides and vice-versa
SDEALLOCATE(ElemToSide)
SDEALLOCATE(AnalyzeSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(GlobalUniqueSideID)
! elem-xgp and metrics
!#ifdef CODE_ANALYZE
!#ifndef PARTICLES
!SDEALLOCATE(SideBoundingBoxVolume)
!#endif
!#endif
SDEALLOCATE(ElemToElemGlob)
! mortars
SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)
SDEALLOCATE(MortarSlave2MasterInfo)
! mappings
SDEALLOCATE(DetJac_Ref)
SDEALLOCATE(Vdm_CLNGeo1_CLNGeo)
SDEALLOCATE(wBaryCL_NGeo1)
SDEALLOCATE(XiCL_NGeo1)
SDEALLOCATE(N_Mesh)
MeshInitIsDone = .FALSE.
SDEALLOCATE(ElemBaryNGeo)
SDEALLOCATE(ElemGlobalID)
SDEALLOCATE(myInvisibleRank)
SDEALLOCATE(LostRotPeriodicSides)
SDEALLOCATE(SideToNonUniqueGlobalSide)

! p-adaption
SDEALLOCATE(DG_Elems_master)
SDEALLOCATE(DG_Elems_slave)
SDEALLOCATE(n_surfmesh)

#if USE_MPI

SDEALLOCATE(DataSizeSideSend)
SDEALLOCATE(DataSizeSideRec)

! Min/Max is only required for HDG
SDEALLOCATE(DataSizeSurfSendMax)
SDEALLOCATE(DataSizeSurfSendMin)
SDEALLOCATE(DataSizeSurfRecMax)
SDEALLOCATE(DataSizeSurfRecMin)

! Do not deallocate during load balance here as it needs to be communicated between the processors
#if USE_LOADBALANCE
IF(.NOT.PerformLoadBalance.AND.ALLOCATED(N_DG))THEN
#endif /*USE_LOADBALANCE*/
  ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  CALL UNLOCK_AND_FREE(N_DG_Mapping_Shared_Win)
  ! Then, free the pointers or arrays
  ADEALLOCATE(N_DG_Mapping_Shared)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

SDEALLOCATE(N_DG)
SDEALLOCATE(Xi_NGeo)

! Arrays are being shifted along load balancing, so they need to be kept allocated
#if defined(PARTICLES) && USE_LOADBALANCE
IF (PerformLoadBalance .AND. .NOT.UseH5IOLoadBalance) RETURN
!IF (PerformLoadBalance) RETURN
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/

! geometry information and VDMS
!SDEALLOCATE(DCL_N)
!SDEALLOCATE(DCL_NGeo)
!SDEALLOCATE(Vdm_CLN_GaussN)

!SDEALLOCATE(Vdm_CLNGeo_CLN)
!SDEALLOCATE(Vdm_NGeo_CLNGeo)

SDEALLOCATE(XiCL_NGeo)
SDEALLOCATE(wBaryCL_NGeo)
! particle input/output
!SDEALLOCATE(Vdm_EQ_N)
!SDEALLOCATE(Vdm_N_EQ)
!SDEALLOCATE(Vdm_N_GL)
! BCS
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
! elem-xgp and metrics
SDEALLOCATE(XCL_NGeo)  ! MPI communication during LB in ExchangeMetrics()
SDEALLOCATE(dXCL_NGeo) ! MPI communication during LB in ExchangeMetrics()
! VolMesh
SDEALLOCATE(N_VolMesh)  ! MPI communication during LB in ExchangeVolMesh()
SDEALLOCATE(N_VolMesh2) ! MPI communication during LB in ExchangeMetrics()

END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
