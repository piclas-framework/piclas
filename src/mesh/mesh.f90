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
PUBLIC::InitMesh
PUBLIC::FinalizeMesh
PUBLIC::GetMeshMinMaxBoundaries
PUBLIC::DefineParametersMesh
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for Mesh
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools    ,ONLY: prms,addStrListEntry
USE MOD_Mesh_pAdaption ,ONLY: PRM_P_ADAPTION_ZERO,PRM_P_ADAPTION_RDN,PRM_P_ADAPTION_NPB,PRM_P_ADAPTION_HH
USE MOD_Mesh_pAdaption ,ONLY: PRM_P_ADAPTION_LVL_MINTWO,PRM_P_ADAPTION_LVL_MINONE,PRM_P_ADAPTION_LVL_DEFAULT,PRM_P_ADAPTION_LVL_TWO
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
#if USE_FV
CALL prms%CreateLogicalOption( 'meshCheckRef-FV'        , "Flag if the mesh Jacobians should be checked in the reference system in addition to the computational system."                      , '.TRUE.')
CALL prms%CreateIntArrayOption('BoundaryType-FV'     , "Type of boundary conditions for FV to be set. Format: (BC_TYPE, BC_STATE)"                                                                 , multiple=.TRUE. , no=2)
#endif
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
#if USE_FV
USE MOD_Mesh_Vars_FV
USE MOD_Metrics_FV             ,ONLY: CalcMetrics_PP_1,CalcSurfMetrics_PP_1
#endif /*FV*/
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
#if USE_FV
USE MOD_LoadBalance_Metrics_FV ,ONLY: ExchangeVolMesh_FV,ExchangeMetrics_FV
#endif /*USE_FV*/
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
USE MOD_DSMC_Vars              ,ONLY: DoRadialWeighting, DoLinearWeighting, DoCellLocalWeighting
USE MOD_Particle_Vars          ,ONLY: usevMPF
#endif
#if USE_HDG && USE_LOADBALANCE
USE MOD_Mesh_Tools             ,ONLY: BuildSideToNonUniqueGlobalSide
#endif /*USE_HDG && USE_LOADBALANCE*/
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars                ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_Particle_Mesh_Vars     ,ONLY: meshScale
USE MOD_Mesh_Vars              ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_pAdaption         ,ONLY: InitpAdaption
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
INTEGER             :: Nloc,iSide,NSideMin,N_max
LOGICAL             :: validMesh,ReadNodes
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

! initialize flag for gradient calculations (filled in setLocalSideIDs)
#if USE_FV
ALLOCATE(IsPeriodicSide(1:nSides))
IsPeriodicSide(:)=.FALSE.
#endif /*USE_FV*/

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
#if USE_FV
  CALL ExchangeVolMesh_FV()
#endif
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

#if !(PP_TimeDiscMethod==700)
CALL DG_ProlongDGElemsToFace() ! Builds DG_Elems_master and ,DG_Elems_slave, requires SideToElem()
#endif /*!(PP_TimeDiscMethod==700)*/

! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----

IF (ABS(meshMode).GT.1) THEN
#if USE_LOADBALANCE
  IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  !IF (PerformLoadBalance) THEN
    ! Shift metric arrays during load balance
    CALL ExchangeMetrics()
#if USE_FV
    CALL ExchangeMetrics_FV()
#endif
  ELSE
#endif /*USE_LOADBALANCE*/

    NGeoRef=3*NGeo ! build jacobian at higher degree
    ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

    ! volume data
    DO iElem = 1, nElems
#if !(PP_TimeDiscMethod==700)
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
#else
      Nloc = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/
      ALLOCATE(N_VolMesh(iElem)%XCL_N(       3,  0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%Metrics_fTilde(3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%Metrics_gTilde(3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%Metrics_hTilde(3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh(iElem)%sJ            (  0:Nloc,0:Nloc,0:Nloc))
    END DO ! iElem = 1, nElems

    ! volume data
    DO iElem = 1, nElems
#if !(PP_TimeDiscMethod==700)
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
#else
      Nloc = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/
      ALLOCATE(N_VolMesh2(iElem)%dXCL_N(     3,3,0:Nloc,0:Nloc,0:Nloc))
      ALLOCATE(N_VolMesh2(iElem)%JaCL_N(     3,3,0:Nloc,0:Nloc,0:Nloc))
    END DO ! iElem = 1, nElems

#if USE_FV
    SDEALLOCATE(       XCL_N_PP_1)
    SDEALLOCATE(      dXCL_N_PP_1)
    SDEALLOCATE(      JaCL_N_PP_1)

    ALLOCATE(       XCL_N_PP_1(  3,0:PP_1,0:PP_1,0:PP_1,nElems)) ! built in CalcMetrics(), required in CalcSurfMetrics()
    ALLOCATE(      dXCL_N_PP_1(3,3,0:PP_1,0:PP_1,0:PP_1,nElems)) ! built in CalcMetrics()
    ALLOCATE(      JaCL_N_PP_1(3,3,0:PP_1,0:PP_1,0:PP_1,nElems))
    JaCL_N_PP_1 = 0.

    ! Vandermonde
    ALLOCATE(Vdm_CLN_N_PP_1(   0:PP_1,0:PP_1))
#endif /*FV*/
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

#if USE_FV
#if defined(PARTICLES)
    ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
    dXCL_NGeo = 0.
    CALL CalcMetrics_PP_1(XCL_NGeo_Out=XCL_NGeo,dXCL_NGeo_Out=dXCL_NGeo)
#else
    ! CALL InitGetCNElemID() ! TODO: is this required here?
    CALL CalcMetrics_PP_1(XCL_NGeo_Out=XCL_NGeo)
#endif /*PARTICLES*/
#endif /*USE_FV*/
#if !(USE_FV) || (USE_HDG)
#ifdef PARTICLES
    IF(.NOT.ALLOCATED(dXCL_NGeo)) ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
    dXCL_NGeo = 0.
    CALL CalcMetrics(XCL_NGeo_Out=XCL_NGeo,dXCL_NGeo_Out=dXCL_NGeo)
#else
    CALL CalcMetrics(XCL_NGeo_Out=XCL_NGeo)
#endif /*PARTICLES*/
#endif /*!(USE_FV) || (USE_HDG)*/
#if USE_LOADBALANCE
  END IF ! PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)
#endif /*USE_LOADBALANCE*/

  ! surface data
  ALLOCATE(N_SurfMesh(1:nSides))
  DO iSide = 1, nSides

    ! Allocate with max. polynomial degree of the two master-slave sides
#if !(PP_TimeDiscMethod==700)
    N_max    = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
    NSideMin = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
#else
    N_max    = PP_N
    NSideMin = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/
    N_SurfMesh(iSide)%NSide = NSideMin
    ALLOCATE(N_SurfMesh(iSide)%Face_xGP (3,0:N_max,0:N_max))
    ALLOCATE(N_SurfMesh(iSide)%NormVec  (3,0:N_max,0:N_max))
    ALLOCATE(N_SurfMesh(iSide)%TangVec1 (3,0:N_max,0:N_max))
    ALLOCATE(N_SurfMesh(iSide)%TangVec2 (3,0:N_max,0:N_max))
    ALLOCATE(N_SurfMesh(iSide)%SurfElem (  0:N_max,0:N_max))
    !!!ALLOCATE(N_SurfMesh(iSide)%SurfElemMin(0:NSideMin,0:NSideMin))
    N_SurfMesh(iSide)%Face_xGP = 0.
    N_SurfMesh(iSide)%NormVec  = 0.
    N_SurfMesh(iSide)%TangVec1 = 0.
    N_SurfMesh(iSide)%TangVec2 = 0.
    N_SurfMesh(iSide)%SurfElem = 0.
    !!!N_SurfMesh(iSide)%SurfElemMin= 0.
  END DO ! iSide = 1, nSides

#if !(USE_FV) || (USE_HDG)
  ! Due to possible load balance, this is done outside of CalcMetrics() now
  DO iElem=1,nElems
    CALL CalcSurfMetrics(iElem)
  END DO ! iElem = 1, nElems
#endif

#if USE_FV
  ALLOCATE(Face_xGP_FV      (3,1:nSides))
  ALLOCATE(NormVec_FV       (3,1:nSides))
  ALLOCATE(TangVec1_FV      (3,1:nSides))
  ALLOCATE(TangVec2_FV      (3,1:nSides))
  ALLOCATE(SurfElem_FV      (  1:nSides))
  ALLOCATE(Ja_Face_FV       (3,3,1:nSides)) ! temp
  Face_xGP_FV=0.
  NormVec_FV=0.
  TangVec1_FV=0.
  TangVec2_FV=0.
  SurfElem_FV=0.
  Ja_Face_FV=0.

  ALLOCATE(Face_xGP_PP_1   (3,0:PP_1,0:PP_1,1:nSides))
  ALLOCATE(NormVec_PP_1    (3,0:PP_1,0:PP_1,1:nSides))
  ALLOCATE(TangVec1_PP_1   (3,0:PP_1,0:PP_1,1:nSides))
  ALLOCATE(TangVec2_PP_1   (3,0:PP_1,0:PP_1,1:nSides))
  ALLOCATE(SurfElem_PP_1     (0:PP_1,0:PP_1,1:nSides))
  ALLOCATE(Ja_Face_PP_1  (3,3,0:PP_1,0:PP_1,1:nSides))
  Face_xGP_PP_1       = 0.
  NormVec_PP_1        = 0.
  TangVec1_PP_1       = 0.
  TangVec2_PP_1       = 0.
  SurfElem_PP_1       = 0.
  Ja_Face_PP_1        = 0.

  DO iElem=1,nElems
    CALL CalcSurfMetrics_PP_1(JaCL_N_PP_1(:,:,:,:,:,iElem),iElem)
  END DO

  ! PP_1 metrics to global ones
  Face_xGP_FV(:,:)   = SUM(SUM(Face_xGP_PP_1      (:,:,:,:),3),2)/4. !average over (PP_1+1)^2 points
  NormVec_FV (:,:)   = NormVec_PP_1       (:,0,0,:)
  TangVec1_FV(:,:)   = TangVec1_PP_1      (:,0,0,:)
  TangVec2_FV(:,:)   = TangVec2_PP_1      (:,0,0,:)
  SurfElem_FV(:)     = SUM(SUM(SurfElem_PP_1(:,:,:),2),1)
  Ja_Face_FV (:,:,:) = SUM(SUM(Ja_Face_PP_1(:,:,:,:,:),4),3)
#endif /*FV*/

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
#endif /*PARTICLES*/
#if USE_FV
  DEALLOCATE(Ja_Face_PP_1)
  DEALLOCATE(Ja_Face_FV)
#endif /*FV*/

  IF((ABS(meshMode).NE.3).AND.(meshMode.GT.1))THEN
#ifdef PARTICLES
    IF(DoRadialWeighting.OR.DoLinearWeighting.OR.DoCellLocalWeighting) THEN
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


#if !(PP_TimeDiscMethod==700)
!==================================================================================================================================
!> Set DG_Elems_slave and DG_Elems_master information
!==================================================================================================================================
SUBROUTINE DG_ProlongDGElemsToFace()
! MODULES
USE MOD_PreProc
USE MOD_GLobals
USE MOD_DG_Vars   ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave,N_DG_Mapping!,pAdaptionType
USE MOD_Mesh_Vars ,ONLY: SideToElem,nSides,nBCSides, offSetElem

USE MOD_Mesh_Vars ,ONLY: firstMortarInnerSide,lastMortarInnerSide,MortarType,MortarInfo
#if USE_MPI
USE MOD_Mesh_Vars ,ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_MPI       ,ONLY: StartExchange_DG_Elems,FinishExchangeMPIData
USE MOD_MPI_Vars  ,ONLY: DataSizeSideSend,DataSizeSideRec,nNbProcs,nMPISides_rec,nMPISides_send,OffsetMPISides_rec
USE MOD_MPI_Vars  ,ONLY: OffsetMPISides_send
USE MOD_MPI_Vars  ,ONLY: DataSizeSurfSendMax,DataSizeSurfRecMax, DataSizeSurfSendMin,DataSizeSurfRecMin
#if !(USE_HDG)
USE MOD_MPI_Vars  ,ONLY: DataSizeSideSendMaster,DataSizeSideRecMaster
#endif /*not USE_HDG*/
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
TYPE(MPI_Request), DIMENSION(nNbProcs) :: RecRequest_U,SendRequest_U,RecRequest_U2,SendRequest_U2
#endif /*USE_MPI*/
!==================================================================================================================================
! Side containers
ALLOCATE(DG_Elems_master(1:nSides))
ALLOCATE(DG_Elems_slave (1:nSides))
! Initialize with element-local N
!IF(pAdaptionType.EQ.0)THEN
  DG_Elems_master = -1
  DG_Elems_slave  = -1
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
ALLOCATE(DataSizeSideSend(nNbProcs,2), DataSizeSideRec(nNbProcs,2))
! Min/Max is only required for HDG
ALLOCATE(DataSizeSurfSendMax(nNbProcs,2), DataSizeSurfRecMax(nNbProcs,2))
ALLOCATE(DataSizeSurfSendMin(nNbProcs,2), DataSizeSurfRecMin(nNbProcs,2))
#if !(USE_HDG)
! Master is only required for Maxwell (with Dielectric)
ALLOCATE(DataSizeSideSendMaster(nNbProcs,2), DataSizeSideRecMaster(nNbProcs,2))
DataSizeSideSendMaster=0
DataSizeSideRecMaster=0
#endif /*not USE_HDG*/
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
#if !(USE_HDG)
      ! Master is only required for Maxwell (with Dielectric)
      Nloc = DG_Elems_master(iSide) ! polynomial degree of the sending slave side
      DataSizeSideRecMaster(iNbProc,1) = DataSizeSideRecMaster(iNbProc,1)  + (Nloc+1)**2
#endif /*not USE_HDG*/
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
#if !(USE_HDG)
      ! Master is only required for Maxwell (with Dielectric)
      Nloc = DG_Elems_master(iSide) ! polynomial degree of the sending slave side
      DataSizeSideSendMaster(iNbProc,1) = DataSizeSideSendMaster(iNbProc,1)  + (Nloc+1)**2
#endif /*not USE_HDG*/
    END DO ! iSide = 1, nMPISides_rec(iNbProc,1)
END IF

  ! Important: Keep symmetry
  DataSizeSideSend(iNbProc,2) = DataSizeSideRec( iNbProc,1)
  DataSizeSideRec( iNbProc,2) = DataSizeSideSend(iNbProc,1)
  DataSizeSurfSendMax(iNbProc,2) = DataSizeSurfRecMax( iNbProc,1)
  DataSizeSurfRecMax( iNbProc,2) = DataSizeSurfSendMax(iNbProc,1)
  DataSizeSurfSendMin(iNbProc,2) = DataSizeSurfRecMin( iNbProc,1)
  DataSizeSurfRecMin( iNbProc,2) = DataSizeSurfSendMin(iNbProc,1)
#if !(USE_HDG)
  ! Master is only required for Maxwell (with Dielectric)
  DataSizeSideSendMaster(iNbProc,2) = DataSizeSideRecMaster( iNbProc,1)
  DataSizeSideRecMaster( iNbProc,2) = DataSizeSideSendMaster(iNbProc,1)
#endif /*not USE_HDG*/
END DO ! iNbProc = 1, nNbProcs

#endif /*USE_MPI*/
END SUBROUTINE DG_ProlongDGElemsToFace
#endif /*!(PP_TimeDiscMethod==700)*/


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
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_Mesh_Vars          ,ONLY: nElems,N_VolMesh, offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: LocalVolume,MeshVolume
USE MOD_Particle_Mesh_Vars ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
USE MOD_ReadInTools
#if USE_MPI
USE mpi_f08
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
#if (USE_FV) && !(USE_HDG)
! Can something be checked in this case?
#else
! Sanity check
IF(.NOT.ALLOCATED(N_Inter)) CALL abort(__STAMP__,'Error in InitElemVolumes(): N_Inter is not allocated')
#endif /*(USE_FV) && !(USE_HDG)*/

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
#if (USE_FV) && !(USE_HDG)
    ElemVolume_Shared(CNElemID) = 1./N_VolMesh(iElem)%sJ(0,0,0)
#else
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    ElemVolume_Shared(CNElemID) = ElemVolume_Shared(CNElemID) + &
                                  N_Inter(Nloc)%wGP(i)*N_Inter(Nloc)%wGP(j)*N_Inter(Nloc)%wGP(k)/N_VolMesh(iElem)%sJ(i,j,k)
  END DO; END DO; END DO
#endif /*(USE_FV) && !(USE_HDG)*/
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
USE mpi_f08
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
#if USE_FV
USE MOD_Mesh_Vars_FV
#endif /*USE_FV*/
#if defined(PARTICLES) && USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars          ,ONLY: DG_Elems_master,DG_Elems_slave, N_DG
#endif /*!(PP_TimeDiscMethod==700)*/
#if USE_MPI
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars          ,ONLY: N_DG_Mapping_Shared, N_DG_Mapping_Shared_Win
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_MPI_Shared_Vars  ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_MPI_Vars         ,ONLY: DataSizeSideSend,DataSizeSurfSendMax,DataSizeSurfSendMin
USE MOD_MPI_Vars         ,ONLY: DataSizeSideRec,DataSizeSurfRecMax,DataSizeSurfRecMin
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance!,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
#if !(USE_HDG)
USE MOD_MPI_Vars         ,ONLY: DataSizeSideSendMaster,DataSizeSideRecMaster
#endif /*not USE_HDG*/
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
! mortars
SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)
SDEALLOCATE(MortarSlave2MasterInfo)
! mappings
SDEALLOCATE(DetJac_Ref)
#if USE_FV
SDEALLOCATE(NormVec_FV)
SDEALLOCATE(TangVec1_FV)
SDEALLOCATE(TangVec2_FV)
SDEALLOCATE(SurfElem_FV)
SDEALLOCATE(Face_xGP_FV)
SDEALLOCATE(NormVec_PP_1)
SDEALLOCATE(TangVec1_PP_1)
SDEALLOCATE(TangVec2_PP_1)
SDEALLOCATE(SurfElem_PP_1)
SDEALLOCATE(Face_xGP_PP_1)
SDEALLOCATE(IsPeriodicSide)
#endif
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

SDEALLOCATE(n_surfmesh)

! p-adaption
#if !(PP_TimeDiscMethod==700)
SDEALLOCATE(DG_Elems_master)
SDEALLOCATE(DG_Elems_slave)

#if USE_MPI
SDEALLOCATE(DataSizeSideSend)
SDEALLOCATE(DataSizeSideRec)

! Min/Max is only required for HDG
SDEALLOCATE(DataSizeSurfSendMax)
SDEALLOCATE(DataSizeSurfSendMin)
SDEALLOCATE(DataSizeSurfRecMax)
SDEALLOCATE(DataSizeSurfRecMin)

#if !(USE_HDG)
! Master is only required for Maxwell (with Dielectric)
SDEALLOCATE(DataSizeSideSendMaster)
SDEALLOCATE(DataSizeSideRecMaster)
#endif /*not USE_HDG*/

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
#endif /*!(PP_TimeDiscMethod==700)*/
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

#if USE_FV
SDEALLOCATE(BoundaryType_FV)
SDEALLOCATE(Elem_xGP_FV)
SDEALLOCATE(Elem_xGP_PP_1)
SDEALLOCATE(JaCL_N_PP_1)
SDEALLOCATE(       XCL_N_PP_1)
SDEALLOCATE(      dXCL_N_PP_1)
#endif
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh