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

INTERFACE SwapMesh
  MODULE PROCEDURE SwapMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

INTERFACE GetMeshMinMaxBoundaries
  MODULE PROCEDURE GetMeshMinMaxBoundaries
END INTERFACE

PUBLIC::InitMesh
PUBLIC::SwapMesh
PUBLIC::FinalizeMesh
PUBLIC::GetMeshMinMaxBoundaries
!===================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS

!==================================================================================================================================
!> Define parameters for Mesh
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateLogicalOption( 'DoSwapMesh'          , "Flag to swap mesh for calculation."                                                                                                 , '.FALSE.')
CALL prms%CreateStringOption(  'SwapMeshExePath'     , "(relative) path to swap-meshfile (mandatory).")
CALL prms%CreateIntOption(     'SwapMeshLevel'       , "0: initial grid\n1: first swap mesh\n2: second swap mesh\n"                                                                         , '0')
CALL prms%CreateStringOption(  'MeshFile'            , "(relative) path to meshfile (mandatory)\n(HALOWIKI:) usually located in directory of project.ini")
CALL prms%CreateLogicalOption( 'useCurveds'          , "Controls usage of high-order information in mesh. Turn off to discard high-order data and treat curved meshes as linear meshes."    , '.FALSE.')
CALL prms%CreateRealOption(    'meshScale'           , "Scale the mesh by this factor (shrink/enlarge)."                                                                                    , '1.0')
CALL prms%CreateLogicalOption( 'meshdeform'          , "Apply simple sine-shaped deformation on cartesion mesh (for testing)."                                                              , '.FALSE.')
CALL prms%CreateLogicalOption( 'meshCheckRef'        , "Flag if the mesh Jacobians should be checked in the reference system in addition to the computational system."                      , '.TRUE.')
CALL prms%CreateLogicalOption( 'CalcMeshInfo'        , 'Calculate and output elem data for myrank, ElemID and tracking info to ElemData'                                                    , '.FALSE.')
CALL prms%CreateLogicalOption( 'crossProductMetrics' , "Compute mesh metrics using cross product form. Caution: in this case free-stream preservation is only guaranteed for N=3*NGeo."     , '.FALSE.')
CALL prms%CreateStringOption(  'BoundaryName'        , "Names of boundary conditions to be set (must be present in the mesh!). For each BoundaryName a BoundaryType needs to be specified." , multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType'        , "Type of boundary conditions to be set. Format: (BC_TYPE, BC_STATE)"                                                                 , multiple=.TRUE. , no=2)
CALL prms%CreateLogicalOption( 'writePartitionInfo'  , "Write information about MPI partitions into a file."                                                                                , '.FALSE.')

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
USE MOD_Metrics                ,ONLY: BuildElem_xGP,CalcMetrics
USE MOD_Prepare_Mesh           ,ONLY: setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools            ,ONLY: PrintOption
USE MOD_ReadInTools            ,ONLY: GETLOGICAL,GETSTR,GETREAL,GETINT,GETREALARRAY
#if USE_MPI
USE MOD_Prepare_Mesh           ,ONLY: exchangeFlip
#endif
#if USE_LOADBALANCE
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
USE MOD_DG_Vars                ,ONLY: N_DG,DG_Elems_master,DG_Elems_slave
USE MOD_Particle_Mesh_Vars     ,ONLY: meshScale
! IMPLICIT VARIABLE HANDLING
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !<  0: only read and build Elem_xGP,
                               !< -1: as 0 + build connectivity and read node info (automatically read for PARTICLES=ON)
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
!CHARACTER(32)       :: hilf2
CHARACTER(LEN=255)  :: FileName
INTEGER             :: Nloc,iSide
REAL                :: RandVal
REAL                :: x1,r
LOGICAL             :: validMesh,ExistFile,ReadNodes
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

! SwapMesh: either supply the path to the swapmesh binary or place the binary into the current working directory
DoSwapMesh=GETLOGICAL('DoSwapMesh','.FALSE.')
IF(DoSwapMesh)THEN
  SwapMeshExePath=GETSTR('SwapMeshExePath','')
  INQUIRE(File=SwapMeshExePath,EXIST=ExistFile)
  IF(.NOT.ExistFile)THEN ! no path to binary found, look for binary in current directory
    FileName='./swapmesh'
    INQUIRE(File=FileName,EXIST=ExistFile)
    IF(.NOT.ExistFile) THEN
      LBWRITE(UNIT_stdOut,'(A)') ' ERROR: no swapmesh binary found'
      LBWRITE(UNIT_stdOut,'(A,A)') ' FileName:             ',TRIM(FileName)
      LBWRITE(UNIT_stdOut,'(A,L1)') ' ExistFile:            ',ExistFile
      DoSwapMesh=.FALSE.
    ELSE
      SwapMeshExePath=FileName
    END IF
  END IF
  SwapMeshLevel=GETINT('SwapMeshLevel','0')
  IF((SwapMeshLevel.LT.0).OR.(SwapMeshLevel.GT.99))THEN
    CALL abort(__STAMP__,'SwapMeshLEvel<0 or SwapMeshLEvel>99, this is invalid!')
  END IF
END IF

! prepare pointer structure (get nElems, etc.)
IF (PRESENT(MeshFile_IN)) THEN
  MeshFile = MeshFile_IN
ELSE
  MeshFile = GETSTR('MeshFile')
END IF
validMesh = ISVALIDMESHFILE(MeshFile)
IF(.NOT.validMesh) &
    CALL CollectiveStop(__STAMP__,'ERROR - Mesh file ['//TRIM(MeshFile)//'] is not a valid HDF5 mesh.')


useCurveds=GETLOGICAL('useCurveds')
#if USE_LOADBALANCE
IF ( (DoLoadBalance.OR.DoInitialAutoRestart) .AND. (.NOT.DoWriteStateToHDF5) .AND. UseH5IOLoadBalance) THEN
  DoWriteStateToHDF5=.TRUE.
  CALL PrintOption('Loadbalancing or InitialAutoRestart enabled: DoWriteStateToHDF5','INFO',LogOpt=DoWriteStateToHDF5)
END IF
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)
  LBWRITE(UNIT_stdOut,'(A67,I2.0)') ' |                           NGeo |                                ', NGeo
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
    LBWRITE(*,*) ' WARNING: Using linear elements although NGeo>1!'
  END IF
END IF

meshScale = GETREAL('meshScale') ! default is 1.0
! Sanity check
IF(ABS(meshScale).LE.0.) CALL abort(__STAMP__,'meshScale is zero')
CALL ReadMesh(MeshFile,ReadNodes) !set nElems

!schmutz fink
PP_nElems=nElems

Coords=>NodeCoords
nElemsLoc=nElems

! scale and deform mesh if desired (warning: no mesh output!)
IF(ABS(meshScale-1.).GT.1e-14) Coords = Coords*meshScale

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  DO iElem=1,nElems
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PI*x(1))*SIN(PI*x(2))*SIN(PI*x(3))
    END DO; END DO; END DO;
  END DO
END IF

! allocate arrays and initialize local polynomial degree
ALLOCATE(N_DG(nElems))
! By default, the initial degree is set to PP_N
!N_DG = 1
!N_DG(1) = PP_N
N_DG = PP_N
!N_DG(1) = 5

!             DO iElem=1,nElems
!               !CALL RANDOM_NUMBER(RandVal)
!               !N_DG (iElem) = 1+INT(RandVal*Nmax)
!               !N_DG (iElem) = MAX(N_DG (iElem),4)
!                 kLoop: DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
!                   x1 = coords(1,i,j,k,iElem)
!                   r = SQRT(coords(1,i,j,k,iElem)**2+&
!                            coords(2,i,j,k,iElem)**2+&
!                            coords(3,i,j,k,iElem)**2  )
!             
!               IF(r.lt.1.0 .and. x1.lt.0)THEN
!                 N_DG(iElem) = PP_N
!                 EXIT kLoop
!               ELSE
!                 N_DG(iElem) = PP_N-6-1
!                 !EXIT kLoop
!               END IF ! r.le.1.0 .and. x.le.0
!                 END DO ; END DO; END DO kLoop;
!             END DO

! Build Elem_xGP
ALLOCATE(N_VolMesh(1:nElems))
CALL BuildElem_xGP(NodeCoords)

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
  CALL fillMeshInfo(ABS(meshMode))
#else
  CALL fillMeshInfo()
#endif /*USE_HDG && USE_LOADBALANCE*/

  ! build necessary mappings
  ALLOCATE(N_Mesh(Nmin:Nmax))
  DO Nloc = Nmin, Nmax
    CALL InitMappings(Nloc, N_Mesh(Nloc)%VolToSideA, N_Mesh(Nloc)%VolToSideIJKA, N_Mesh(Nloc)%FS2M)
  END DO ! Nloc = Nmin, Nmax

END IF ! meshMode.GT.0

CALL InitpAdaption()

IF (ABS(meshMode).GT.1) THEN

  ! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----

  NGeoRef=3*NGeo ! build jacobian at higher degree
  ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

  ! volume data
  DO iElem = 1, nElems
    Nloc = N_DG(iElem)
    ALLOCATE(N_VolMesh(iElem)%dXCL_N(      3,3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%Metrics_fTilde(3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%Metrics_gTilde(3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%Metrics_hTilde(3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%sJ            (  0:Nloc,0:Nloc,0:Nloc))
  END DO ! iElem = 1, nElems

  ! surface data
  ALLOCATE(N_SurfMesh(1:nSides))
  DO iSide = 1, nSides

    ! Allocate with max. polynomial degree of the two master-slave sides
    Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
    ALLOCATE(N_SurfMesh(iSide)%Face_xGP (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%NormVec  (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%TangVec1 (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%TangVec2 (3,0:Nloc,0:Nloc))
    ALLOCATE(N_SurfMesh(iSide)%SurfElem (  0:Nloc,0:Nloc))
    N_SurfMesh(iSide)%Face_xGP = 0.
    N_SurfMesh(iSide)%NormVec  = 0.
    N_SurfMesh(iSide)%TangVec1 = 0.
    N_SurfMesh(iSide)%TangVec2 = 0.
    N_SurfMesh(iSide)%SurfElem = 0.
  END DO ! iSide = 1, nSides

#ifdef maxwell
#if defined(ROS) || defined(IMPA)
  ALLOCATE(nVecLoc(1:3,0:PP_N,0:PP_N,1:6,PP_nElems))
  ALLOCATE(SurfLoc(0:PP_N,0:PP_N,1:6,PP_nElems))
  nVecLoc=0.
  SurfLoc=0.
#endif /*ROS or IMPA*/
#endif /*maxwell*/


! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
! assign 1/detJ (sJ)
! assign normal and tangential vectors and surfElems on faces

  crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
  LBWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."

  ! get XCL_NGeo
  ALLOCATE(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
  XCL_NGeo = 0.
#if defined(PARTICLES)
  ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
  dXCL_NGeo = 0.
#endif /*defined(PARTICLES)*/

#ifdef PARTICLES
  CALL CalcMetrics(XCL_NGeo_Out=XCL_NGeo, dXCL_NGeo_Out=dXCL_NGeo)
#else
  CALL CalcMetrics(XCL_NGeo_Out=XCL_NGeo)
#endif

#if defined(PARTICLES) || USE_HDG
  ! Compute element bary and element radius for processor-local elements (without halo region)
  ALLOCATE(ElemBaryNGeo(1:3,1:nElems))
  CALL BuildElementOrigin()
#endif /*defined(PARTICLES) || USE_HDG*/

#ifndef PARTICLES
  ! dealloacte pointers
  LBWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
  CALL deleteMeshPointer()
#endif

  ! Initialize element volumes and characteristic lengths
  CALL InitElemVolumes()

#ifndef PARTICLES
  DEALLOCATE(NodeCoords)
#endif

  IF(ABS(meshMode).NE.3)THEN
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

MeshInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


!===================================================================================================================================
!> Get local N for each element and side
!===================================================================================================================================
SUBROUTINE InitpAdaption()
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY: DG_Elems_master,DG_Elems_slave,N_DG
USE MOD_IO_HDF5   ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars ,ONLY: nSides,nElems
!USE MOD_DG        ,ONLY: DG_ProlongDGElemsToFace
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Read p-adaption specific input data
!pAdaption    = GETLOGICAL('pAdaption','.FALSE.')

! add array containing the local polynomial degree to the hdf5 output
CALL AddToElemData(ElementOut,'Nloc',IntArray=N_DG)
ALLOCATE(DG_Elems_master(1:nSides))
ALLOCATE(DG_Elems_slave (1:nSides))

! Set polynomial degree at the element sides
DG_Elems_master = PP_N
DG_Elems_slave  = PP_N

CALL DG_ProlongDGElemsToFace()

END SUBROUTINE InitpAdaption


!==================================================================================================================================
!> Set DG_Elems_slave and DG_Elems_master information
!==================================================================================================================================
SUBROUTINE DG_ProlongDGElemsToFace()
! MODULES
USE MOD_GLobals
USE MOD_DG_Vars   ,ONLY: N_DG,DG_Elems_master,DG_Elems_slave
USE MOD_Mesh_Vars ,ONLY: SideToElem,nSides,nBCSides
#if USE_MPI
USE MOD_MPI       ,ONLY: StartExchange_DG_Elems,FinishExchangeMPIData
USE MOD_MPI_Vars  ,ONLY: DataSizeSideSend,DataSizeSideRec,nNbProcs,nMPISides_rec,nMPISides_send,OffsetMPISides_rec
USE MOD_MPI_Vars  ,ONLY: OffsetMPISides_send,DGExchange
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iSide,ElemID,nbElemID
#if USE_MPI
INTEGER                      :: iNbProc,Nloc
INTEGER, DIMENSION(nNbProcs) :: RecRequest_U,SendRequest_U,RecRequest_U2,SendRequest_U2
#endif /*USE_MPI*/
!==================================================================================================================================
! set information which polynomial degree elements adjacent to a side have
DO iSide = 1,nSides
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,iSide)
  ! Master sides
  IF(ElemID  .GT.0) DG_Elems_master(iSide) = N_DG(ElemID)
  ! Slave side (ElemID,locSide and flip =-1 if not existing)
  IF(nbElemID.GT.0) DG_Elems_slave( iSide) = N_DG(nbElemID)
  ! Boundaries
  IF(iSide.LE.nBCSides) DG_Elems_slave( iSide) = DG_Elems_master(iSide)
END DO

#if USE_MPI
! Exchange element local polynomial degree (N_LOC)
CALL StartExchange_DG_Elems(DG_Elems_slave ,1,nSides,SendRequest_U ,RecRequest_U ,SendID=2)  ! RECEIVE MINE, SEND YOUR / DG_Elems_slave:  slave  -> master                                                                                       ! Send MINE, receive YOUR / DG_Elems_slave : slave  -> master
CALL StartExchange_DG_Elems(DG_Elems_master,1,nSides,SendRequest_U2,RecRequest_U2,SendID=1)  ! RECEIVE YOUR, SEND MINE / DG_Elems_master: master -> slave    
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U ,RecRequest_U ,SendID=2) ! Send YOUR - receive MINE
CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=1) ! Send YOUR - receive MINE
! Initialize the send/rec face sizes for master/slave communication
IF(.NOT.ALLOCATED(DataSizeSideSend)) ALLOCATE(DataSizeSideSend(nNbProcs,2), DataSizeSideRec(nNbProcs,2))
DataSizeSideSend=0
DataSizeSideRec=0
DO iNbProc = 1, nNbProcs

  ! 1: Set number of sides and offset for SEND MINE - RECEIVE YOUR case
  IF(nMPISides_rec(iNbProc,1).GT.0) THEN
    DO iSide = OffsetMPISides_rec(iNbProc-1,1)+1, OffsetMPISides_rec(iNbProc,1)
      Nloc = DG_Elems_slave(iSide) ! polynomial degree of the sending slave side
      DataSizeSideRec(iNbProc,1) = DataSizeSideRec(iNbProc,1)  + (Nloc+1)**2
    END DO ! iSide = 1, nMPISides_rec(iNbProc,1)
  END IF

  IF(nMPISides_send(iNbProc,1).GT.0) THEN
    DO iSide = OffsetMPISides_send(iNbProc-1,1)+1, OffsetMPISides_send(iNbProc,1)
      Nloc = DG_Elems_slave(iSide) ! polynomial degree of the sending slave side
      DataSizeSideSend(iNbProc,1) = DataSizeSideSend(iNbProc,1)  + (Nloc+1)**2
    END DO ! iSide = 1, nMPISides_rec(iNbProc,1)
  END IF

  ! Important: Keep symmetry
  DataSizeSideSend(iNbProc,2) = DataSizeSideRec( iNbProc,1)
  DataSizeSideRec( iNbProc,2) = DataSizeSideSend(iNbProc,1)
END DO ! iNbProc = 1, nNbProcs

! Create unrolled arrays
ALLOCATE(DGExchange(nNbProcs))
DO iNbProc=1,nNbProcs
  ALLOCATE(DGExchange(iNbProc)%FaceDataRecv(PP_nVar, MAXVAL(DataSizeSideRec(iNbProc,:))))
  ALLOCATE(DGExchange(iNbProc)%FaceDataSend(PP_nVar, MAXVAL(DataSizeSideSend(iNbProc,:))))
  IPWRITE(*,*)'send1',DataSizeSideSend(:,1)
  IPWRITE(*,*)'send2',DataSizeSideSend(:,2)

  IPWRITE(*,*)'recv1',DataSizeSideRec(:,1)
  IPWRITE(*,*)'recv2',DataSizeSideRec(:,2)
END DO !iProc=1,nNBProcs
#endif /*USE_MPI*/
END SUBROUTINE DG_ProlongDGElemsToFace


SUBROUTINE SwapMesh()
!============================================================================================================================
! use the posti tool swapmesh in order to map the DG solution as well as particles into a new state file with a different
! mesh file
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars     ,ONLY: ProjectName
USE MOD_Preproc
USE MOD_Mesh_Vars        ,ONLY: SwapMeshExePath,SwapMeshLevel,MeshFile
USE MOD_Restart_Vars     ,ONLY: RestartFile
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
CHARACTER(LEN=255)  :: LocalName                   ! local project name
CHARACTER(LEN=255)  :: FileName                    ! names/paths to files
CHARACTER(LEN=255)  :: NewFolderName,OldFolderName ! names/paths to folders
LOGICAL             :: objExist                    ! file or folder exist variable
LOGICAL             :: CreateSwapScript            ! temporary varaible for creating a swapmesh.sh shell script
LOGICAL             :: LogSwapMesh                 ! create log file for swapmesh: log.swapmesh
LOGICAL             :: CleanUp                     ! rm old state and mesh file in swapmesh folder
LOGICAL             :: KeepSwapFile                ! if true, do not remove created swap file, e.g., "PlasmaPlume_NewMesh_State...."
CHARACTER(LEN=22)   :: ParameterFile               ! swapmesh parameter file containing all conversion information
CHARACTER(LEN=3)    :: hilf,hilf2                  ! auxiliary variable for INTEGER -> CHARACTER conversion
CHARACTER(LEN=255)  :: SYSCOMMAND,SWITCHFOLDER     ! string to fit the system command
INTEGER             :: iSTATUS                     ! error status return code
INTEGER             :: StartIndex                  ! get string index (position) of certain substring
CHARACTER(len=255)  :: cwd                         ! current working directory
INTEGER             :: cwdLen                      ! string length of cwd
!============================================================================================================================
ParameterFile='parameter_swapmesh.ini'
IF(SwapMeshLevel.GT.0)THEN
  LocalName=TRIM(ProjectName(1:LEN(TRIM(ProjectName))-3))
ELSE
  LocalName=TRIM(ProjectName)
END IF

! current mesh folder
WRITE(UNIT=hilf,FMT='(I3)') SwapMeshLevel
IF(SwapMeshLevel.GT.9)THEN
  OldFolderName='mesh'//TRIM(ADJUSTL(hilf))
ELSE
  OldFolderName='mesh0'//TRIM(ADJUSTL(hilf))
ENDIF

CALL getcwd(cwd)
cwdLen=LEN(TRIM(cwd))
!WRITE(*,*) TRIM(cwd)
!print*,TRIM(cwd(cwdLen-5:cwdLen)) ! e.g. : PlasmaPlume_cylinder/linear/mesh02 -> mesh02
!print*,TRIM(OldFolderName)        ! e.g. : mesh02
IF(TRIM(cwd(cwdLen-LEN(TRIM(OldFolderName))+1:cwdLen)).NE.TRIM(OldFolderName))THEN
  SWRITE(UNIT_stdOut,'(A)')   ' ERROR: something is wrong with the directory or SwapMeshLevel'
  SWRITE(UNIT_stdOut,'(A,A)') ' cwd:           ',TRIM(cwd)
  SWRITE(UNIT_stdOut,'(A,I3,A,A)') ' SwapMeshLevel: ',SwapMeshLevel,' -> ',OldFolderName
  CALL abort(__STAMP__,'Abort. Check SwapMeshLevel.',999,999.)
END IF

! check if next level mesh folder exists
WRITE(UNIT=hilf,FMT='(I3)') SwapMeshLevel+1
IF(SwapMeshLevel+1.GT.9)THEN
  NewFolderName='mesh'//TRIM(ADJUSTL(hilf))
ELSE
  NewFolderName='mesh0'//TRIM(ADJUSTL(hilf))
ENDIF
INQUIRE(File='../'//TRIM(NewFolderName),EXIST=objExist)
IF(.NOT.objExist)THEN ! no swapmesh folder found
  LBWRITE(UNIT_stdOut,'(A)')   ' ERROR: no swapmesh folder found. Cannot perform swapmesh'
  LBWRITE(UNIT_stdOut,'(A,A)') ' objName:             ',TRIM(NewFolderName)
  LBWRITE(UNIT_stdOut,'(A,L1)') ' objExist:            ',objExist
ELSE
  FileName='../'//TRIM(NewFolderName)//'/'//ParameterFile
  INQUIRE(File=FileName,EXIST=objExist)
  IF(.NOT.objExist) THEN
    LBWRITE(UNIT_stdOut,'(A)')   ' ERROR: no '//ParameterFile//' found. Cannot perform swapmesh'
    LBWRITE(UNIT_stdOut,'(A,A)') ' objName:             ',TRIM(FileName)
    LBWRITE(UNIT_stdOut,'(A,L1)') ' objExist:            ',objExist
  ELSE
    !NewFolderName=FileName

    SWITCHFOLDER='cd ../'//TRIM(NewFolderName)//' && '
    ! print LocalName and SwapMeshLevel as new mesh file name to parameter_swapmesh.ini
    ! cd ../meshXX &&  sed -i -e "s/.*MeshFileNew.*/MeshFileNew=PlasmaPlume_01_mesh.h5/" parameter_swapmesh.ini
    SYSCOMMAND=' sed -i -e "s/.*MeshFileNew.*/MeshFileNew='//&
               TRIM(LocalName)//'_'//TRIM(NewFolderName(5:6))//'_mesh.h5/" '//ParameterFile
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)

    ! print current polynomial degree to parameter_swapmesh.ini
    WRITE(UNIT=hilf,FMT='(I3)') PP_N
    ! cd ../meshXX && sed -i -e "s/.*NNew.*/NNew=2/" parameter_swapmesh.ini
    SYSCOMMAND=' sed -i -e "s/.*NNew.*/NNew='//&
               TRIM(ADJUSTL(hilf))//'/" '//ParameterFile
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)

    ! create symbolic link to old mesh file and restart file
    ! (delete the old links if they exist, they might point to the wrong file)
    ! cd ../meshXX && rm PlasmaPlume_State_000.000000000000100.h5 > /dev/null 2>&1
    SYSCOMMAND=' rm '//TRIM(RestartFile)//' > /dev/null 2>&1'
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    ! cd ../meshXX && ln -s ../mesh00/PlasmaPlume_State_000.000000000000100.h5 .
    SYSCOMMAND=' ln -s ../'//TRIM(OldFolderName)//'/'//TRIM(RestartFile)//' .'
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)

    ! cd ../meshXX && rm PlasmaPlume_mesh.h5 > /dev/null 2>&1
    SYSCOMMAND=' rm '//TRIM(MeshFile)//' > /dev/null 2>&1'
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    ! cd ../meshXX && ln -s ../mesh00/PlasmaPlume_mesh.h5 .
    SYSCOMMAND=' ln -s ../'//TRIM(OldFolderName)//'/'//TRIM(MeshFile)//' .'
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)

    ! call swapmesh and create new state file (delete all '_NewMesh_State*' state file of new mesh)
    ! cd ../meshXX && rm PlasmaPlume_NewMesh_State* > /dev/null 2>&1
    SYSCOMMAND=' rm '//TRIM(ProjectName)//'_NewMesh_State* > /dev/null 2>&1'
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    SYSCOMMAND=' echo "#!/bin/bash" > swapmesh.sh'
    print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
    CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)

    ! run swapmesh executable
    ! =======================================================================================================================
    ! CAUTION: CURRENTLY ONLY WORKS IN SINGLE EXECUTION!!!! WHEN STARTED WITH MPIRUN AND PERFORMED BY MPIROOT NOTHING OCCURS!
    ! =======================================================================================================================
    CreateSwapScript=.FALSE.
    IF(CreateSwapScript)THEN
      ! cd ../mesh01 && echo "/home/stephen/Flexi/ParaViewPlugin_newest_version/build_hdf16/bin/swapmesh parameter_swapmesh.ini
      !                       PlasmaPlume_State_000.000000000000100.h5" >> swapmesh.sh && chmod +x swapmesh.sh
      SYSCOMMAND=' '//ParameterFile//' '//TRIM(RestartFile)//'" >> swapmesh.sh && chmod +x swapmesh.sh'
      print*,                   TRIM(SWITCHFOLDER)//' echo "'//TRIM(SwapMeshExePath)//TRIM(SYSCOMMAND)
      CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//' echo "'//TRIM(SwapMeshExePath)//TRIM(SYSCOMMAND), &
      WAIT=.TRUE., EXITSTAT=iSTATUS)
      !SYSCOMMAND=' nohup ./swapmesh.sh > log.swapmesh &'
      SYSCOMMAND=' ./swapmesh.sh > log.swapmesh '
      print*,TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
      CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    ELSE
      ! cd ../meshXX && ~/build_hdf16/bin/swapmesh parameter_swapmesh.ini PlasmaPlume_State_000.000000000000100.h5
      LogSwapMesh=.TRUE.
      IF(LogSwapMesh)THEN
        SYSCOMMAND=' '//ParameterFile//' '//TRIM(RestartFile)//' > log.swapmesh'
      ELSE
        SYSCOMMAND=' '//ParameterFile//' '//TRIM(RestartFile)
      END IF
      print*,                   TRIM(SWITCHFOLDER)//' '//TRIM(SwapMeshExePath)//TRIM(SYSCOMMAND)
      CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//' '//TRIM(SwapMeshExePath)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    END IF



    ! cd ../meshXX && mv PlasmaPlume_NewMesh_State_000.000000000000100.h5 PlasmaPlume_State_000.000000000000100.h5
    WRITE(UNIT=hilf,FMT='(I3)') SwapMeshLevel+1
    IF(SwapMeshLevel+1.GT.9)THEN
      hilf2=TRIM(ADJUSTL(hilf))
    ELSE
      hilf2='0'//TRIM(ADJUSTL(hilf))
    ENDIF


    KeepSwapFile=.FALSE.
    IF(.NOT.KeepSwapFile)THEN
      ! move the created state file (remove the "NewMesh" suffix and replace with SwapMeshLevel)
      StartIndex=INDEX(TRIM(RestartFile),'_State')+7
      ! cd ../mesh01 && mv PlasmaPlume_01_NewMesh_State_000.000000000002000.h5 PlasmaPlume_01_State_000.000000000002000.h5
      SYSCOMMAND=' mv '//TRIM(ProjectName)//'_NewMesh_State_'//TRIM(RestartFile(StartIndex:LEN(RestartFile)))//' '//&
                        TRIM(LocalName)//'_'//TRIM(ADJUSTL(hilf2))//'_State_'//TRIM(RestartFile(StartIndex:LEN(RestartFile)))
      print*,TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
      CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    END IF
    IF(iSTATUS.NE.0)THEN
      CALL abort(&
      __STAMP__&
      ,'new swapmesh state file could not be created.',999,999.)
    END IF


    CleanUp=.TRUE.
    IF(CleanUp)THEN
      ! Clean up: rm old state and mesh file !
      ! cd ../meshXX && rm PlasmaPlume_State_000.000000000000100.h5 > /dev/null 2>&1
      SYSCOMMAND=' rm '//TRIM(RestartFile)//' > /dev/null 2>&1'
      print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
      CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
      ! cd ../mesh02 && rm PlasmaPlume_01_mesh.h5 > /dev/null 2>&1
      SYSCOMMAND=' rm '//TRIM(MeshFile)//' > /dev/null 2>&1'
      print*,                   TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND)
      CALL EXECUTE_COMMAND_LINE(TRIM(SWITCHFOLDER)//TRIM(SYSCOMMAND), WAIT=.TRUE., EXITSTAT=iSTATUS)
    END IF

!PlasmaPlume_State_






  END IF
END IF

END SUBROUTINE SwapMesh


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
   CALL MPI_ALLREDUCE(xmin_loc(1:3),xmin(1:3), 3, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERROR)
   ! Find global max
   CALL MPI_ALLREDUCE(xmax_loc(1:3),xmax(1:3), 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)

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
USE MOD_Globals            ,ONLY: UNIT_StdOut
USE MOD_Interpolation_Vars ,ONLY: N_Inter
USE MOD_DG_Vars            ,ONLY: N_DG
USE MOD_Mesh_Vars          ,ONLY: nElems,N_VolMesh
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
  Nloc = N_DG(iElem)
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
CNVolume = SUM(ElemVolume_Shared(:))
IF (myComputeNodeRank.EQ.0) THEN
  ! All-reduce between node leaders
  CALL MPI_ALLREDUCE(CNVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
! Broadcast from node leaders to other processors on the same node
CALL MPI_BCAST(MeshVolume,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iERROR)
#else
! In this case, no shared array is created and all arrays are processor-local
CALL MPI_ALLREDUCE(LocalVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
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
USE MOD_Globals   ,ONLY: myrank
USE MOD_Globals   ,ONLY: iError,MPI_COMM_WORLD
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
CALL MPI_ALLREDUCE(nUniqueSides,nGlobalUniqueSides,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
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
USE MOD_LoadBalance_Vars     ,ONLY: PerformLoadBalance
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
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
! geometry information and VDMS
SDEALLOCATE(Xi_NGeo)
SDEALLOCATE(DCL_N)
SDEALLOCATE(DCL_NGeo)
SDEALLOCATE(VdM_CLN_GaussN)
SDEALLOCATE(Vdm_CLNGeo_CLN)
SDEALLOCATE(Vdm_NGeo_CLNgeo)
SDEALLOCATE(Vdm_EQ_N)
SDEALLOCATE(Vdm_N_EQ)
SDEALLOCATE(Vdm_GL_N)
SDEALLOCATE(Vdm_N_GL)
! mapping from elems to sides and vice-versa
SDEALLOCATE(ElemToSide)
SDEALLOCATE(AnalyzeSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(GlobalUniqueSideID)
! elem-xgp and metrics
#ifdef maxwell
#if defined(ROS) || defined(IMPA)
SDEALLOCATE(nVecLoc)
SDEALLOCATE(SurfLoc)
#endif /*ROS or IMPA*/
#endif /*maxwell*/
!#ifdef CODE_ANALYZE
!#ifndef PARTICLES
!SDEALLOCATE(SideBoundingBoxVolume)
!#endif
!#endif
SDEALLOCATE(ElemToElemGlob)
SDEALLOCATE(XCL_NGeo)
!SDEALLOCATE(dXCL_NGeo)
SDEALLOCATE(wbaryCL_NGeo)
SDEALLOCATE(XiCL_NGeo)
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

#if defined(PARTICLES) && USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
! BCS
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)

END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
