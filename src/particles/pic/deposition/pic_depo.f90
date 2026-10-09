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

MODULE MOD_PICDepo
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
IMPLICIT NONE
PRIVATE

TYPE NodeDepoMapping
  INTEGER                                     :: NodeID
  TYPE (NodeDepoMapping), POINTER             :: next => NULL()
END TYPE
!===================================================================================================================================
PUBLIC:: Deposition, InitializeDeposition, FinalizeDeposition, DefineParametersPICDeposition
PUBLIC:: InitDepoSurfNodes
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for PIC Deposition
!==================================================================================================================================
SUBROUTINE DefineParametersPICDeposition()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools      ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%CreateStringOption( 'PIC-TimeAverageFile'      , 'Read charge density from .h5 file and save to PartSource\n'//&
                                                           'WARNING: Currently not correctly implemented for shared memory', 'none')
CALL prms%CreateLogicalOption('PIC-RelaxDeposition'      , 'Relaxation of current PartSource with RelaxFac\n'//&
                                                           'into PartSourceOld', '.FALSE.')
CALL prms%CreateRealOption(   'PIC-RelaxFac'             , 'Relaxation factor of current PartSource with RelaxFac\n'//&
                                                           'into PartSourceOld', '0.001')

CALL prms%CreateRealOption(   'PIC-shapefunction-radius'             , 'Radius of shape function'   , '1.')
CALL prms%CreateIntOption(    'PIC-shapefunction-alpha'              , 'Exponent of shape function' , '2')
CALL prms%CreateIntOption(    'PIC-shapefunction-dimension'          , '1D, 2D or 3D shape function', '3')
CALL prms%CreateIntOption(    'PIC-shapefunction-direction'          , &
    'Only required for PIC-shapefunction-dimension 1 or 2: Shape function direction for 1D (the direction in which the charge '//&
    'will be distributed) and 2D (the direction in which the charge will be constant)', '1')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-3D-deposition' ,'Deposit the charge over volume (3D)\n'//&
                                                                   ' or over a line (1D)/area(2D)\n'//&
                                                                   '1D shape function: volume or line\n'//&
                                                                   '2D shape function: volume or area', '.TRUE.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'Minimum shape function radius (for cylindrical and spherical)', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'Scaling factor of shape function radius '//&
                                                             '(for cylindrical and spherical)', '0.')
CALL prms%CreateRealOption(     'PIC-shapefunction-adaptive-DOF'  ,'Average number of DOF in shape function radius (assuming a '//&
    'Cartesian grid with equal elements). Only implemented for PIC-Deposition-Type = shape_function_adaptive (2). The maximum '//&
    'number of DOF is limited by the polynomial degree and depends on PIC-shapefunction-dimension=1, 2 or 3\n'//&
   '1D: 2*(N+1)\n'//&
   '2D: Pi*(N+1)^2\n'//&
   '3D: (4/3)*Pi*(N+1)^3\n')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-adaptive-smoothing', 'Enable smooth transition of element-dependent radius when'//&
                                                                      ' using shape_function_adaptive.', '.FALSE.')

END SUBROUTINE DefineParametersPICDeposition


!===================================================================================================================================
!> Initialize the deposition variables first
!===================================================================================================================================
SUBROUTINE InitializeDeposition()
! MODULES
USE MOD_Globals
USE MOD_Basis                  ,ONLY: BarycentricWeights,InitializeVandermonde
USE MOD_Basis                  ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_Interpolation_Vars     ,ONLY: N_Inter
USE MOD_Mesh_Vars              ,ONLY: nElems,N_VolMesh,offSetElem
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueGlobalNodes, GEO
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Tools          ,ONLY: CalcCellLocNodeVolumes,ReadTimeAverage
USE MOD_PICInterpolation_Vars  ,ONLY: InterpolationType
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID, GetCNElemID
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Symmetry_Vars          ,ONLY: Symmetry
#if USE_MPI
USE MOD_PICDepo_MPI            ,ONLY: InitDepoNodesMPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nProcessors_Global
! USE MOD_MPI_Shared
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Interpolation_Vars     ,ONLY: Nmin,Nmax
USE MOD_DG_Vars                ,ONLY: N_DG_Mapping
USE MOD_Particle_Boundary_Vars ,ONLY: Do2DSurfaceCharge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ALLOCSTAT, iElem, i, j, k, kk, ll, mm, Nloc
CHARACTER(255)            :: TimeAverageFile
#if USE_MPI
LOGICAL,ALLOCATABLE       :: DoNodeMapping(:), SendNode(:)
#else
INTEGER                   :: iNode
#endif
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE DEPOSITION...'

IF(.NOT.DoDeposition) THEN
  ! fill deposition type with empty string
  DepositionType='NONE'
  OutputSource=.FALSE.
  RelaxDeposition=.FALSE.
  RETURN
END IF

#if USE_LOADBALANCE
! Not "LB via MPI" means during 1st initialisation
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  ALLOCATE(PS_N(nElems))
  !--- Allocate arrays for charge density collection and initialize
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ALLOCATE(PS_N(iElem)%PartSource(1:4,0:Nloc,0:Nloc,0:Nloc))
    PS_N(iElem)%PartSource = 0.0
  END DO ! iElem = 1, nElems
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

!--- check if relaxation of current PartSource with RelaxFac into PartSourceOld
RelaxDeposition = GETLOGICAL('PIC-RelaxDeposition','F')
IF (RelaxDeposition) THEN
  RelaxFac     = GETREAL('PIC-RelaxFac','0.001')
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
#if ((USE_HDG) && (PP_nVar==1))
    ALLOCATE(PS_N(iElem)%PartSourceOld(1  ,1:2,0:Nloc,0:Nloc,0:Nloc),STAT=ALLOCSTAT)
#else
    ALLOCATE(PS_N(iElem)%PartSourceOld(1:4,1:2,0:Nloc,0:Nloc,0:Nloc),STAT=ALLOCSTAT)
#endif
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in pic_depo.f90: Cannot allocate PartSourceOld!')
    PS_N(iElem)%PartSourceOld = 0.
  END DO ! iElem = 1, nElems
  OutputSource = .TRUE.
ELSE
  OutputSource = GETLOGICAL('PIC-OutputSource','F')
END IF

!--- check if charge density is computed from TimeAverageFile
TimeAverageFile = GETSTR('PIC-TimeAverageFile','none')
IF (TRIM(TimeAverageFile).NE.'none') THEN
  CALL abort(__STAMP__,'This feature is currently not working! PartSource must be correctly handled in shared memory context.')
  CALL ReadTimeAverage(TimeAverageFile)
  IF (.NOT.RelaxDeposition) THEN
  !-- switch off deposition: use only the read PartSource
    DoDeposition=.FALSE.
    DepositionType='constant'
    RETURN
  ELSE
  !-- use read PartSource as initialValue for relaxation
  !-- CAUTION: will be overwritten by DG_Source if present in restart-file!
    DO iElem = 1, nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      DO kk = 0, Nloc
        DO ll = 0, Nloc
          DO mm = 0, Nloc
#if ((USE_HDG) && (PP_nVar==1))
            PS_N(iElem)%PartSourceOld(1  ,1,mm,ll,kk) = PS_N(iElem)%PartSource(4  ,mm,ll,kk)
            PS_N(iElem)%PartSourceOld(1  ,2,mm,ll,kk) = PS_N(iElem)%PartSource(4  ,mm,ll,kk)
#else
            PS_N(iElem)%PartSourceOld(1:4,1,mm,ll,kk) = PS_N(iElem)%PartSource(1:4,mm,ll,kk)
            PS_N(iElem)%PartSourceOld(1:4,2,mm,ll,kk) = PS_N(iElem)%PartSource(1:4,mm,ll,kk)
#endif
          END DO !mm
        END DO !ll
      END DO !kk
    END DO !iElem
  END IF
END IF

!--- init DepositionType-specific vars
SELECT CASE(TRIM(DepositionType))
! ------------------------------------------------
CASE('cell_volweight')
! ------------------------------------------------
  ALLOCATE(CellVolWeight(Nmin:Nmax))
  ALLOCATE(CellVolWeight_Volumes(0:1,0:1,0:1,nElems))

  DO Nloc=Nmin,Nmax
    ALLOCATE(CellVolWeight(Nloc)%Fac(0:Nloc))
    CellVolWeight(Nloc)%Fac(0:Nloc) = N_Inter(Nloc)%xGP(0:Nloc)
    CellVolWeight(Nloc)%Fac(0:Nloc) = (CellVolWeight(Nloc)%Fac(0:Nloc)+1.0)/2.0
  END DO

  CellVolWeight_Volumes=0.0
  DO iElem=1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    DO i=0,Nloc;DO j=0,Nloc;DO k=0,Nloc
      ASSOCIATE( sJ => N_VolMesh(iElem)%sJ , xGP => N_Inter(Nloc)%xGP , wGP => N_Inter(Nloc)%wGP , &
               CVWV => CellVolWeight_Volumes(:,:,:,iElem) )
        ! CVWV cannot be accessed here with "0" because of the associate construct!
        CVWV(1,1,1) = CVWV(1,1,1) + 1./sJ(i,j,k)*((1.-xGP(i))*(1.-xGP(j))*(1.-xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(1,1,2) = CVWV(1,1,2) + 1./sJ(i,j,k)*((1.-xGP(i))*(1.-xGP(j))*(1.+xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(1,2,1) = CVWV(1,2,1) + 1./sJ(i,j,k)*((1.-xGP(i))*(1.+xGP(j))*(1.-xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(1,2,2) = CVWV(1,2,2) + 1./sJ(i,j,k)*((1.-xGP(i))*(1.+xGP(j))*(1.+xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(2,1,1) = CVWV(2,1,1) + 1./sJ(i,j,k)*((1.+xGP(i))*(1.-xGP(j))*(1.-xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(2,1,2) = CVWV(2,1,2) + 1./sJ(i,j,k)*((1.+xGP(i))*(1.-xGP(j))*(1.+xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(2,2,1) = CVWV(2,2,1) + 1./sJ(i,j,k)*((1.+xGP(i))*(1.+xGP(j))*(1.-xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
        CVWV(2,2,2) = CVWV(2,2,2) + 1./sJ(i,j,k)*((1.+xGP(i))*(1.+xGP(j))*(1.+xGP(k))*wGP(i)*wGP(j)*wGP(k)/8.)
      END ASSOCIATE
    END DO; END DO; END DO
  END DO
! ------------------------------------------------
CASE('cell_volweight_mean')
! ------------------------------------------------
#if USE_MPI
  ALLOCATE(DoNodeMapping(0:nProcessors_Global-1),SendNode(1:nUniqueGlobalNodes))
  DoNodeMapping = .FALSE.
  SendNode = .FALSE.
#endif
  IF ((TRIM(InterpolationType).NE.'cell_volweight')) THEN
    ALLOCATE(CellVolWeight(Nmin:Nmax))
    DO Nloc=Nmin,Nmax
      ALLOCATE(CellVolWeight(Nloc)%Fac(0:Nloc))
      CellVolWeight(Nloc)%Fac(0:Nloc) = N_Inter(Nloc)%xGP(0:Nloc)
      CellVolWeight(Nloc)%Fac(0:Nloc) = (CellVolWeight(Nloc)%Fac(0:Nloc)+1.0)/2.0
    END DO
  END IF

  ALLOCATE(NodeSource(1:4,1:nUniqueGlobalNodes))
  NodeSource=0.0
  IF(DoDielectricSurfaceCharge)THEN
    ALLOCATE(NodeSourceExt(1:nUniqueGlobalNodes))
    NodeSourceExt    = 0.
  END IF ! DoDielectricSurfaceCharge

  IF (GEO%nPeriodicVectors.GT.0) CALL InitializePeriodicNodes(&
#if USE_MPI
                                        DoNodeMapping,SendNode&
#endif /*USE_MPI*/
      )

#if USE_MPI
  CALL InitDepoNodesMPI(DoNodeMapping,SendNode)
#else
  nDepoNodes      = nUniqueGlobalNodes
  nDepoNodesTotal = nDepoNodes
  ALLOCATE(DepoNodetoGlobalNode(1:nDepoNodesTotal))
  DO iNode=1, nUniqueGlobalNodes
    DepoNodetoGlobalNode(iNode) = iNode
  END DO
#endif /*USE_MPI*/


! Initialize sub-cell volumes around nodes
  CALL CalcCellLocNodeVolumes()

! ------------------------------------------------
CASE('shape_function', 'shape_function_cc', 'shape_function_adaptive')
! ------------------------------------------------
  !--- Allocate arrays for shape function charge density collection and initialize
  ALLOCATE(N_ShapeTmp(NMax))
  DO Nloc = 1, NMax
    ALLOCATE(N_ShapeTmp(Nloc)%PartSource(1:4,0:Nloc,0:Nloc,0:Nloc))
    N_ShapeTmp(Nloc)%PartSource = 0.0
  END DO ! iElem = 1, nElems

#if USE_MPI
  ALLOCATE(RecvRequest(nShapeExchangeProcs),SendRequest(nShapeExchangeProcs))
#endif
  ! --- Set shape function radius in each cell when using adaptive shape function
  IF(TRIM(DepositionType).EQ.'shape_function_adaptive') CALL InitShapeFunctionAdaptive()

  ! --- Set integration factor only for uncorrected shape function methods
  SELECT CASE(TRIM(DepositionType))
  CASE('shape_function_cc', 'shape_function_adaptive')
    w_sf  = 1.0 ! set dummy value
  CASE('shape_function')
    IF(Symmetry%axisymmetric) CALL abort(__STAMP__,'Axisymmetric simulations only with shape_function_cc or shape_function_adaptive!')
  END SELECT

  ! --- Set periodic case matrix for shape function deposition (virtual displacement of particles in the periodic directions)
  CALL InitAxisymmetrySF()
  CALL InitPeriodicSFCaseMatrix()

  ! --- Set element flag for cycling already completed elements
#if USE_MPI
  ALLOCATE(ChargeSFDone(1:nComputeNodeTotalElems))
#else
  ALLOCATE(ChargeSFDone(1:nElems))
#endif /*USE_MPI*/
! ------------------------------------------------
CASE('cell_mean')
! ------------------------------------------------
CASE DEFAULT
! ------------------------------------------------
  CALL abort(__STAMP__,'Unknown DepositionType in pic_depo.f90')
END SELECT

! Surface charge model
IF (Do2DSurfaceCharge) THEN
  CALL InitDepoSurfNodes() ! Get nDepoSurfNodes
END IF ! DoSurfaceCharge

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition


!===================================================================================================================================
!> Find all surface nodes connected to a BC where surface deposition is active
!>
!> 1. Loop over the processor-local elements
!> 2. Loop over the corner vertices of the element
!> 3. Use VertexConnectInfo to get the neighbour element index and vertex for all possible connection (also periodic)
!> 4. Check the sides of connected to the neighbour node and find out if the side is a BC side
!===================================================================================================================================
SUBROUTINE InitDepoSurfNodes()
! MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars ,ONLY: nNonUniqueGlobalNodes
USE MOD_Mesh_Vars          ,ONLY: readFEMconnectivity,nGlobalElems
USE MOD_Particle_Mesh_Vars ,ONLY: VertexConnectInfo_shared
USE MOD_Mesh_Vars          ,ONLY: NGeo,NonUniqueGlobalSideIDToNonUniqueGlobalNodeID!,SideToNonUniqueGlobalSide
USE MOD_Mesh_Vars          ,ONLY: BoundaryType,nFEMVertices,NonUniqueGlobalNodeIDToFEMVertexID,nBCSides,BC,nSides
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,SideInfo_Shared,ElemInfo_Shared,VertexInfo_Shared
USE MOD_Mesh_pAdaption     ,ONLY: getlocsidelist
USE MOD_Mesh_Tools         ,ONLY: GetCornerNodeMapCGNS,GetCNElemID,GetGlobalElemID
USE MOD_Particle_Mesh_Vars ,ONLY: nNonUniqueGlobalSides,ElemSideNodeID_Shared
USE MOD_Interpolation_Vars ,ONLY: Nmax
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVISU,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
#if USE_MPI
USE MOD_PICDepo_MPI        ,ONLY: InitDepoSurfNodesMPI,CollectSurfNodeAreaOnMPIRoot,ReverseExchangeSurfNodeArea
USE MOD_MPI_Shared_Vars    ,ONLY: myComputeNodeRank,nComputeNodeTotalElems,nComputeNodeProcessors
USE MOD_MPI_Shared_vars    ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared         ,ONLY: BARRIER_AND_SYNC
USE MOD_Particle_Mesh_Vars ,ONLY: VertexInfo_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
USE MOD_PICDepo_MPI        ,ONLY: LBReverseExchangeSurfNodeSource
#endif /*USE_LOADBALANCE*/
USE MOD_Restart_Vars       ,ONLY: DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
LOGICAL             :: InitializeSurfNodeArrays
#else
LOGICAL,PARAMETER   :: InitializeSurfNodeArrays=.TRUE.
#endif /*USE_LOADBALANCE*/
INTEGER :: BCType,NonUniqueGlobalSideID,NonUniqueGlobalNbSideID,iGlobalElemID,BCIndex,ElemType,SideID
INTEGER :: iVertexConnect,GlobalNbElemID,NbLocVertexID,LocSideList(3),iNeighbourLocSideList,iNeighbourLocSide
INTEGER :: iLocSideList,iLocSide
INTEGER :: FirstVertexInd,LastVertexInd,FirstVertexConnectInd,LastVertexConnectInd
INTEGER :: FEMVertexID,iVertexInd,NonUniqueNodeID,CNS(8),iNode
INTEGER :: CNElemID,LocSideID
INTEGER :: FirstCNElemID,LastCNElemID,iCNELemID,localVertexID
REAL    :: StartT,EndT
INTEGER, PARAMETER :: MaxAllowedSymmetries=6 ! The number 6 is chosen at random to limit the maximum number of symmetries
REAL    :: SubSideAreaEquiN1(0:1,0:1)
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A,I0,A)',ADVANCE='NO') ' | Initializing node mappings for 2D surface deposition...'
GETTIME(StartT)
! OPTIMIZE: Can the mappings that are created here be stored in .h5 for restart purposes and when running piclas2vtk to save time?
! Sanity check: This routine requires FEM connectivity
IF(.NOT.readFEMconnectivity) CALL CollectiveStop(__STAMP__,'Error in surface deposition init: readFEMconnectivity=T is required')

#if USE_LOADBALANCE
! Set flag when performinggg load balancing: The MPIRoot keeps all arrays and does not deallocate them as it has the global mappings
InitializeSurfNodeArrays = .FALSE.
IF (.NOT.PerformLoadBalance.OR.(PerformLoadBalance.AND.(.NOT.MPIRoot))) InitializeSurfNodeArrays = .TRUE.
#endif /*USE_LOADBALANCE*/

! Surface mapping from p,q-system to iNode (node coord system)
ALLOCATE(pq2iNode(0:1,0:1,1:nSides))
pq2iNode = 0


IF (InitializeSurfNodeArrays) THEN
  ! Flag the unique deposition nodes per processor
  nDepoSurfNodes = 0
  ALLOCATE(IsDepoSurfNode(1:nFEMVertices))
  IsDepoSurfNode = .FALSE.

  ! Mapping from NonuniqueGlobalNodeID to FEMVertexID
  ! OPTIMIZE: Make this array SHM
  ALLOCATE(NonUniqueGlobalNodeIDToFEMVertexID(1:nNonUniqueGlobalNodes))
  NonUniqueGlobalNodeIDToFEMVertexID = 0

  ! For counting the number of visualisation sides
  ! OPTIMIZE: Make these arrays SHM
  nDepoSurfSides = 0
  ALLOCATE(IsDepoSurfSide(1:nNonUniqueGlobalSides))
  IsDepoSurfSide = .FALSE.
  ! 1-4: NodeIDs
  ALLOCATE(NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(1:4,1:nNonUniqueGlobalSides))
  NonUniqueGlobalSideIDToNonUniqueGlobalNodeID = 0

  ! The cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
  CALL GetCornerNodeMapCGNS(NGeo,CornerNodesCGNS = CNS)
END IF ! InitializeSurfNodeArrays

! Consider all compute-node elements where deposition takes place (this leads to sending to processes in the halo region)
! and which are required for the field solver source terms (element-local vertices, which can lead to receiving from other
! processes)
#if USE_MPI
! With SHM array: Divide the loop into sections to split the workload
FirstCNElemID = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
LastCNElemID  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
! Without SHM array: Every process loops over all elements on the CN node
FirstCNElemID = 1
LastCNElemID  = nComputeNodeTotalElems
#else
FirstCNElemID = 1
LastCNElemID  = nElems
#endif

! The MPIRoot misuses the following two variables as it needs all FEMVertexIDs for outptut to .h5 because all processes send their
! FEMVertexIDs
IF (MPIRoot) THEN
  IF (InitializeSurfNodeArrays) THEN
    FirstCNElemID = 1
    LastCNElemID  = nGlobalElems
  ELSE
    FirstCNElemID = 1
    LastCNElemID  = -1
  END IF ! InitializeSurfNodeArrays
END IF ! MPIRoot

! 1. Identify all (FEMVertexID) nodes and (NonUniqueGlobalSideID or NonUniqueGlobalNbSideID) sides that are needed for deposition
! Loop over the process-local global elements indices
DO iCNELemID = FirstCNElemID, LastCNElemID
  ! Get global element index
  IF (MPIRoot) THEN
    ! Little hack: switch CN and global variable name in loop
    iGlobalElemID = iCNELemID
  ELSE
    iGlobalElemID = GetGlobalElemID(iCNELemID)
  END IF ! MPIRoot
  ! iElem = iGlobalElemID - offsetElem
  ElemType = ElemInfo_Shared(ELEM_TYPE,iGlobalElemID)
  ! Sanity check: currently only hexahedral elements are implemented
  SELECT CASE(ElemType)
  CASE(108,118,208)
    ! Hexahedral elements
  CASE DEFAULT
    CALL abort(__STAMP__,'InitDepoSurfNodes(): Element type not implemented, ElemType =',IntInfoOpt=ElemType)
  END SELECT
  ! Get local FEMElemInfo of current element
  FirstVertexInd = ElemInfo_Shared(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1 ! this comes from FEMElemInfo() from mesh.h5
  LastVertexInd  = ElemInfo_Shared(ELEM_LASTVERTEXIND,iGlobalElemID)    ! this comes from FEMElemInfo() from mesh.h5
  ! Sanity check
  IF (FirstVertexInd.GE.LastVertexInd) THEN
    IPWRITE(*,*) 'iGlobalElemID :', iGlobalElemID
    IPWRITE(*,*) 'FirstVertexInd:', FirstVertexInd
    IPWRITE(*,*) 'LastVertexInd :', LastVertexInd
    CALL abort(__STAMP__,' FirstVertexInd >= LastVertexInd')
  END IF ! FirstVertexInd.GE.LastVertexInd
  ! Loop over all non-unique vertices (the total number via iGlobalElemID and iVertexInd corresponds to nVertices in .h5)
  DO iVertexInd = FirstVertexInd,LastVertexInd
    ! Get topologically unique global vertex ID (via VertexInfo from mesh.h5), includes periodicity (needed for a FEM solver
    FEMVertexID = VertexInfo_Shared(VERTEX_FEMID,iVertexInd)
    ! Get the local vertex index
    localVertexID = iVertexInd-FirstVertexInd+1
    ! Get the non-unique node index
    NonUniqueNodeID = CNS(localVertexID) + FirstVertexInd - 1
    ! Store mapping iVertexInd to NonUniqueNodeID
    VertexInfo_Shared(VERTEX_NONUNIQUENODEID,iVertexInd) = NonUniqueNodeID
    ! Mapping from NonUniqueNodeID to FEMVertexID
    NonUniqueGlobalNodeIDToFEMVertexID(NonUniqueNodeID) = FEMVertexID
    ! Get local vertex connectivity: First and Last connected vertex index
    FirstVertexConnectInd = VertexInfo_Shared(VERTEX_FIRSTCONNECTIND,iVertexInd)+1
    LastVertexConnectInd  = VertexInfo_Shared(VERTEX_LASTCONNECTIND,iVertexInd)
    ! Check nodes without connections
    IF (FirstVertexConnectInd.GT.LastVertexConnectInd) THEN ! Vertex has no neighbours (solo vertex)
      ! Check if any of the three connected sides is a deposition side
      ! Set sides depending on the element type: Only implemented for Hexahedral elements
      CALL GetLocSideList(ElemType,localVertexID,LocSideList)
      ! Loop over the three connected sides of the element
      iSide: DO iLocSideList = 1, 3
        ! Check if current element has already been flagged
        iLocSide = LocSideList(iLocSideList)
        ! Get non-unique global side index of the element that is connected to the FEMVertexID/NonUniqueNodeID
        NonUniqueGlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobalElemID) + iLocSide
        ! Get boundary condition index
        BCIndex = SideInfo_Shared(SIDE_BCID,NonUniqueGlobalSideID)
        IF(BCIndex.LE.0) CYCLE iSide ! Skip inner sides
        ! Get boundary condition type
        BCType = BoundaryType(BCIndex,BC_TYPE)
        ! FEATURE: Implement inner BCs for surface charge deposition
        ! IF(BCType.EQ.100) CALL abort(__STAMP__,'InitDepoSurfNodes(): Inner BCs not implemented for surface charge deposition')
        ! FEATURE: Define a list of all BCType numbers that allow surface deposition
        IF(BCType.NE.30) CYCLE iSide ! Skip non-DCBC sides
        ! Depo node/side found
        IsDepoSurfNode(FEMVertexID) = .TRUE.
        IsDepoSurfSide(NonUniqueGlobalSideID) = .TRUE.
      END DO iSide ! iLocSideList = 1, 3
    ELSE ! Vertex has neibouring vertices
      DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
        ! Get neighbour infos. Note the ABS() for +/- master/slave notation
        GlobalNbElemID = ABS(VertexConnectInfo_Shared(VERTEXCONNECT_NBELEMID   ,iVertexConnect))
        NbLocVertexID  =     VertexConnectInfo_Shared(VERTEXCONNECT_NBLOCNODEID,iVertexConnect)
        ! Set sides depending on the element type: Only implemented for Hexahedral elements
        CALL GetLocSideList(ElemType,NbLocVertexID,LocSideList)
        ! Loop over the three connected sides of the neighbour element, which is connected with a corner to iVertexConnect
        iNbSide: DO iNeighbourLocSideList = 1, 3
          ! Check if current element has already been flagged
          iNeighbourLocSide = LocSideList(iNeighbourLocSideList)
          ! Get non-unique global side index of the neighbouring element that is connected to the FEMVertexID/NonUniqueNodeID
          NonUniqueGlobalNbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalNbElemID) + iNeighbourLocSide
          ! Get boundary condition index
          BCIndex = SideInfo_Shared(SIDE_BCID,NonUniqueGlobalNbSideID)
          IF(BCIndex.LE.0) CYCLE iNbSide ! Skip inner sides
          ! Get boundary condition type
          BCType = BoundaryType(BCIndex,BC_TYPE)
          ! FEATURE: Implement inner BCs for surface charge deposition
          ! IF(BCType.EQ.100) CALL abort(__STAMP__,'InitDepoSurfNodes(): Inner BCs not implemented for surface charge deposition')
          ! FEATURE: define a list of all BCType numbers that allow surface deposition
          IF(BCType.NE.30) CYCLE iNbSide ! Skip non-DCBC sides
          ! Depo node/side found
          IsDepoSurfNode(FEMVertexID) = .TRUE.
          IsDepoSurfSide(NonUniqueGlobalNbSideID) = .TRUE.
        END DO iNbSide ! iNeighbourLocSideList = 1, 3
      END DO ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
    END IF ! FirstVertexConnectInd.GT.LastVertexConnectInd
  END DO ! iVertexInd = iFirstVertexInd,LastVertexInd
END DO !  iCNELemID = FirstCNElemID, LastCNElemID
#if USE_MPI
! NOTE: VERTEX_NONUNIQUENODEID is not filled for elements that are not on the compute node. Only the CN with MPIRoot has
! all info, because the MPIRoot processes loops over all global elements
CALL BARRIER_AND_SYNC(VertexInfo_Shared_Win,MPI_COMM_SHARED)

#endif /*USE_MPI*/

IF (InitializeSurfNodeArrays) THEN
  ! 2. Create a mapping that returns the four (NonUniqueNodeID) nodes for a (NonUniqueGlobalSideID) side
  ! Separate loop for setting the node IDs is needed because setting them in the loop above does not work
  ! Additionally, the scaling factor is determined
  DO NonUniqueGlobalSideID = 1,nNonUniqueGlobalSides
    ! Check if the side has charge deposition activated
    IF(IsDepoSurfSide(NonUniqueGlobalSideID))THEN
      ! Get global element index
      iGlobalElemID = SideInfo_Shared(SIDE_ELEMID,NonUniqueGlobalSideID)
      ! Get compute node element index of the side
      CNElemID  = GetCNElemID(iGlobalElemID)
      ! Skip global elements, which are not on the compute node and not inside the halo region by checking if the CN element index is
      IF(CNElemID.EQ.-1) CYCLE
      ! Get the local side index (1-6)
      LocSideID = SideInfo_Shared(SIDE_LOCALID,NonUniqueGlobalSideID)
      ! Loop over all 4 node of the side
      DO iNode = 1, 4
        ! Get the non-unique global side index of the node/local side ID/compute element ID
        NonUniqueNodeID = ElemSideNodeID_Shared(iNode,LocSideID,CNElemID) + 1
        ! Store the non-unique node index for the current non-unique global side index
        NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(iNode,NonUniqueGlobalSideID) = NonUniqueNodeID
      END DO ! iNode = 1, 4
    END IF
  END DO ! NonUniqueGlobalSideID = 1,nNonUniqueGlobalSides
  ! END DO

  ! Sanity check: Loop over all deposition surface side IDs and make sure the mapping is correct
  DO NonUniqueGlobalSideID = 1,nNonUniqueGlobalSides
    ! Only check surfaces that are marked for deposition
    IF(IsDepoSurfSide(NonUniqueGlobalSideID))THEN
      ! Get global element index
      iGlobalElemID = SideInfo_Shared(SIDE_ELEMID,NonUniqueGlobalSideID)
      ! Get compute node element index of the side
      CNElemID = GetCNElemID(iGlobalElemID)
      ! Skip global elements, which are not on the compute node and not inside the halo region by checking if the CN element index
      IF(CNElemID.EQ.-1) CYCLE
      ! Check if any connected unique node ID is zero, which is impossible
      IF (ANY(NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(:,NonUniqueGlobalSideID).EQ.0)) THEN
        IPWRITE(*,*) 'NonUniqueGlobalSideID,NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(:,NonUniqueGlobalSideID):',&
                      NonUniqueGlobalSideID,NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(:,NonUniqueGlobalSideID)
        CALL abort(__STAMP__,'Wrong NonUniqueNodeID encountered in InitDepoSurfNodes()')
      END IF ! ANY()
    END IF
  END DO ! NonUniqueGlobalSideID = startVar,nVar

  ! Count the number of unique deposition nodes per processor
  nDepoSurfNodes = COUNT(IsDepoSurfNode)
  ! Count the number of unique deposition sides per processor
  nDepoSurfSides = COUNT(IsDepoSurfSide)
  ! DEALLOCATE(IsDepoSurfSide)

  ! Build Mappings between FEM vertices and surface deposition node IDs
  nDepoSurfNodesTotal = nDepoSurfNodes
  ALLOCATE(DepoSurfNodeID2FEMVertexID(1:nDepoSurfNodesTotal))
  DepoSurfNodeID2FEMVertexID = -1
  ALLOCATE(FEMVertexID2DepoSurfNodeID(1:nFEMVertices))
  FEMVertexID2DepoSurfNodeID = 0
  nDepoSurfNodesTotal = 0
  DO FEMVertexID=1, nFEMVertices
    IF (IsDepoSurfNode(FEMVertexID)) THEN
      nDepoSurfNodesTotal = nDepoSurfNodesTotal + 1
      DepoSurfNodeID2FEMVertexID(nDepoSurfNodesTotal) = FEMVertexID
      FEMVertexID2DepoSurfNodeID(FEMVertexID) = nDepoSurfNodesTotal
    ELSE
      FEMVertexID2DepoSurfNodeID(FEMVertexID) = -1
    END IF
  END DO

  ! Surface area associated with each deposition FEM vertex
  ALLOCATE(SurfNodeArea(1:nDepoSurfNodesTotal))
  SurfNodeArea = 0.
END IF ! InitializeSurfNodeArrays

! Build Vandermonde mapping from NodeType to NodeTypeVISU (equidistant with N=1)
CALL BuildSurfVdm(Nmax)

! FEATURE: [Surface Charge] Add contribution where inner BCs are used in combination with surface charging
! Loop over all boundary condition sides
DO SideID=1,nBCSides
  ! Get BC type
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  SELECT CASE(BCType)
#if USE_HDG
  CASE(HDGDIRICHLETBCSIDEIDS) ! HDG Dirichlet BC Side IDs: BCType = BoundaryType(BC(SideID),BC_TYPE)
    ! Skip
  CASE(10,11,12) !Neumann,
    ! Skip
  CASE(20) ! Conductor: Floating Boundary Condition (FPC)
    ! Skip
  CASE(30) ! Distributed Capacitance
    ! Sum up all contributions of surface area to the FEM vertices of each side
    CALL Buildpq2iNode(SideID,SubSideAreaEquiN1)
    IF(.NOT.DoRestart) CALL CalculateSurfNodeArea(SideID)
  CASE(40) ! Circuit model boundary condition
    ! Skip
#else
  CASE DEFAULT ! unknown BCType
    CALL CollectiveStop(__STAMP__,' unknown BC Type in hdg.f90!',IntInfo=BCType)
#endif /*USE_HDG*/
  END SELECT ! BCType
END DO
DEALLOCATE(Vdm_N_EQ)

! DEALLOCATE(IsDepoSurfNode)
#if USE_MPI
CALL InitDepoSurfNodesMPI() ! Initialize MPI communicator for surface node communication
! Collect surface area contributions SurfNodeArea(iDepoSurfNodeID) of all processes on the MPIRoot
IF(.NOT.PerformLoadBalance) CALL CollectSurfNodeAreaOnMPIRoot()
! Initialize the the SurfNodeArea(iDepoSurfNodeID) container on all processes except MPIRoot, which distribtues the data to all others
CALL ReverseExchangeSurfNodeArea()
#endif /*USE_MPI*/
! All processes, except the MPIRoot re-allocate the array during load balance and the MPIRoot sends each process the surface charge
! they need
IF (InitializeSurfNodeArrays) THEN
  ALLOCATE(SurfNodeSource(1:nDepoSurfNodesTotal))
  SurfNodeSource=0.0
END IF ! InitializeSurfNodeArrays
#if USE_LOADBALANCE
! MPIRoot sends SurfNodeSource to all processes
IF(PerformLoadBalance) CALL LBReverseExchangeSurfNodeSource()
#endif /*USE_LOADBALANCE*/

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!',DisplayLine=.FALSE.)

InitDepoSurfNodesIsDone = .TRUE.

END SUBROUTINE InitDepoSurfNodes


!===================================================================================================================================
!> Build Vandermondes
!> a) Vdm_EQ_N: Vandermonde for mapping from N=1 (equidistant) to N=Nloc (Gauss/Gauss-Lobatto)
!> b) Vdm_N_EQ: Map G/GL (current node type) to equidistant distribution with N=1
!===================================================================================================================================
SUBROUTINE BuildSurfVdm(Nmax)
! MODULES
USE MOD_Preproc
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_PICDepo_Vars       ,ONLY: Vdm_EQ_N,Vdm_N_EQ
USE MOD_Interpolation_Vars ,ONLY: Nmin
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVISU,NodeType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Nmax !< Maximum polynomial degree vor Vandermonde mappings
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: Nloc
!===================================================================================================================================
! Build Vandermonde Vdm_EQ_N for mapping from N=1 (equidistant) to N=Nloc (Gauss/Gauss-Lobatto)
ALLOCATE(Vdm_EQ_N(Nmin:Nmax))
DO Nloc = Nmin, Nmax
  ALLOCATE(Vdm_EQ_N(Nloc)%Vdm(0:Nloc,0:1))
  CALL GetVandermonde(1, NodeTypeVISU, Nloc, NodeType, Vdm_EQ_N(Nloc)%Vdm(0:Nloc,0:1), modal=.FALSE.)
END DO ! Nloc = Nmin, Nmax

! Build Vandermonde Vdm_N_EQ for mapping G/GL (current node type) to equidistant distribution with N=1
ALLOCATE(Vdm_N_EQ(Nmin:Nmax))
DO Nloc = Nmin, Nmax
  ! Allocate and determine Vandermonde mapping from NodeType to equidistant (visu) node set
  ALLOCATE(Vdm_N_EQ(Nloc)%Vdm(0:1,0:Nloc))
  CALL GetVandermonde(Nloc, NodeType, 1, NodeTypeVISU, Vdm_N_EQ(Nloc)%Vdm(0:1,0:Nloc), modal=.FALSE.)
  ! Required only for integration
  ALLOCATE(Vdm_N_EQ(Nloc)%xIP_VISU(0:Nloc))
  ALLOCATE(Vdm_N_EQ(Nloc)%wIP_VISU(0:Nloc))
  CALL GetNodesAndWeights(Nloc, NodeTypeVISU, xIP = Vdm_N_EQ(Nloc)%xIP_VISU, wIP = Vdm_N_EQ(Nloc)%wIP_VISU)
END DO ! Nloc = Nmin, Nmax

END SUBROUTINE BuildSurfVdm


!===================================================================================================================================
!> Determine the surface area associated with a FEM vertex
!===================================================================================================================================
SUBROUTINE Buildpq2iNode(SideID,SubSideAreaEquiN1)
! MODULES
#if USE_MPI
USE MOD_Globals            ,ONLY: myrank
#endif /*USE_MPI*/
USE MOD_Preproc
USE MOD_Globals            ,ONLY: UNIT_stdOut,abort,VECNORM3D
USE MOD_PICDepo_Vars       ,ONLY: IsDepoSurfSide,Vdm_N_EQ,pq2iNode
USE MOD_Interpolation_Vars ,ONLY: N_Inter,Nmax
USE MOD_Mesh_Vars          ,ONLY: N_SurfMesh
USE MOD_Mesh_Vars          ,ONLY: SideToNonUniqueGlobalSide
USE MOD_Mesh_Vars          ,ONLY: NonUniqueGlobalSideIDToNonUniqueGlobalNodeID
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID !< Local side index
REAL,INTENT(OUT)   :: SubSideAreaEquiN1(0:1,0:1) !< Side areas
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,p,q,Nloc,NonUniqueGlobalSideID,NonUniqueNodeID,iERROR,i
INTEGER,PARAMETER:: NSideN1 = 1 ! Set polynomial for equidistant basis
REAL    :: SideAreaNloc,SideAreaEquiN1
REAL    :: SurfElemEquiN1(0:1,0:1),Face_xGPEquiN1(3,0:1,0:1)
REAL    :: tmp(1:3,0:Nmax,0:Nmax),tmp2(1:3,0:Nmax,0:Nmax)
!===================================================================================================================================
! Get non-unique global side index from local side index
NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,SideID) ! Get global side index

! Get polynomial degree of side (can be inner side or boundary side)
Nloc = N_SurfMesh(SideID)%NSide

! IPWRITE(*,*) 'Nloc,NSideN1:', Nloc,NSideN1
! Get SurfElemEquiN1: Surface area elements on equidistant basis with N=1
! Check if the polynomial degree is one
IF(Nloc.EQ.NSideN1)THEN
  ! Store in temp array for switching from NodeType to NodeTypeVISU
  tmp2(1,0:NSideN1,0:NSideN1) = N_SurfMesh(SideID)%SurfElem(0:NSideN1,0:NSideN1)
ELSE
  ! From high to low
  ! Transform side keeping the same degree: switch to Legendre basis
  CALL ChangeBasis2D(1, Nloc, Nloc, N_Inter(Nloc)%sVdm_Leg, N_SurfMesh(SideID)%SurfElem(0:Nloc,0:Nloc) ,&
                                                                                  tmp(1,0:Nloc,0:Nloc) )
  ! Switch back to nodal basis
  CALL ChangeBasis2D(1, NSideN1, NSideN1, N_Inter(NSideN1)%Vdm_Leg , tmp(1,0:NSideN1,0:NSideN1) ,&
                                                                    tmp2(1,0:NSideN1,0:NSideN1) )
END IF ! Nloc.EQ.NSideN1

! Swtich from NodeType (N=1) to NodeTypeVISU (N=1)
CALL ChangeBasis2D(1, 1, 1, Vdm_N_EQ(NSideN1)%Vdm, tmp2(1,0:NSideN1,0:NSideN1), SurfElemEquiN1(0:NSideN1,0:NSideN1) )

! Sanity check: Make sure that the sum of the sub areas does not change when switching from G/GL with N=Nloc to equidistant with N=1
! Calculate side area: G/GL with N=Nloc
SideAreaNloc = 0
DO q=0,Nloc; DO p=0,Nloc
  SideAreaNloc = SideAreaNloc + N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)*N_SurfMesh(SideID)%SurfElem(p,q)
END DO; END DO ! p,q

! Calculate side area: equidistant visu nodes with N=1
SideAreaEquiN1 = 0
DO q=0,NSideN1; DO p=0,NSideN1
  SubSideAreaEquiN1(p,q) = Vdm_N_EQ(NSideN1)%wIP_VISU(p)*Vdm_N_EQ(NSideN1)%wIP_VISU(q)*SurfElemEquiN1(p,q)
  SideAreaEquiN1 = SideAreaEquiN1 + SubSideAreaEquiN1(p,q)
END DO; END DO ! p,q

! Sanity check: SideAreaNloc and SideAreaEquiN1 must only differ relatively by 1e-3
IF (.NOT.ALMOSTEQUALRELATIVE(SideAreaNloc, SideAreaEquiN1, 1e-3)) THEN
  IPWRITE(*,*) 'Error in area calculation for surface charge deposition with N=1 and linear weighting'
  IPWRITE(*,*) 'SideID,NonUniqueGlobalSideID,Nloc,SideAreaNloc  :', SideID,NonUniqueGlobalSideID,Nloc,SideAreaNloc
  IPWRITE(*,*) 'SideID,NonUniqueGlobalSideID,Nloc,SideAreaEquiN1:', SideID,NonUniqueGlobalSideID,NSideN1,SideAreaEquiN1
  CALL abort(__STAMP__,' Sum of the sub areas has changed  when switching from G/GL with N=Nloc to equidistant with N=1', IERROR)
END IF ! .NOT.ALMOSTEQUALRELATIVE(SideAreaNloc, SideAreaEquiN1, 1edd-3)

! Get Face_xGPEquiN1: xGP mapped to equidistant basis with N=1
! Check if the polynomial degree is one
IF(Nloc.EQ.NSideN1)THEN
  ! Store in temp array for switching from NodeType to NodeTypeVISU
  tmp2(1:3,0:NSideN1,0:NSideN1) = N_SurfMesh(SideID)%Face_xGP(1:3,0:NSideN1,0:NSideN1)
ELSE
  ! From high to low
  ! Transform side keeping the same degree: switch to Legendre basis
  CALL ChangeBasis2D(3, Nloc, Nloc, N_Inter(Nloc)%sVdm_Leg, N_SurfMesh(SideID)%Face_xGP(1:3,0:Nloc,0:Nloc) ,&
                                                                                    tmp(1:3,0:Nloc,0:Nloc) )
  ! Switch back to nodal basis
  CALL ChangeBasis2D(3, NSideN1, NSideN1, N_Inter(NSideN1)%Vdm_Leg , tmp(1:3,0:NSideN1,0:NSideN1) ,&
                                                                    tmp2(1:3,0:NSideN1,0:NSideN1) )
END IF ! Nloc.EQ.NSideN1

! Swtich from NodeType (N=1) to NodeTypeVISU (N=1)
CALL ChangeBasis2D(3, 1, 1, Vdm_N_EQ(NSideN1)%Vdm, tmp2(1:3,0:NSideN1,0:NSideN1), Face_xGPEquiN1(1:3,0:NSideN1,0:NSideN1) )

! Map surface charge from vertices to SideID surface with N=1
! Note that the loop runs in the p-q-oriented system
DO q=0,1; DO p=0,1
  ! Get local node index by checking the distance of the four cornder nodes
  CALL GetClosestNode(NonUniqueGlobalSideID,Face_xGPEquiN1(1:3,p,q),iNode)
! IPWRITE(*,*) 'NonUniqueGlobalSideID,p,q,iNode:', NonUniqueGlobalSideID,p,q,iNode
! IPWRITE(*,*) 'Face_xGPEquiN1(1:3,p,q)                             :', Face_xGPEquiN1(1:3,p,q)
! IPWRITE(*,*) 'BezierControlPoints3D(1:3,p,q,NonUniqueGlobalSideID):', BezierControlPoints3D(1:3,p,q,NonUniqueGlobalSideID)
  ! Set mapping
  pq2iNode(p,q,SideID) = iNode
  ! IPWRITE(*,*) 'p,q,iNode,2*q + p + 1:', p,q,iNode,2*q + p + 1
  ! Mapping from non-unique global side index to non-unique global node index
  NonUniqueNodeID = NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(iNode,NonUniqueGlobalSideID)
! IPWRITE(*,*) 'NodeCoords_Shared(1:3,NonUniqueNodeID)              :', NodeCoords_Shared(1:3,NonUniqueNodeID)
  ! Sanity check
  IF (NonUniqueNodeID.LE.0) THEN
    IPWRITE(*,*) 'NonUniqueNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID:',&
                  NonUniqueNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID
    CALL abort(__STAMP__,' NonUniqueNodeID <= 0')
  END IF ! NonUniqueNodeID.LE.0
END DO; END DO ! q=0,1; DO p=0,1

! Sanity check: Make sure that each node (1,2,3 and 4) occur exactly once
i=0
DO q=0,1; DO p=0,1
  i=i+1
  IF (COUNT(pq2iNode(:,:,SideID).EQ.i).GT.1) CALL abort(__STAMP__,' Error in Buildpq2iNode. Equal indices in pq2iNode(:,:,SideID)')
  IF (COUNT(pq2iNode(:,:,SideID).EQ.i).LT.1) CALL abort(__STAMP__,' Error in Buildpq2iNode. Index missing in pq2iNode(:,:,SideID)')
END DO; END DO ! q=0,1; DO p=0,clas2vtk1

END SUBROUTINE Buildpq2iNode


!===================================================================================================================================
!> Determine the surface area associated with a FEM vertex
!===================================================================================================================================
SUBROUTINE CalculateSurfNodeArea(SideID)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals            ,ONLY: UNIT_stdOut,abort,VECNORM3D,myrank
USE MOD_PICDepo_Vars       ,ONLY: FEMVertexID2DepoSurfNodeID,IsDepoSurfSide,SurfNodeArea,pq2iNode
USE MOD_Mesh_Vars          ,ONLY: SideToNonUniqueGlobalSide
USE MOD_Mesh_Vars          ,ONLY: N_SurfMesh
USE MOD_Mesh_Vars          ,ONLY: NonUniqueGlobalSideIDToNonUniqueGlobalNodeID,NonUniqueGlobalNodeIDToFEMVertexID
USE MOD_Interpolation_Vars ,ONLY: N_Inter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: SideID   !< Local side index
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,FEMVertexID,iDepoSurfNodeID,p,q,NonUniqueGlobalSideID,NonUniqueNodeID
INTEGER :: NodIndx(4), k, j, Nloc
!===================================================================================================================================
! Get non-unique global side index from local side index
NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,SideID) ! Get global side index
! Map surface charge from vertices to SideID surface with N=1
! Note that the loop runs in the p-q-oriented system
DO q=0,1; DO p=0,1
  ! Get local node index by checking the distance of the four cornder nodes
  ! Set mapping p,q -> iNode
  iNode = pq2iNode(p,q,SideID)
  ! IPWRITE(*,*) 'p,q,iNode,2*q + p + 1:', p,q,iNode,2*q + p + 1
  ! Mapping from non-unique global side index to non-unique global node index
  NonUniqueNodeID = NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(iNode,NonUniqueGlobalSideID)
  ! Sanity check
  IF (NonUniqueNodeID.LE.0) THEN
    IPWRITE(*,*) 'NonUniqueNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID:',&
                  NonUniqueNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID
    CALL abort(__STAMP__,' NonUniqueNodeID <= 0')
  END IF ! NonUniqueNodeID.LE.0
  ! Mapping from NonUniqueNodeID to FEMVertexID
  FEMVertexID = NonUniqueGlobalNodeIDToFEMVertexID(NonUniqueNodeID)
  ! Get surface deposition node index
  iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
  ! Add contribution to the FEM vertex (note that double periodicity collapses all four corner nodes into a single vertex index)
  !SurfNodeArea(iDepoSurfNodeID) = SurfNodeArea(iDepoSurfNodeID) + SubSideAreaEquiN1(p,q)
  ! Note that on inner BC this might not work because the side is not always oriented in the master ordering: iNode = 2*q + p + 1
  NodIndx(2*q + p + 1) = iDepoSurfNodeID
END DO; END DO ! q=0,1; DO p=0,1

! Note that on inner BC this might not work because the side is not always oriented in the master ordering
Nloc = N_SurfMesh(SideID)%NSide
DO j=0,Nloc;DO k=0,Nloc
  ASSOCIATE(xGP => N_Inter(Nloc)%xGP, wGP => N_Inter(Nloc)%wGP)
    ! CVWV cannot be accessed here with "0" because of the associate construct!
    SurfNodeArea(NodIndx(2)) = SurfNodeArea(NodIndx(2)) + wGP(j)*wGP(k)*( (1.-xGP(j)) * (1.-xGP(k) )*N_SurfMesh(SideID)%SurfElem(j,k) )/4.
    SurfNodeArea(NodIndx(1)) = SurfNodeArea(NodIndx(1)) + wGP(j)*wGP(k)*( (1.-xGP(j)) * (1.+xGP(k) )*N_SurfMesh(SideID)%SurfElem(j,k) )/4.
    SurfNodeArea(NodIndx(4)) = SurfNodeArea(NodIndx(4)) + wGP(j)*wGP(k)*( (1.+xGP(j)) * (1.-xGP(k) )*N_SurfMesh(SideID)%SurfElem(j,k) )/4.
    SurfNodeArea(NodIndx(3)) = SurfNodeArea(NodIndx(3)) + wGP(j)*wGP(k)*( (1.+xGP(j)) * (1.+xGP(k) )*N_SurfMesh(SideID)%SurfElem(j,k) )/4.
  END ASSOCIATE
END DO; END DO

END SUBROUTINE CalculateSurfNodeArea


!===================================================================================================================================
!> Determine the closest NodeCoords node closest to the vector x(1:3) and return the node index (1, 2, 3 or 4)
!===================================================================================================================================
SUBROUTINE GetClosestNode(NonUniqueGlobalSideID,x,NodeIndex)
! MODULES
USE MOD_Preproc
USE MOD_Globals            ,ONLY: UNIT_stdOut,abort,VECNORM3D
USE MOD_Mesh_Vars          ,ONLY: NonUniqueGlobalSideIDToNonUniqueGlobalNodeID
USE MOD_Particle_Mesh_Vars ,ONLY: NodeCoords_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: NonUniqueGlobalSideID !< Non-unique global side index
REAL,INTENT(IN)     :: x(1:3)                !< Coordinate of the face interpolation point on the mesh surface corners
INTEGER,INTENT(OUT) :: NodeIndex             !< Node index (1 to 4) of the NodeCoords node closest to the vector x(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,NonUniqueNodeID
REAL    :: norm,PartDistDepo(4)
!===================================================================================================================================
! Loop over the four corner nodes
NODELOOP: DO iNode = 1, 4
  ! Get the non-unique node index
  NonUniqueNodeID = NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(iNode,NonUniqueGlobalSideID)
  ! IPWRITE(*,*) 'NonUniqueNodeID:', NonUniqueNodeID
  ! Sanity check
  IF(NonUniqueNodeID.LE.0) CALL abort(__STAMP__,'Wrong NonUniqueNodeID encountered in surface charge deposition init')
  ! Calculate the distance
  norm = VECNORM3D(NodeCoords_Shared(1:3,NonUniqueNodeID) - x(1:3))
  ! IPWRITE(*,*) 'iNode,NodeCoords_Shared(1:3,NonUniqueNodeID) - x(1:3),norm:',&
  !               iNode,NodeCoords_Shared(1:3,NonUniqueNodeID) - x(1:3),norm
  ! Check if the distance is greater than zero
  IF(norm.GT.0.)THEN
    PartDistDepo(iNode) = 1./norm
  ELSE
    PartDistDepo(:) = 0.
    PartDistDepo(iNode) = 1.0
    EXIT NODELOOP
  END IF ! norm.GT.0.
END DO NODELOOP ! iNode = 1, 4

! Get index of maximum location
NodeIndex = MAXLOC(PartDistDepo,DIM=1)
! IPWRITE(*,*) 'PartDistDepo,NodeIndex:', PartDistDepo,NodeIndex

END SUBROUTINE GetClosestNode


!===================================================================================================================================
!> Calculate the shape function radius for each element depending on the neighbouring element sizes and the own element size
!===================================================================================================================================
SUBROUTINE InitShapeFunctionAdaptive()
! MODULES
USE MOD_Preproc
USE MOD_Globals                     ,ONLY: UNIT_stdOut,abort
USE MOD_PICDepo_Vars                ,ONLY: SFAdaptiveDOF,SFAdaptiveSmoothing,SFElemr2_Shared,dim_sf,dimFactorSF
USE MOD_ReadInTools                 ,ONLY: GETREAL,GETLOGICAL
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,BoundsOfElem_Shared
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID
USE MOD_Globals_Vars                ,ONLY: PI
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemMidPoint_Shared,ElemToElemMapping,ElemToElemInfo
USE MOD_Mesh_Tools                  ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars          ,ONLY: NodeCoords_Shared
USE MOD_PICDepo_Shapefunction_Tools ,ONLY: SFNorm
USE MOD_Interpolation_Vars          ,ONLY: Nmax
USE MOD_DG_Vars                     ,ONLY: N_DG_Mapping
#if USE_MPI
USE MOD_Globals                     ,ONLY: IERROR
USE MOD_PICDepo_Vars                ,ONLY: SFElemr2_Shared_Win
USE MOD_Globals                     ,ONLY: MPIRoot
USE MOD_MPI_Shared_Vars             ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank,MPI_COMM_SHARED
USE MOD_MPI_Shared
#else
USE MOD_Mesh_Vars                   ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars            ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: UniqueNodeID,NonUniqueNodeID,iNode,NeighUniqueNodeID,GlobalElemID
REAL                           :: SFDepoScaling
LOGICAL                        :: ElemDone
INTEGER                        :: ppp,globElemID
REAL                           :: r_sf_tmp,SFAdaptiveDOFDefault,DOFMax
INTEGER                        :: iCNElem,firstElem,lastElem,jNode,NbElemID,NeighNonUniqueNodeID,minN_PP
CHARACTER(32)                  :: hilf2,hilf3
#if USE_MPI
#endif /*USE_MPI*/
REAL                           :: CharacteristicLength,BoundingBoxVolume
!===================================================================================================================================
! Set the number of DOF/SF
! Check which shape function dimension is used and set default value
SELECT CASE(dim_sf)
CASE(1)
  SFAdaptiveDOFDefault=2.0*(1.+1.)
  DOFMax = 2.0*(REAL(Nmax)+1.) ! Max. DOF per element in 1D
  hilf2 = '2*(N+1)'            ! Max. DOF per element in 1D for abort message
CASE(2)
  SFAdaptiveDOFDefault=PI*(1.+1.)**2
  DOFMax = PI*(REAL(Nmax)+1.)**2 ! Max. DOF per element in 2D
  hilf2 = 'PI*(N+1)**2'          ! Max. DOF per element in 1D for abort message
CASE(3)
  SFAdaptiveDOFDefault=(4./3.)*PI*(1.+1.)**3
  DOFMax = (4./3.)*PI*(REAL(Nmax)+1.)**3 ! Max. DOF per element in 2D
  hilf2 = '(4/3)*PI*(N+1)**3'            ! Max. DOF per element in 1D for abort message
END SELECT
WRITE(UNIT=hilf3,FMT='(G0)') SFAdaptiveDOFDefault
SFAdaptiveDOF = GETREAL('PIC-shapefunction-adaptive-DOF',TRIM(hilf3))

LBWRITE(UNIT_StdOut,'(A,F10.2)') "         PIC-shapefunction-adaptive-DOF :", SFAdaptiveDOF
LBWRITE(UNIT_StdOut,'(A,A19,A,F10.2,A,I0,A)') " Maximum allowed is ",TRIM(hilf2)," :", DOFMax," (calculated for N=",Nmax,")"
LBWRITE(UNIT_StdOut,*) "Set a value lower or equal to than the maximum for a given polynomial degree N\n"
LBWRITE(UNIT_StdOut,*) "              N:     1      2      3      4      5       6       7"
LBWRITE(UNIT_StdOut,*) "  ----------------------------------------------------------------"
LBWRITE(UNIT_StdOut,*) "           | 1D:     4      6      8     10     12      14      16"
LBWRITE(UNIT_StdOut,*) "  Max. DOF | 2D:    12     28     50     78    113     153     201"
LBWRITE(UNIT_StdOut,*) "           | 3D:    33    113    268    523    904    1436    2144"
LBWRITE(UNIT_StdOut,*) "  ----------------------------------------------------------------"

IF(SFAdaptiveDOF.GT.DOFMax)THEN
  SWRITE(UNIT_StdOut,*) "Reduce the number of DOF/SF in order to have no DOF outside of the deposition range (neighbour elems)"
  CALL abort(__STAMP__,'PIC-shapefunction-adaptive-DOF > '//TRIM(hilf2)//' is not allowed')
ELSE
  ! Check which shape function dimension is used
  SELECT CASE(dim_sf)
  CASE(1)
    SFDepoScaling  = SFAdaptiveDOF/2.0
  CASE(2)
    SFDepoScaling  = SQRT(SFAdaptiveDOF/PI)
  CASE(3)
    SFDepoScaling  = (3.*SFAdaptiveDOF/(4.*PI))**(1./3.)
  END SELECT
END IF

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

CALL Allocate_Shared((/2,nComputeNodeTotalElems/),SFElemr2_Shared_Win,SFElemr2_Shared)
CALL MPI_WIN_LOCK_ALL(0,SFElemr2_Shared_Win,IERROR)
#else
ALLOCATE(SFElemr2_Shared(1:2,1:nElems))
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif
  SFElemr2_Shared = HUGE(1.)
#if USE_MPI
END IF
CALL BARRIER_AND_SYNC(SFElemr2_Shared_Win,MPI_COMM_SHARED)
#endif
DO iCNElem = firstElem,lastElem
  ElemDone = .FALSE.
  minN_PP = N_DG_Mapping(2,GetGlobalElemID(iCNElem))
  DO ppp = 1,ElemToElemMapping(2,iCNElem)
    ! Get neighbour global element ID
    globElemID = GetGlobalElemID(ElemToElemInfo(ElemToElemMapping(1,iCNElem)+ppp))
    NbElemID = GetCNElemID(globElemID)
    minN_PP = MIN(minN_PP, N_DG_Mapping(2,globElemID))
    ! Loop neighbour nodes
    Nodeloop: DO jNode = 1, 8
      NeighNonUniqueNodeID = ElemNodeID_Shared(jNode,NbElemID)
      NeighUniqueNodeID = NodeInfo_Shared(NeighNonUniqueNodeID)
      ! Loop my nodes
      DO iNode = 1, 8
        NonUniqueNodeID = ElemNodeID_Shared(iNode,iCNElem)
        UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
        IF (UniqueNodeID.EQ.NeighUniqueNodeID) CYCLE Nodeloop ! Skip coinciding nodes of my and my neighbours element
      END DO
      ElemDone =.TRUE.

      ! Measure distance from my corner nodes to neighbour elem corner nodes
      DO iNode = 1, 8
        NonUniqueNodeID = ElemNodeID_Shared(iNode,iCNElem)

        ! Only measure distances in the dimension in which the nodes to not coincide (i.e. they are not projected onto each other
        ! in 1D or 2D deposition)
        ASSOCIATE( v1 => NodeCoords_Shared(1:3,NonUniqueNodeID)      ,&
                   v2 => NodeCoords_Shared(1:3,NeighNonUniqueNodeID) )
          IF(SFMeasureDistance(v1,v2)) THEN
            r_sf_tmp = SFNorm(v1-v2)
            IF (r_sf_tmp.LT.SFElemr2_Shared(1,iCNElem)) SFElemr2_Shared(1,iCNElem) = r_sf_tmp
          END IF ! SFMeasureDistance(v1,v2)
        END ASSOCIATE
      END DO ! iNode = 1, 8

    END DO Nodeloop
  END DO

  ! Sanity check if no neighbours are present
  IF (.NOT.ElemDone) THEN
    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iCNElem)
      r_sf_tmp = SFNorm(ElemMidPoint_Shared(1:3,iCNElem)-NodeCoords_Shared(1:3,NonUniqueNodeID))
      IF (r_sf_tmp.LT.SFElemr2_Shared(1,iCNElem)) SFElemr2_Shared(1,iCNElem) = r_sf_tmp
    END DO
  END IF

  ! Because ElemVolume_Shared(CNElemID) is not available for halo elements, the bounding box volume is used as an approximate
  ! value for the element volume from which the characteristic length of the element is calculated
  GlobalElemID = GetGlobalElemID(iCNElem)
  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,GlobalElemID) ) ! 1-2: Min, Max value; 1-3: x,y,z
    BoundingBoxVolume = (Bounds(2,1)-Bounds(1,1)) * (Bounds(2,2)-Bounds(1,2)) * (Bounds(2,3)-Bounds(1,3))
  END ASSOCIATE
  ! Check which shape function dimension is used
  SELECT CASE(dim_sf)
  CASE(1)
    !CharacteristicLength = ElemVolume_Shared(iCNElem) / dimFactorSF
    CharacteristicLength = BoundingBoxVolume / dimFactorSF
  CASE(2)
    !CharacteristicLength = SQRT(ElemVolume_Shared(iCNElem) / dimFactorSF)
    CharacteristicLength = SQRT(BoundingBoxVolume / dimFactorSF)
  CASE(3)
    !CharacteristicLength = ElemCharLength_Shared(iCNElem)
    CharacteristicLength = BoundingBoxVolume**(1./3.)
  END SELECT

  ! Check characteristic length of cell (or when using SFAdaptiveSmoothing)
  IF(CharacteristicLength.LT.SFElemr2_Shared(1,iCNElem).OR.SFAdaptiveSmoothing)THEN
    SFElemr2_Shared(1,iCNElem) = (SFElemr2_Shared(1,iCNElem) + CharacteristicLength)/2.0
  END IF

  ! Scale the radius so that it reaches at most the neighbouring cells but no further (all neighbours of the 8 corner nodes)
  SFElemr2_Shared(1,iCNElem) = SFElemr2_Shared(1,iCNElem) * SFDepoScaling / (minN_PP+1.)
  SFElemr2_Shared(2,iCNElem) = SFElemr2_Shared(1,iCNElem)**2

  ! Sanity checks
  IF(SFElemr2_Shared(1,iCNElem).LE.0.0)      CALL abort(__STAMP__,'Shape function radius <= zero!')
  IF(SFElemr2_Shared(1,iCNElem).GE.HUGE(1.)) CALL abort(__STAMP__,'Shape function radius >= HUGE(1.)!')
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(SFElemr2_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE InitShapeFunctionAdaptive


!===================================================================================================================================
!> Fill PeriodicSFCaseMatrix when using shape function deposition in combination with periodic boundaries
!===================================================================================================================================
SUBROUTINE InitPeriodicSFCaseMatrix()
! MODULES
USE MOD_Globals            ,ONLY: UNIT_StdOut
#if USE_MPI
USE MOD_Globals            ,ONLY: MPIRoot
#endif /*USE_MPI*/
USE MOD_Particle_Mesh_Vars ,ONLY: PeriodicSFCaseMatrix,NbrOfPeriodicSFCases
USE MOD_PICDepo_Vars       ,ONLY: dim_sf,dim_periodic_vec1,dim_periodic_vec2,dim_sf_dir1,dim_sf_dir2
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_ReadInTools        ,ONLY: PrintOption
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: I,J
!===================================================================================================================================
IF (GEO%nPeriodicVectors.LE.0) THEN

  ! Set defaults and return in non-periodic case
  NbrOfPeriodicSFCases = 0
  ALLOCATE(PeriodicSFCaseMatrix(1:1,1:3))
  PeriodicSFCaseMatrix(:,:) = 0

ELSE

  ! Build case matrix:
  ! Particles may move in more periodic directions than their charge is deposited, e.g., fully periodic in combination with
  ! 1D shape function
  NbrOfPeriodicSFCases = 3**dim_sf

  ALLOCATE(PeriodicSFCaseMatrix(1:NbrOfPeriodicSFCases,1:3))
  PeriodicSFCaseMatrix(:,:) = 0
  IF (dim_sf.EQ.1) THEN
    PeriodicSFCaseMatrix(1,1) = 1
    PeriodicSFCaseMatrix(3,1) = -1
  END IF
  IF (dim_sf.EQ.2) THEN
    PeriodicSFCaseMatrix(1:3,1) = 1
    PeriodicSFCaseMatrix(7:9,1) = -1
    DO I = 1,3
      PeriodicSFCaseMatrix(I*3-2,2) = 1
      PeriodicSFCaseMatrix(I*3,2) = -1
    END DO
  END IF
  IF (dim_sf.EQ.3) THEN
    PeriodicSFCaseMatrix(1:9,1) = 1
    PeriodicSFCaseMatrix(19:27,1) = -1
    DO I = 1,3
      PeriodicSFCaseMatrix(I*9-8:I*9-6,2) = 1
      PeriodicSFCaseMatrix(I*9-2:I*9,2) = -1
      DO J = 1,3
        PeriodicSFCaseMatrix((J*3-2)+(I-1)*9,3) = 1
        PeriodicSFCaseMatrix((J*3)+(I-1)*9,3) = -1
      END DO
    END DO
  END IF

  ! Define which of the periodic vectors are used for 2D shape function and display info
  IF(dim_sf.EQ.2)THEN
    IF(GEO%nPeriodicVectors.EQ.1)THEN
      dim_periodic_vec1 = 1
      dim_periodic_vec2 = 0
    ELSEIF(GEO%nPeriodicVectors.EQ.2)THEN
      dim_periodic_vec1 = 1
      dim_periodic_vec2 = 2
    ELSEIF(GEO%nPeriodicVectors.EQ.3)THEN
      dim_periodic_vec1 = dim_sf_dir1
      dim_periodic_vec2 = dim_sf_dir2
    END IF ! GEO%nPeriodicVectors.EQ.1
    CALL PrintOption('Dimension of 1st periodic vector for 2D shape function','INFO',IntOpt=dim_periodic_vec1)
    LBWRITE(UNIT_StdOut,*) "1st PeriodicVector =", GEO%PeriodicVectors(1:3,dim_periodic_vec1)
    CALL PrintOption('Dimension of 2nd periodic vector for 2D shape function','INFO',IntOpt=dim_periodic_vec2)
    LBWRITE(UNIT_StdOut,*) "2nd PeriodicVector =", GEO%PeriodicVectors(1:3,dim_periodic_vec2)
  END IF ! dim_sf.EQ.2

END IF

END SUBROUTINE InitPeriodicSFCaseMatrix


!===================================================================================================================================
!> Fill PeriodicSFCaseMatrix when using shape function deposition in combination with periodic boundaries
!===================================================================================================================================
SUBROUTINE InitAxisymmetrySF()
! MODULES
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_Symmetry_Vars          ,ONLY: Symmetry
USE MOD_Particle_Mesh_Vars     ,ONLY: AxisymmetricSF
USE MOD_ReadInTools            ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPartBound
!===================================================================================================================================
AxisymmetricSF = .FALSE.
IF (Symmetry%Axisymmetric) THEN
  ! Check boundaries for symmetry axis
  DO iPartBound=1,nPartBound
    IF(PartBound%TargetBoundCond(iPartBound).EQ.PartBound%SymmetryAxis) THEN
      AxisymmetricSF = .TRUE.
      CALL PrintOption('Found symmetry axis for shape function deposition','INFO',LogOpt=AxisymmetricSF)
      RETURN
    END IF
  END DO
END IF
END SUBROUTINE InitAxisymmetrySF


SUBROUTINE Deposition(doParticle_In, stage_opt, skipVerifyCharge_opt)
!============================================================================================================================
! This subroutine performs the deposition of the particle charge and current density to the grid
! following list of distribution methods are implemented
! - shape function       (only one type implemented)
! useVMPF added, therefore, this routine contains automatically the use of variable mpfs
!============================================================================================================================
! USE MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Analyze_Vars ,ONLY: DoVerifyCharge,PartAnalyzeStep
USE MOD_Particle_Vars
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Method        ,ONLY: DepositionMethod
USE MOD_PIC_Analyze           ,ONLY: VerifyDepositedCharge
USE MOD_TimeDisc_Vars         ,ONLY: iter
USE MOD_Mesh_Vars             ,ONLY: nElems
#if USE_MPI
USE MOD_MPI_Shared            ,ONLY: BARRIER_AND_SYNC
#endif  /*USE_MPI*/
#if USE_HDG
USE MOD_HDG_Vars              ,ONLY: HDGSkip, HDGSkipInit, HDGSkip_t0
USE MOD_TimeDisc_Vars         ,ONLY: time
#endif  /*USE_HDG*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT variable declaration
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength) ! Marked particles for deposition
INTEGER,INTENT(IN),OPTIONAL   :: stage_opt ! Input: Runge-Kutta stage info
LOGICAL,INTENT(IN),OPTIONAL   :: skipVerifyCharge_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration
INTEGER                       :: stage,iElem
LOGICAL                       :: skipVerifyCharge
!===================================================================================================================================
! Return, if no deposition is required
IF(.NOT.DoDeposition) RETURN

! Set stage when using Runge-Kutta methods
IF (PRESENT(stage_opt)) THEN
  stage = stage_opt
ELSE
  stage = 0
END IF

#if USE_HDG
! HDGSkip Check whether the deposition should be skipped in this iteration
IF (iter.GT.0 .AND. HDGSkip.NE.0) THEN
  IF (time.LT.HDGSkip_t0) THEN
    IF (MOD(iter,INT(HDGSkipInit,8)).NE.0) RETURN
  ELSE
    IF (MOD(iter,INT(HDGSkip,8)).NE.0) RETURN
  END IF
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  IF (stage.GT.1) THEN
    RETURN
  END IF
#endif
END IF
#endif /*USE_HDG*/

! Nullify the charge container
IF((stage.EQ.0).OR.(stage.EQ.1))THEN
  DO iElem = 1, nElems
    PS_N(iElem)%PartSource = 0.0
  END DO ! iElem = 1, nElems
END IF

! Check whether only specific particles are to be deposited using the logical vector doParticle_In(1:PDM%ParticleVecLength)
IF(PRESENT(doParticle_In)) THEN
  CALL DepositionMethod(doParticle_In, stage_opt=stage)
ELSE
  CALL DepositionMethod(stage_opt=stage)
END IF

! Deposited charge verification: Compare the deposited charge in PS_N(:)%PartSource(4,:,:,:) with the sum of all particles charges
! divided by the system volume
IF((stage.EQ.0).OR.(stage.EQ.4)) THEN
  IF(MOD(iter,PartAnalyzeStep).EQ.0) THEN
    IF (PRESENT(skipVerifyCharge_opt)) THEN
      skipVerifyCharge = skipVerifyCharge_opt
    ELSE
      skipVerifyCharge = .FALSE.
    END IF
    IF(DoVerifyCharge.AND.(.NOT.skipVerifyCharge)) CALL VerifyDepositedCharge()
  END IF
END IF

END SUBROUTINE Deposition


PPURE LOGICAL FUNCTION SFMeasureDistance(v1,v2)
!============================================================================================================================
! Check if the two position vectors coincide in the 1D or 2D projection. If yes, then return .FALSE., else return .TRUE.
! If two points coincide in the direction in which the shape function is not deposited, they are ignored (coincide means that the
! real values are equal up to relative precision of 1e-5)
!============================================================================================================================
USE MOD_PICDepo_Vars ,ONLY: dim_sf,dim_sf_dir,dim_sf_dir1,dim_sf_dir2
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: v1(1:3) !< Input vector 1
REAL, INTENT(IN) :: v2(1:3) !< Input vector 2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SFMeasureDistance = .TRUE. ! Default, also used for dim_sf=3 (3D case)

! Depending on the dimensionality
SELECT CASE (dim_sf)
CASE (1)
  SFMeasureDistance = MERGE(.FALSE. , .TRUE. , ALMOSTEQUALRELATIVE(v1(dim_sf_dir) , v2(dim_sf_dir) , 1e-6))
CASE (2)
  SFMeasureDistance = MERGE(.FALSE. , .TRUE. , ALMOSTEQUALRELATIVE(v1(dim_sf_dir1) , v2(dim_sf_dir1) , 1e-6) .AND. &
                                               ALMOSTEQUALRELATIVE(v1(dim_sf_dir2) , v2(dim_sf_dir2) , 1e-6)       )
END SELECT

END FUNCTION SFMeasureDistance


!===================================================================================================================================
!> Find the corresponding node neighbours across one or more periodic boundaries
!===================================================================================================================================
SUBROUTINE InitializePeriodicNodes(&
#if USE_MPI
  DoNodeMapping,SendNode&
#endif /*USE_MPI*/
)
! MODULES
USE MOD_Globals
USE MOD_Basis                  ,ONLY: BarycentricWeights,InitializeVandermonde
USE MOD_Basis                  ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_Mesh_Vars              ,ONLY: nElems,BoundaryType
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueGlobalNodes, GEO, NodeCoords_Shared, SideInfo_Shared,ElemSideNodeID_Shared
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_PICDepo_Vars
USE MOD_PICDepo_Tools          ,ONLY: CalcCellLocNodeVolumes,ReadTimeAverage
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID, GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: NodeInfo_Shared
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: ELEM_RANK
USE MOD_Particle_Mesh_Vars     ,ONLY: NodeToElemInfo,NodeToElemMapping
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: nProcessors_Global,MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
#if USE_MPI
LOGICAL,INTENT(INOUT) :: DoNodeMapping(0:nProcessors_Global-1)
LOGICAL,INTENT(INOUT) :: SendNode(1:nUniqueGlobalNodes)
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iNode, CNElemID, CNNbElemID, iLocSide, jNode, locElemID, kNode
INTEGER                   :: NbElemID, NbLocSide, PeriodicNode, PVID, SideID, NbSideID, NumPerioNodes, iPeriodNode, zNode,zGlobalNode
REAL                      :: tmpDist, Dist
LOGICAL                   :: NodeDone(4)
INTEGER, ALLOCATABLE      :: PeriodicNodeSourceMap(:,:)
TYPE tPeriodicNodeMap
  INTEGER                       :: nPeriodicNodes
  INTEGER,ALLOCATABLE           :: Mapping(:)
  INTEGER,ALLOCATABLE           :: Rank(:)
END TYPE
TYPE(tPeriodicNodeMap), ALLOCATABLE :: PeriodicNodeMap(:)
INTEGER,ALLOCATABLE       :: PeriodicNodesPerNode(:)
INTEGER                   :: UniqueNodeID
#if USE_MPI
LOGICAL                   :: NoBCSideOnNode,FoundOwnNode
INTEGER                   :: TestNode, minRank, elemCount,jElem,TestElemID
INTEGER, ALLOCATABLE      :: SendPeriodicNodes(:), iSendNode(:), RecvPeriodicNodes(:)
INTEGER                   :: GlobalElemRank, iProc
INTEGER                   :: iRank
! Non-symmetric particle exchange
TYPE(MPI_Request)         :: SendRequestNonSymDepo(0:nProcessors_Global-1)      , RecvRequestNonSymDepo(0:nProcessors_Global-1)

TYPE tPeriodicSendRecv
  INTEGER, ALLOCATABLE    :: Send(:,:)
  INTEGER, ALLOCATABLE    :: Recv(:,:)
  INTEGER, ALLOCATABLE    :: SendNodes(:)
  INTEGER, ALLOCATABLE    :: RecvNodes(:)
END TYPE
TYPE(tPeriodicSendRecv), ALLOCATABLE :: PeriodicSendRecv(:)

TYPE tNodewoBCSide
  INTEGER, ALLOCATABLE    :: RankID(:)
END TYPE
TYPE(tNodewoBCSide), ALLOCATABLE :: NodewoBCSide(:)
#endif
REAL              :: StartT,EndT
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A,I0,A)',ADVANCE='NO') ' | Initializing periodic nodes for cell_volweight_mean...'
GETTIME(StartT)

ALLOCATE(PeriodicNodeSourceMap(1:2*GEO%nPeriodicVectors,1:nUniqueGlobalNodes))
ALLOCATE(PeriodicNodeMap(1:nUniqueGlobalNodes))
PeriodicNodeMap(:)%nPeriodicNodes = 0
PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,:) = 0
PeriodicNodeSourceMap(GEO%nPeriodicVectors+1:2*GEO%nPeriodicVectors,:) = -1
DO locElemID=1, nElems
  DO iLocSide = 1, 6
    SideID=GetGlobalNonUniqueSideID(offsetElem+locElemID,iLocSide)
    CNElemID = GetCNElemID(locElemID+offsetElem)
    IF (SideInfo_Shared(SIDE_BCID,SideID).EQ.0) CYCLE
    IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%PeriodicBC) THEN
      PVID = BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)
      NodeDone = .FALSE.
      NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
      CNNbElemID = GetCNElemID(NbElemID)
      DO NbLocSide = 1, 6
        NbSideID = GetGlobalNonUniqueSideID(NbElemID,NbLocSide)
        IF (SideInfo_Shared(SIDE_BCID,NbSideID).GT.0) THEN
          IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,NbSideID))).EQ.PartBound%PeriodicBC) THEN
            IF (PVID.EQ.-BoundaryType(SideInfo_Shared(SIDE_BCID,NbSideID),BC_ALPHA)) EXIT
          END IF
        END IF
      END DO
      DO iNode=1,4
        IF (PeriodicNodeSourceMap(ABS(PVID),NodeInfo_Shared(ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)).GT.0) CYCLE
        PeriodicNode= 0
        Dist = HUGE(Dist)

        DO jNode = 1,4
          IF (NodeDone(jNode)) CYCLE
          tmpDist = VECNORM3D(NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1) + SIGN( GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID)) &
                -NodeCoords_Shared(1:3,ElemSideNodeID_Shared(jNode,NbLocSide,CNNbElemID)+1))
          IF (tmpDist.LT.Dist) THEN
            PeriodicNode = jNode
            Dist = tmpDist
          END IF
        END DO ! jNode = 1,4
        IF (PeriodicNode.EQ.0) CALL abort(__STAMP__,'Cannot find all periodic nodes for CVWM!')
        NodeDone(PeriodicNode) = .TRUE.
        UniqueNodeID = NodeInfo_Shared(ElemSideNodeID_Shared(PeriodicNode,NbLocSide,CNNbElemID)+1)
        PeriodicNodeSourceMap(ABS(PVID),NodeInfo_Shared(ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)) = UniqueNodeID
#if USE_MPI
        GlobalElemRank = ElemInfo_Shared(ELEM_RANK,NbElemID)
        IF (GlobalElemRank.NE.myRank) THEN
          PeriodicNodeSourceMap(GEO%nPeriodicVectors+ABS(PVID),NodeInfo_Shared(ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)) = GlobalElemRank
        END IF ! GlobalElemRank.NE.myRank
#endif
      END DO ! iNode=1,4
    END IF ! PartBound%PeriodicBC
  END DO ! iLocSide = 1, 6
END DO ! locElemID=1, nElems

#if USE_MPI
! Find periodic nodes which do not have a corresponding BC side to which they are connected
IF (ANY(PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,:).GT.0)) THEN
  ALLOCATE(NodewoBCSide(nUniqueGlobalNodes))
END IF
DO kNode = 1, nUniqueGlobalNodes
  NumPerioNodes = COUNT(PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,kNode).GT.0)
  IF (NumPerioNodes.GT.0) THEN
    minRank = myRank
    ALLOCATE(NodewoBCSide(kNode)%RankID(NodeToElemMapping(2,kNode)))
    NodewoBCSide(kNode)%RankID(:) = -1
    elemCount = 0
    DO jElem = NodeToElemMapping(1,kNode) + 1, NodeToElemMapping(1,kNode) + NodeToElemMapping(2,kNode)
      elemCount = elemCount  + 1
      TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
      IF (ElemInfo_Shared(ELEM_RANK,TestElemID).EQ.myRank) CYCLE
      NoBCSideOnNode = .TRUE.
      LocSideLoop: DO iLocSide = 1, 6
        SideID=GetGlobalNonUniqueSideID(TestElemID,iLocSide)
        CNElemID = GetCNElemID(TestElemID)
        FoundOwnNode = .FALSE.
        IF (SideInfo_Shared(SIDE_BCID,SideID).EQ.0) CYCLE LocSideLoop
        IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%PeriodicBC) THEN
          NoBCSideOnNode = .FALSE.
          PVID = BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)
          IF (PeriodicNodeSourceMap(ABS(PVID),kNode).GT.0) CYCLE LocSideLoop
          NodeCycle: DO iNode=1,4
            IF (NodeInfo_Shared(ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1).EQ.kNode) THEN
              FoundOwnNode=.TRUE.
              EXIT NodeCycle
            END IF
          END DO NodeCycle
          IF (.NOT.FoundOwnNode) CYCLE LocSideLoop
          NodeDone = .FALSE.
          NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
          CNNbElemID = GetCNElemID(NbElemID)
          DO NbLocSide = 1, 6
            NbSideID = GetGlobalNonUniqueSideID(NbElemID,NbLocSide)
            IF (SideInfo_Shared(SIDE_BCID,NbSideID).GT.0) THEN
              IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,NbSideID))).EQ.PartBound%PeriodicBC) THEN
                IF (PVID.EQ.-BoundaryType(SideInfo_Shared(SIDE_BCID,NbSideID),BC_ALPHA)) EXIT
              END IF
            END IF
          END DO
          PeriodicNode= 0
          Dist = HUGE(Dist)

          DO jNode = 1,4
            IF (NodeDone(jNode)) CYCLE
            tmpDist = VECNORM3D(NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1) + SIGN( GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID)) &
                  -NodeCoords_Shared(1:3,ElemSideNodeID_Shared(jNode,NbLocSide,CNNbElemID)+1))
            IF (tmpDist.LT.Dist) THEN
              PeriodicNode = jNode
              Dist = tmpDist
            END IF
          END DO ! jNode = 1,4
          IF (PeriodicNode.EQ.0) CALL abort(__STAMP__,'Cannot find all periodic nodes for CVWM!')
          NodeDone(PeriodicNode) = .TRUE.
          UniqueNodeID = NodeInfo_Shared(ElemSideNodeID_Shared(PeriodicNode,NbLocSide,CNNbElemID)+1)
          PeriodicNodeSourceMap(ABS(PVID),NodeInfo_Shared(ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)) = UniqueNodeID
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,NbElemID)
          IF (GlobalElemRank.NE.myRank) THEN
            PeriodicNodeSourceMap(GEO%nPeriodicVectors+ABS(PVID),NodeInfo_Shared(ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)) = GlobalElemRank
          END IF ! GlobalElemRank.NE.myRank
        END IF ! PartBound%PeriodicBC
      END DO LocSideLoop ! iLocSide = 1, 6
      IF (NoBCSideOnNode) THEN
        NodewoBCSide(kNode)%RankID(elemCount) = ElemInfo_Shared(ELEM_RANK,TestElemID)
      ELSE
        minRank = MIN(minRank,ElemInfo_Shared(ELEM_RANK,TestElemID))
      END IF
    END DO
    IF (minRank.NE.myRank) NodewoBCSide(kNode)%RankID(:) = -1
  END IF
END DO



ALLOCATE(SendPeriodicNodes(0:nProcessors_Global-1), PeriodicSendRecv(0:nProcessors_Global-1))
ALLOCATE(iSendNode(0:nProcessors_Global-1),RecvPeriodicNodes(0:nProcessors_Global-1))

iSendNode = 0
SendPeriodicNodes = 0; RecvPeriodicNodes =0
IF (ANY(PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,:).GT.0)) THEN
  DO iNode = 1, nUniqueGlobalNodes
    IF (ALLOCATED(NodewoBCSide(iNode)%RankID)) THEN
      NumPerioNodes = COUNT(NodewoBCSide(iNode)%RankID(:) .NE.-1)
      IF (NumPerioNodes.GT.0) THEN
        DO jElem = 1,  NodeToElemMapping(2,iNode)
          IF (NodewoBCSide(iNode)%RankID(jElem).NE.-1) THEN
            SendPeriodicNodes(NodewoBCSide(iNode)%RankID(jElem)) = SendPeriodicNodes(NodewoBCSide(iNode)%RankID(jElem)) + 1
          END IF
        END DO
      END IF
    END IF
  END DO

  DO iRank= 0, nProcessors_Global-1
    IF (iRank.EQ.myRank) CYCLE
    IF (SendPeriodicNodes(iRank).GT.0) THEN
      ALLOCATE(PeriodicSendRecv(iRank)%Send(2*GEO%nPeriodicVectors+1,SendPeriodicNodes(iRank)))
      PeriodicSendRecv(iRank)%Send = 0
    END IF
  END DO
  DO iNode = 1, nUniqueGlobalNodes
   IF (ALLOCATED(NodewoBCSide(iNode)%RankID)) THEN
      NumPerioNodes = COUNT(NodewoBCSide(iNode)%RankID(:) .NE.-1)
      IF (NumPerioNodes.GT.0) THEN
        DO jElem = 1,  NodeToElemMapping(2,iNode)
          iRank = NodewoBCSide(iNode)%RankID(jElem)
          IF (iRank.NE.-1) THEN
            iSendNode(iRank) = iSendNode(iRank) + 1
            PeriodicSendRecv(iRank)%Send(1:2*GEO%nPeriodicVectors,iSendNode(iRank)) = PeriodicNodeSourceMap(1:2*GEO%nPeriodicVectors, iNode)
            PeriodicSendRecv(iRank)%Send(2*GEO%nPeriodicVectors+1,iSendNode(iRank)) = iNode
          END IF
        END DO
      END IF
    END IF
  END DO
END IF

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_IRECV( RecvPeriodicNodes(iProc)     &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1667                         &
                , MPI_COMM_PICLAS              &
                , RecvRequestNonSymDepo(iProc) &
                , IERROR)
  CALL MPI_ISEND( SendPeriodicNodes(iProc)     &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1667                         &
                , MPI_COMM_PICLAS              &
                , SendRequestNonSymDepo(iProc) &
                , IERROR)
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  IF (RecvPeriodicNodes(iProc).NE.0) THEN
    ALLOCATE(PeriodicSendRecv(iProc)%Recv(2*GEO%nPeriodicVectors+1,RecvPeriodicNodes(iProc)))
    CALL MPI_IRECV( PeriodicSendRecv(iProc)%Recv(:,:)                   &
                  , RecvPeriodicNodes(iProc)*(2*GEO%nPeriodicVectors+1) &
                  , MPI_INTEGER                                         &
                  , iProc                                               &
                  , 667                                                 &
                  , MPI_COMM_PICLAS                                     &
                  , RecvRequestNonSymDepo(iProc)                        &
                  , IERROR)
  END IF
  IF (SendPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_ISEND( PeriodicSendRecv(iProc)%Send                        &
                  , SendPeriodicNodes(iProc)*(2*GEO%nPeriodicVectors+1) &
                  , MPI_INTEGER                                         &
                  , iProc                                               &
                  , 667                                                 &
                  , MPI_COMM_PICLAS                                     &
                  , SendRequestNonSymDepo(iProc)                        &
                  , IERROR)
  END IF
END DO
! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  IF (RecvPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
  IF (SendPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (RecvPeriodicNodes(iRank).GT.0) THEN
    DO iNode = 1, RecvPeriodicNodes(iRank)
      zGlobalNode = PeriodicSendRecv(iRank)%Recv(2*GEO%nPeriodicVectors+1,iNode)
      DO jNode = 1, GEO%nPeriodicVectors
        IF((PeriodicNodeSourceMap(jNode,zGlobalNode).EQ.0).AND.(PeriodicSendRecv(iRank)%Recv(jNode,iNode).NE.0))THEN
          PeriodicNodeSourceMap(jNode,zGlobalNode) = PeriodicSendRecv(iRank)%Recv(jNode,iNode)
          PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,zGlobalNode) = PeriodicSendRecv(iRank)%Recv(jNode+GEO%nPeriodicVectors,iNode)
        ELSEIF((PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,zGlobalNode).LT.0).AND.(PeriodicSendRecv(iRank)%Recv(jNode+GEO%nPeriodicVectors,iNode).GE.0))THEN
          PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,zGlobalNode) = PeriodicSendRecv(iRank)%Recv(jNode+GEO%nPeriodicVectors,iNode)
        END IF ! (PeriodicNodeSourceMap(jNode,zGlobalNode).EQ.0)
      END DO ! jNode = 1, GEO%nPeriodicVectors
    END DO
    DEALLOCATE(PeriodicSendRecv(iRank)%Recv)
  END IF
  IF (SendPeriodicNodes(iRank).GT.0) THEN
    DEALLOCATE(PeriodicSendRecv(iRank)%Send)
  END IF
END DO

IF (ALLOCATED(NodewoBCSide)) THEN
  DO iNode = 1, nUniqueGlobalNodes
    SDEALLOCATE(NodewoBCSide(iNode)%RankID)
  END DO
  DEALLOCATE(NodewoBCSide)
END IF

iSendNode = 0
SendPeriodicNodes = 0; RecvPeriodicNodes =0
DO iNode = 1, nUniqueGlobalNodes
  NumPerioNodes = COUNT(PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,iNode).GT.0)
  IF (NumPerioNodes.GT.0) THEN
    DO jNode= 1, GEO%nPeriodicVectors
      IF (PeriodicNodeSourceMap(jNode,iNode).NE.0) THEN
        IF(PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode).NE.myrank)THEN
          IF (PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode).NE.-1) THEN
            SendPeriodicNodes(PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode)) = &
              SendPeriodicNodes(PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode)) + 1
          END IF
        END IF ! PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode).NE.myrank
      END IF
    END DO
  END IF
END DO

DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (SendPeriodicNodes(iRank).GT.0) THEN
    ALLOCATE(PeriodicSendRecv(iRank)%RecvNodes(SendPeriodicNodes(iRank)))

    PeriodicSendRecv(iRank)%RecvNodes = 0
  END IF
END DO
DO iNode = 1, nUniqueGlobalNodes
  NumPerioNodes = COUNT(PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,iNode).GT.0)
  IF (NumPerioNodes.GT.0) THEN
    DO jNode= 1, GEO%nPeriodicVectors
      IF (PeriodicNodeSourceMap(jNode,iNode).NE.0) THEN
        iRank = PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode)
        IF (iRank.EQ.myRank) CYCLE
        IF (iRank.NE.-1) THEN
          iSendNode(iRank) = iSendNode(iRank) + 1
          PeriodicSendRecv(iRank)%RecvNodes(iSendNode(iRank)) = PeriodicNodeSourceMap(jNode,iNode)
        END IF
      END IF
    END DO
  END IF
END DO

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_IRECV( RecvPeriodicNodes(iProc)     &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1667                         &
                , MPI_COMM_PICLAS              &
                , RecvRequestNonSymDepo(iProc) &
                , IERROR)
  CALL MPI_ISEND( SendPeriodicNodes(iProc)     &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1667                         &
                , MPI_COMM_PICLAS              &
                , SendRequestNonSymDepo(iProc) &
                , IERROR)
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  IF (RecvPeriodicNodes(iProc).NE.0) THEN
    ALLOCATE(PeriodicSendRecv(iProc)%SendNodes(RecvPeriodicNodes(iProc)))
    CALL MPI_IRECV( PeriodicSendRecv(iProc)%SendNodes(:) &
                  , RecvPeriodicNodes(iProc)             &
                  , MPI_INTEGER                          &
                  , iProc                                &
                  , 667                                  &
                  , MPI_COMM_PICLAS                      &
                  , RecvRequestNonSymDepo(iProc)         &
                  , IERROR)
  END IF
  IF (SendPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_ISEND( PeriodicSendRecv(iProc)%RecvNodes(:) &
                  , SendPeriodicNodes(iProc)             &
                  , MPI_INTEGER                          &
                  , iProc                                &
                  , 667                                  &
                  , MPI_COMM_PICLAS                      &
                  , SendRequestNonSymDepo(iProc)         &
                  , IERROR)
  END IF
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  IF (RecvPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
  IF (SendPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

iSendNode = 0
DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (RecvPeriodicNodes(iRank).GT.0) THEN
    ALLOCATE(PeriodicSendRecv(iRank)%Send(GEO%nPeriodicVectors+1,RecvPeriodicNodes(iRank)))
    DO iNode = 1, RecvPeriodicNodes(iRank)
      zGlobalNode = PeriodicSendRecv(iRank)%SendNodes(iNode)
      iSendNode(iRank) = iSendNode(iRank) + 1
      PeriodicSendRecv(iRank)%Send(1:GEO%nPeriodicVectors,iSendNode(iRank)) = PeriodicNodeSourceMap(1:GEO%nPeriodicVectors, zGlobalNode)
      PeriodicSendRecv(iRank)%Send(GEO%nPeriodicVectors+1,iSendNode(iRank)) = zGlobalNode
    END DO
  END IF
END DO

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  IF (SendPeriodicNodes(iProc).NE.0) THEN
    ALLOCATE(PeriodicSendRecv(iProc)%Recv(GEO%nPeriodicVectors+1,SendPeriodicNodes(iProc)))
    CALL MPI_IRECV( PeriodicSendRecv(iProc)%Recv(:,:)                 &
                  , SendPeriodicNodes(iProc)*(GEO%nPeriodicVectors+1) &
                  , MPI_INTEGER                                       &
                  , iProc                                             &
                  , 667                                               &
                  , MPI_COMM_PICLAS                                   &
                  , RecvRequestNonSymDepo(iProc)                      &
                  , IERROR)
  END IF
  IF (RecvPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_ISEND( PeriodicSendRecv(iProc)%Send                      &
                  , RecvPeriodicNodes(iProc)*(GEO%nPeriodicVectors+1) &
                  , MPI_INTEGER                                       &
                  , iProc                                             &
                  , 667                                               &
                  , MPI_COMM_PICLAS                                   &
                  , SendRequestNonSymDepo(iProc)                      &
                  , IERROR)
  END IF
END DO
! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  IF (SendPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
  IF (RecvPeriodicNodes(iProc).NE.0) THEN
    CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (SendPeriodicNodes(iRank).GT.0) THEN
    DO iNode = 1, SendPeriodicNodes(iRank)
      zGlobalNode = PeriodicSendRecv(iRank)%Recv(GEO%nPeriodicVectors+1,iNode)
      DO jNode = 1, GEO%nPeriodicVectors
        IF((PeriodicNodeSourceMap(jNode,zGlobalNode).EQ.0).AND.(PeriodicSendRecv(iRank)%Recv(jNode,iNode).NE.0))THEN
          PeriodicNodeSourceMap(jNode,zGlobalNode) = PeriodicSendRecv(iRank)%Recv(jNode,iNode)
        END IF ! (PeriodicNodeSourceMap(jNode,zGlobalNode).EQ.0)
      END DO ! jNode = 1, GEO%nPeriodicVectors
    END DO
    DEALLOCATE(PeriodicSendRecv(iRank)%Recv)
    DEALLOCATE(PeriodicSendRecv(iRank)%RecvNodes)
  END IF
  IF (RecvPeriodicNodes(iRank).GT.0) THEN
    DEALLOCATE(PeriodicSendRecv(iRank)%Send)
    DEALLOCATE(PeriodicSendRecv(iRank)%SendNodes)
  END IF
END DO
#endif
DO iNode = 1, nUniqueGlobalNodes
  NumPerioNodes = COUNT(PeriodicNodeSourceMap(1:GEO%nPeriodicVectors,iNode).GT.0)
  IF (NumPerioNodes.GT.1) NumPerioNodes = 2**NumPerioNodes - 1
  IF (NumPerioNodes.NE.0) THEN
    PeriodicNodeMap(iNode)%nPeriodicNodes = NumPerioNodes
    ALLOCATE(PeriodicNodeMap(iNode)%Mapping(NumPerioNodes), PeriodicNodeMap(iNode)%Rank(NumPerioNodes))
    PeriodicNodeMap(iNode)%Mapping = 0
    PeriodicNodeMap(iNode)%Rank = -1
    iPeriodNode = 0
    DO jNode = 1, GEO%nPeriodicVectors
      IF (PeriodicNodeSourceMap(jNode,iNode).NE.0) THEN
        iPeriodNode = iPeriodNode + 1
        PeriodicNodeMap(iNode)%Mapping(iPeriodNode) = PeriodicNodeSourceMap(jNode,iNode)
        PeriodicNodeMap(iNode)%Rank(iPeriodNode) = PeriodicNodeSourceMap(jNode+GEO%nPeriodicVectors,iNode)
      END IF
    END DO
    IF (NumPerioNodes.GT.0) THEN
      DO jNode = 1, GEO%nPeriodicVectors
        IF (PeriodicNodeSourceMap(jNode,iNode).NE.0) THEN
          DO zNode = 1, GEO%nPeriodicVectors
            zGlobalNode = PeriodicNodeSourceMap(zNode,PeriodicNodeSourceMap(jNode,iNode))
            IF ((zGlobalNode.NE.0).AND.(zGlobalNode.NE.iNode)) THEN
              IF (.NOT.ANY(PeriodicNodeMap(iNode)%Mapping(:).EQ.zGlobalNode)) THEN
                iPeriodNode = iPeriodNode + 1
                PeriodicNodeMap(iNode)%Mapping(iPeriodNode) = zGlobalNode
                PeriodicNodeMap(iNode)%Rank(iPeriodNode) = PeriodicNodeSourceMap(zNode+GEO%nPeriodicVectors,PeriodicNodeSourceMap(jNode,iNode))
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF
END DO
#if USE_MPI
IF (GEO%nPeriodicVectors.GT.1) THEN
  iSendNode = 0
  SendPeriodicNodes = 0; RecvPeriodicNodes =0
  DO iNode = 1, nUniqueGlobalNodes
    IF (PeriodicNodeMap(iNode)%nPeriodicNodes.GT.0) THEN
      IF (ANY(PeriodicNodeMap(iNode)%Mapping.EQ.0)) THEN
        DO jNode = 1, PeriodicNodeMap(iNode)%nPeriodicNodes
          iRank = PeriodicNodeMap(iNode)%Rank(jNode)
          IF (iRank.EQ.myRank) CYCLE
          IF (iRank.NE.-1) THEN
            SendPeriodicNodes(iRank) = SendPeriodicNodes(iRank) + 1
          END IF
        END DO
      END IF
    END IF
  END DO
  DO iRank= 0, nProcessors_Global-1
    IF (iRank.EQ.myRank) CYCLE
    IF (SendPeriodicNodes(iRank).GT.0) THEN
      ALLOCATE(PeriodicSendRecv(iRank)%RecvNodes(SendPeriodicNodes(iRank)))
      PeriodicSendRecv(iRank)%RecvNodes = 0
    END IF
  END DO
  DO iNode = 1, nUniqueGlobalNodes
    IF (PeriodicNodeMap(iNode)%nPeriodicNodes.GT.0) THEN
      IF (ANY(PeriodicNodeMap(iNode)%Mapping.EQ.0)) THEN
        DO jNode = 1, PeriodicNodeMap(iNode)%nPeriodicNodes
          iRank = PeriodicNodeMap(iNode)%Rank(jNode)
          IF (iRank.EQ.myRank) CYCLE
          IF (iRank.NE.-1) THEN
            iSendNode(iRank) = iSendNode(iRank) + 1
            PeriodicSendRecv(iRank)%RecvNodes(iSendNode(iRank)) = PeriodicNodeMap(iNode)%Mapping(jNode)
          END IF
        END DO
      END IF
    END IF
  END DO

  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    CALL MPI_IRECV( RecvPeriodicNodes(iProc)     &
                  , 1                            &
                  , MPI_INTEGER                  &
                  , iProc                        &
                  , 1667                         &
                  , MPI_COMM_PICLAS              &
                  , RecvRequestNonSymDepo(iProc) &
                  , IERROR)
    CALL MPI_ISEND( SendPeriodicNodes(iProc)     &
                  , 1                            &
                  , MPI_INTEGER                  &
                  , iProc                        &
                  , 1667                         &
                  , MPI_COMM_PICLAS              &
                  , SendRequestNonSymDepo(iProc) &
                  , IERROR)
  END DO

 ! Finish communication
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO

  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    IF (RecvPeriodicNodes(iProc).NE.0) THEN
      ALLOCATE(PeriodicSendRecv(iProc)%SendNodes(RecvPeriodicNodes(iProc)))
      CALL MPI_IRECV( PeriodicSendRecv(iProc)%SendNodes(:) &
                    , RecvPeriodicNodes(iProc)             &
                    , MPI_INTEGER                          &
                    , iProc                                &
                    , 667                                  &
                    , MPI_COMM_PICLAS                      &
                    , RecvRequestNonSymDepo(iProc)         &
                    , IERROR)
    END IF
    IF (SendPeriodicNodes(iProc).NE.0) THEN
      CALL MPI_ISEND( PeriodicSendRecv(iProc)%RecvNodes(:) &
                    , SendPeriodicNodes(iProc)             &
                    , MPI_INTEGER                          &
                    , iProc                                &
                    , 667                                  &
                    , MPI_COMM_PICLAS                      &
                    , SendRequestNonSymDepo(iProc)         &
                    , IERROR)
    END IF
  END DO

  ! Finish communication
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    IF (RecvPeriodicNodes(iProc).NE.0) THEN
      CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END IF
    IF (SendPeriodicNodes(iProc).NE.0) THEN
      CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END IF
  END DO

  iSendNode = 0
  DO iRank= 0, nProcessors_Global-1
    IF (iRank.EQ.myRank) CYCLE
    IF (RecvPeriodicNodes(iRank).GT.0) THEN
      ALLOCATE(PeriodicSendRecv(iRank)%Send(2**GEO%nPeriodicVectors,RecvPeriodicNodes(iRank)))
      PeriodicSendRecv(iRank)%Send = 0
      DO iNode = 1, RecvPeriodicNodes(iRank)
        zGlobalNode = PeriodicSendRecv(iRank)%SendNodes(iNode)
        iSendNode(iRank) = iSendNode(iRank) + 1
        PeriodicSendRecv(iRank)%Send(1:PeriodicNodeMap(zGlobalNode)%nPeriodicNodes,iSendNode(iRank)) &
          = PeriodicNodeMap(zGlobalNode)%Mapping(1:PeriodicNodeMap(zGlobalNode)%nPeriodicNodes)
        PeriodicSendRecv(iRank)%Send(2**GEO%nPeriodicVectors,iSendNode(iRank)) = zGlobalNode
      END DO
    END IF
  END DO

  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    IF (SendPeriodicNodes(iProc).NE.0) THEN
      ALLOCATE(PeriodicSendRecv(iProc)%Recv(2**GEO%nPeriodicVectors,SendPeriodicNodes(iProc)))
      CALL MPI_IRECV( PeriodicSendRecv(iProc)%Recv(:,:)                  &
                    , SendPeriodicNodes(iProc)*(2**GEO%nPeriodicVectors) &
                    , MPI_INTEGER                                        &
                    , iProc                                              &
                    , 667                                                &
                    , MPI_COMM_PICLAS                                    &
                    , RecvRequestNonSymDepo(iProc)                       &
                    , IERROR)
    END IF
    IF (RecvPeriodicNodes(iProc).NE.0) THEN
      CALL MPI_ISEND( PeriodicSendRecv(iProc)%Send                       &
                    , RecvPeriodicNodes(iProc)*(2**GEO%nPeriodicVectors) &
                    , MPI_INTEGER                                        &
                    , iProc                                              &
                    , 667                                                &
                    , MPI_COMM_PICLAS                                    &
                    , SendRequestNonSymDepo(iProc)                       &
                    , IERROR)
    END IF
  END DO
  ! Finish communication
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE
    IF (SendPeriodicNodes(iProc).NE.0) THEN
      CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END IF
    IF (RecvPeriodicNodes(iProc).NE.0) THEN
      CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END IF
  END DO

  DO iRank= 0, nProcessors_Global-1
    IF (iRank.EQ.myRank) CYCLE
    IF (SendPeriodicNodes(iRank).GT.0) THEN
      DO iNode = 1, SendPeriodicNodes(iRank)
        zGlobalNode = PeriodicSendRecv(iRank)%Recv(2**GEO%nPeriodicVectors,iNode)
        PeriodicNodeMap(zGlobalNode)%Mapping(1:PeriodicNodeMap(zGlobalNode)%nPeriodicNodes) &
        = PeriodicSendRecv(iRank)%Recv(1:PeriodicNodeMap(zGlobalNode)%nPeriodicNodes,iNode)
      END DO
      DEALLOCATE(PeriodicSendRecv(iRank)%Recv)
      DEALLOCATE(PeriodicSendRecv(iRank)%RecvNodes)
    END IF
    IF (RecvPeriodicNodes(iRank).GT.0) THEN
      DEALLOCATE(PeriodicSendRecv(iRank)%Send)
      DEALLOCATE(PeriodicSendRecv(iRank)%SendNodes)
    END IF
  END DO

  DEALLOCATE(PeriodicSendRecv, iSendNode, SendPeriodicNodes, RecvPeriodicNodes)
END IF
#endif
DO iNode = 1, nUniqueGlobalNodes
  IF (PeriodicNodeMap(iNode)%nPeriodicNodes.GT.0) THEN
    IF (ANY(PeriodicNodeMap(iNode)%Mapping.EQ.0)) THEN
      DO jNode = 1, PeriodicNodeMap(iNode)%nPeriodicNodes
        IF (PeriodicNodeMap(iNode)%Mapping(jNode).EQ.0) THEN
          DO kNode =1, jNode - 1
            zGlobalNode = PeriodicNodeMap(iNode)%Mapping(kNode)
            DO zNode = 1, PeriodicNodeMap(zGlobalNode)%nPeriodicNodes
              IF ((PeriodicNodeMap(zGlobalNode)%Mapping(zNode).NE.0).AND.(PeriodicNodeMap(zGlobalNode)%Mapping(zNode).NE.iNode)) THEN
                IF (.NOT.ANY(PeriodicNodeMap(iNode)%Mapping(:).EQ.PeriodicNodeMap(zGlobalNode)%Mapping(zNode))) THEN
                  PeriodicNodeMap(iNode)%Mapping(jNode) = PeriodicNodeMap(zGlobalNode)%Mapping(zNode)
                END IF
              END IF
            END DO
          END DO
        END IF
      END DO
    END IF
  END IF
END DO
#if USE_MPI
DO iNode = 1, nUniqueGlobalNodes
  IF (PeriodicNodeMap(iNode)%nPeriodicNodes.GT.0) THEN
    DO jNode = 1, PeriodicNodeMap(iNode)%nPeriodicNodes
      TestNode = PeriodicNodeMap(iNode)%Mapping(jNode)
      DO jElem = NodeToElemMapping(1,TestNode) + 1, NodeToElemMapping(1,TestNode) + NodeToElemMapping(2,TestNode)
        TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
        GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
        IF (GlobalElemRank.NE.myRank) THEN
          SendNode(TestNode) = .TRUE.
          DoNodeMapping(GlobalElemRank) = .TRUE.
        END IF
      END DO
    END DO
  END IF
END DO
#endif /*USE_MPI*/
DEALLOCATE(PeriodicNodeSourceMap)

! FERTIG
#if USE_MPI
CALL Allocate_Shared((/nUniqueGlobalNodes/),Periodic_nNodes_Shared_Win    ,Periodic_nNodes_Shared)
CALL Allocate_Shared((/nUniqueGlobalNodes/),Periodic_offsetNode_Shared_Win,Periodic_offsetNode_Shared)
CALL MPI_WIN_LOCK_ALL(0,Periodic_nNodes_Shared_Win    ,IERROR)
CALL MPI_WIN_LOCK_ALL(0,Periodic_offsetNode_Shared_Win,IERROR)
Periodic_nNodes => Periodic_nNodes_Shared
Periodic_offsetNode => Periodic_offsetNode_Shared
IF (myComputeNodeRank.EQ.0) THEN
  Periodic_nNodes = 0
  Periodic_offsetNode = 0
END IF ! myComputeNodeRank.EQ.0
CALL BARRIER_AND_SYNC(Periodic_nNodes_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Periodic_offsetNode_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(Periodic_nNodes(1:nUniqueGlobalNodes))
Periodic_nNodes = 0
ALLOCATE(Periodic_offsetNode(1:nUniqueGlobalNodes))
Periodic_offsetNode = 0
#endif /*USE_MPI*/

ALLOCATE(PeriodicNodesPerNode(nUniqueGlobalNodes))
DO iNode = 1,nUniqueGlobalNodes
  PeriodicNodesPerNode(iNode) = PeriodicNodeMap(iNode)%nPeriodicNodes
END DO ! iNode = nUniqueGlobalNodes

#if USE_MPI
IF(myComputeNodeRank.EQ.0)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE        ,PeriodicNodesPerNode,nUniqueGlobalNodes,MPI_INTEGER,MPI_MAX,0,MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/
  nTotalPeriodicNodes = 0
  DO iNode = 1,nUniqueGlobalNodes
    Periodic_offsetNode(iNode) = nTotalPeriodicNodes
    Periodic_nNodes(    iNode) = PeriodicNodesPerNode(iNode)
    nTotalPeriodicNodes        = nTotalPeriodicNodes + PeriodicNodesPerNode(iNode)
  END DO ! iNode = nUniqueGlobalNodes
#if USE_MPI
ELSE
  CALL MPI_REDUCE(PeriodicNodesPerNode,0                   ,nUniqueGlobalNodes,MPI_INTEGER,MPI_MAX,0,MPI_COMM_SHARED,IERROR)
END IF
! Root knows the global number, now broadcast to other procs
CALL MPI_BCAST(nTotalPeriodicNodes,1,MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

CALL BARRIER_AND_SYNC(Periodic_nNodes_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Periodic_offsetNode_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

IF(nTotalPeriodicNodes.GT.0) THEN
#if USE_MPI
  CALL Allocate_Shared((/nTotalPeriodicNodes/),Periodic_Nodes_Shared_Win,Periodic_Nodes_Shared)
  CALL MPI_WIN_LOCK_ALL(0,Periodic_Nodes_Shared_Win     ,IERROR)
  Periodic_Nodes => Periodic_Nodes_Shared
  IF (myComputeNodeRank.EQ.0) Periodic_Nodes = 0
  CALL BARRIER_AND_SYNC(Periodic_Nodes_Shared_Win,MPI_COMM_SHARED)
#else
  ALLOCATE(Periodic_Nodes(1:nTotalPeriodicNodes))
  Periodic_Nodes = 0
#endif /*USE_MPI*/

  ! Every processor loops over its own periodic map and fills the Periodic_Nodes_Shared_Win array. MPI_ACCUMULATE ensures that data
  ! is consistent
  DO iNode = 1,nUniqueGlobalNodes
    ASSOCIATE(offset => Periodic_offsetNode(iNode))

      IF (PeriodicNodeMap(iNode)%nPeriodicNodes.GT.0) THEN
#if USE_MPI
        CALL MPI_ACCUMULATE(PeriodicNodeMap(iNode)%Mapping             ,PeriodicNodeMap(iNode)%nPeriodicNodes, MPI_INTEGER, &
                            0    ,INT(offset*SIZE_INT,MPI_ADDRESS_KIND),PeriodicNodeMap(iNode)%nPeriodicNodes, MPI_INTEGER, &
                            MPI_REPLACE                                ,Periodic_Nodes_Shared_Win            , iError)
#else
        Periodic_Nodes(1+offset:offset+PeriodicNodeMap(iNode)%nPeriodicNodes) = PeriodicNodeMap(iNode)%Mapping
#endif /*USE_MPI*/
    END IF ! PeriodicNodeMap(iNode)%nPeriodicNodes.GT.0

    END ASSOCIATE
  END DO ! iNode = nUniqueGlobalNodes
#if USE_MPI
  CALL MPI_WIN_FLUSH(0 ,Periodic_Nodes_Shared_Win,iError)
  CALL BARRIER_AND_SYNC(Periodic_Nodes_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/
END IF

SDEALLOCATE(PeriodicNodesPerNode)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!',DisplayLine=.FALSE.)

END SUBROUTINE InitializePeriodicNodes


SUBROUTINE FinalizeDeposition()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic deposition
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,PeriodicSFCaseMatrix
USE MOD_PICDepo_Vars
#if USE_MPI
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
#if USE_LOADBALANCE
USE MOD_PICDepo_MPI            ,ONLY: ExchangeNodeSourceExtMPI,ExchangeSurfNodeSourceMPI
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
!USE MOD_Particle_Mesh_Vars    ,ONLY: GlobalElem2CNTotalElem,GlobalElem2CNTotalElem_Shared!,GlobalElem2CNTotalElem_Shared_Win
!USE MOD_MPI_Shared_Vars       ,ONLY: nComputeNodeProcessors,nProcessors_Global
USE MOD_LoadBalance_Vars       ,ONLY: NodeSourceExtEquiLB!,PartSourceLB
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: NodeInfo_Shared,ElemNodeID_Shared
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Boundary_Vars ,ONLY: Do2DSurfaceCharge
USE MOD_Mesh_Vars              ,ONLY: NonUniqueGlobalNodeIDToFEMVertexID,NonUniqueGlobalSideIDToNonUniqueGlobalNodeID
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
INTEGER,PARAMETER :: N_variables=1
INTEGER           :: iElem,CNElemID
INTEGER           :: NodeID(1:8)
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
SDEALLOCATE(GaussBorder)
SDEALLOCATE(Vdm_EquiN_GaussN)
SDEALLOCATE(Knots)
SDEALLOCATE(GaussBGMIndex)
SDEALLOCATE(GaussBGMFactor)
SDEALLOCATE(GEO%PeriodicBGMVectors)
SDEALLOCATE(BGMSource)
SDEALLOCATE(GPWeight)
SDEALLOCATE(ElemRadius2_sf)
SDEALLOCATE(Vdm_NDepo_GaussN)
SDEALLOCATE(DDMassInv)
SDEALLOCATE(XiNDepo)
SDEALLOCATE(swGPNDepo)
SDEALLOCATE(wBaryNDepo)
SDEALLOCATE(NDepochooseK)
SDEALLOCATE(tempcharge)
SDEALLOCATE(CellVolWeight)
SDEALLOCATE(CellVolWeight_Volumes)
SDEALLOCATE(ChargeSFDone)
SDEALLOCATE(PeriodicSFCaseMatrix)
SDEALLOCATE(N_ShapeTmp)

#if USE_MPI
SDEALLOCATE(FlagShapeElem)
SDEALLOCATE(SendDofShapeID)
SDEALLOCATE(CNRankToSendRank)

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

IF(DoDeposition)THEN
  ! Deposition-dependent arrays
  SELECT CASE(TRIM(DepositionType))
  CASE('cell_volweight_mean')
    CALL UNLOCK_AND_FREE(NodeVolume_Shared_Win)
    IF(GEO%nPeriodicVectors.GT.0)THEN
      CALL UNLOCK_AND_FREE(Periodic_nNodes_Shared_Win)
      CALL UNLOCK_AND_FREE(Periodic_offsetNode_Shared_Win)
      IF(nTotalPeriodicNodes.GT.0) CALL UNLOCK_AND_FREE(Periodic_Nodes_Shared_Win)
    END IF ! GEO%nPeriodicVectors.GT.0
  CASE('shape_function_adaptive')
    CALL UNLOCK_AND_FREE(SFElemr2_Shared_Win)
  END SELECT

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

  ADEALLOCATE(NodeVolume_Shared)
  ADEALLOCATE(Periodic_Nodes_Shared)
  ADEALLOCATE(Periodic_nNodes_Shared)
  ADEALLOCATE(Periodic_offsetNode_Shared)
END IF ! DoDeposition

! Then, free the pointers or arrays
#endif /*USE_MPI*/

IF(DoDeposition)THEN
  ADEALLOCATE(NodeVolume)
  ADEALLOCATE(Periodic_Nodes)
  ADEALLOCATE(Periodic_nNodes)
  ADEALLOCATE(Periodic_offsetNode)

  ! Deposition-dependent pointers/arrays
  SELECT CASE(TRIM(DepositionType))
    CASE('cell_volweight_mean')
    CASE('shape_function_adaptive')
      ADEALLOCATE(SFElemr2_Shared)
  END SELECT
END IF

#if USE_LOADBALANCE
IF ((PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN

  !IF(DoDeposition)THEN
  !  SDEALLOCATE(PartSourceLB)
  !  ALLOCATE(PartSourceLB(1:4,0:PP_N,0:PP_N,0:PP_N,nElems))
  !  CALL abort(__STAMP__,'not implemented')
  !  !PartSourceLB = PartSource
  !END IF ! DoDeposition

  IF(DoDielectricSurfaceCharge)THEN
    IF(DoDeposition) CALL ExchangeNodeSourceExtMPI()
    !SDEALLOCATE(NodeSourceExtEquiLB)
    !ALLOCATE(NodeSourceExtEquiLB(1:4,0:PP_N,0:PP_N,0:PP_N,nElems))
    ALLOCATE(NodeSourceExtEquiLB(1:N_variables,0:1,0:1,0:1,nElems))
    ! Loop over all elements and store absolute charge values in equidistantly distributed nodes of PP_N=1
    DO iElem=1,PP_nElems
      ! Copy values to equidistant distribution
      CNElemID = GetCNElemID(iElem+offsetElem)
      NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,CNElemID))
      NodeSourceExtEquiLB(1,0,0,0,iElem) = NodeSourceExt(NodeID(1))
      NodeSourceExtEquiLB(1,1,0,0,iElem) = NodeSourceExt(NodeID(2))
      NodeSourceExtEquiLB(1,1,1,0,iElem) = NodeSourceExt(NodeID(3))
      NodeSourceExtEquiLB(1,0,1,0,iElem) = NodeSourceExt(NodeID(4))
      NodeSourceExtEquiLB(1,0,0,1,iElem) = NodeSourceExt(NodeID(5))
      NodeSourceExtEquiLB(1,1,0,1,iElem) = NodeSourceExt(NodeID(6))
      NodeSourceExtEquiLB(1,1,1,1,iElem) = NodeSourceExt(NodeID(7))
      NodeSourceExtEquiLB(1,0,1,1,iElem) = NodeSourceExt(NodeID(8))
    END DO!iElem
  END IF ! DoDielectricSurfaceCharge

  IF (Do2DSurfaceCharge) THEN
    ! Exchange deposited charge if load balance is performed because the next exchange would be in the hdg solver
    CALL ExchangeSurfNodeSourceMPI()
    SDEALLOCATE(SurfNodeSourceMPI)
    SDEALLOCATE(pq2iNode)
    SDEALLOCATE(Vdm_EQ_N)
    SDEALLOCATE(Vdm_N_EQ)
    ! The root process keeps all the data relevant for sending the surface charge to the other processes
    IF (.NOT.MPIRoot) THEN
      SDEALLOCATE(SurfNodeSource)
      SDEALLOCATE(SurfNodeSymmetryFactor)
      SDEALLOCATE(DepoSurfNodeID2FEMVertexID)
      SDEALLOCATE(FEMVertexID2DepoSurfNodeID)
      SDEALLOCATE(NonUniqueGlobalNodeIDToFEMVertexID)
      SDEALLOCATE(NonUniqueGlobalSideIDToNonUniqueGlobalNodeID)
      SDEALLOCATE(IsDepoSurfSide)
      SDEALLOCATE(IsDepoSurfNode)
      SDEALLOCATE(SurfNodeArea)
    END IF ! .NOT.MPIRoot
  END IF ! Do2DSurfaceCharge


  !! Finalize here because GetCNElemID() is required in this routine for load balancing of NodeSourceExtEquiLB = NodeSourceExt
  !IF (nComputeNodeProcessors.NE.nProcessors_Global) THEN
  !  CALL UNLOCK_AND_FREE(GlobalElem2CNTotalElem_Shared_Win)
  !  ADEALLOCATE(GlobalElem2CNTotalElem)
  !  ADEALLOCATE(GlobalElem2CNTotalElem_Shared)
  !END IF ! nComputeNodeProcessors.NE.nProcessors_Global
ELSE
#endif /*USE_LOADBALANCE*/
  IF (Do2DSurfaceCharge) THEN
    SDEALLOCATE(Vdm_EQ_N)
    SDEALLOCATE(Vdm_N_EQ)
    SDEALLOCATE(SurfNodeSource)
    SDEALLOCATE(SurfNodeSymmetryFactor)
    SDEALLOCATE(DepoSurfNodeID2FEMVertexID)
    SDEALLOCATE(FEMVertexID2DepoSurfNodeID)
    SDEALLOCATE(NonUniqueGlobalNodeIDToFEMVertexID)
    SDEALLOCATE(NonUniqueGlobalSideIDToNonUniqueGlobalNodeID)
    SDEALLOCATE(IsDepoSurfSide)
    SDEALLOCATE(IsDepoSurfNode)
    SDEALLOCATE(SurfNodeArea)
    SDEALLOCATE(pq2iNode)
  END IF ! Do2DSurfaceCharge
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

#if USE_LOADBALANCE
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  ! Keep for load balance and deallocate/reallocate after communication
  SDEALLOCATE(PS_N) ! PartSource, PartSourceOld and PartSourceTmp
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

SDEALLOCATE(DepoNodetoGlobalNode)
SDEALLOCATE(NodeSource)
SDEALLOCATE(NodeSourceExt)

#if USE_MPI
SDEALLOCATE(NodeSendDepoRankToGlobalRank)
SDEALLOCATE(NodeRecvDepoRankToGlobalRank)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(SendRequest)
SDEALLOCATE(NodeSourceExtMPI)
SDEALLOCATE(NodeMappingSend)
SDEALLOCATE(NodeMappingRecv)

SDEALLOCATE(SurfNodeSendDepoRankToGlobalRank)
SDEALLOCATE(SurfNodeRecvDepoRankToGlobalRank)
SDEALLOCATE(SurfRecvRequest)
SDEALLOCATE(SurfSendRequest)
! SDEALLOCATE(NodeSourceExtMPI)
SDEALLOCATE(SurfNodeMappingSend)
SDEALLOCATE(SurfNodeMappingRecv)
#endif /*USE_MPI*/

END SUBROUTINE FinalizeDeposition

#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
END MODULE MOD_PICDepo