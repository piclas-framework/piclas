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
#if USE_PETSC
#include "petsc/finclude/petsc.h"
#endif

!===================================================================================================================================
!> Module for the HDG method
!===================================================================================================================================
MODULE MOD_HDG
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
PUBLIC :: InitHDG,FinalizeHDG
PUBLIC :: HDG, RecomputeEFieldHDG
PUBLIC :: DefineParametersHDG
#if defined(PARTICLES)
PUBLIC :: CalculatePhiAndEFieldFromCurrentsVDL
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> Define parameters for HDG (Hybridized Discontinuous Galerkin)
!===================================================================================================================================
SUBROUTINE DefineParametersHDG()
! MODULES
USE MOD_ReadInTools           ,ONLY: prms
#if defined(PARTICLES)
USE MOD_Part_BR_Elecron_Fluid ,ONLY: DefineParametersBR
#endif /*defined(PARTICLES)*/
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection("HDG")

CALL prms%CreateIntOption(    'HDGNonLinSolver'        ,'Select Newton or Fix Point algorithm (default is Newton and Fix Point is not implemented)', '1')
CALL prms%CreateLogicalOption('NewtonExactSourceDeriv' ,'Use exact derivative of exponential function in source term instead of linearized.', '.FALSE.')
CALL prms%CreateIntOption(    'AdaptIterNewton'        ,'Set number of iteration steps after which the system matrix is recomputed', '0')
CALL prms%CreateLogicalOption('NewtonAdaptStartValue'  ,'Initial recomputation of the system matrix with initial guess of the solution', '.FALSE.')
CALL prms%CreateIntOption(    'AdaptIterNewtonToLinear','Maximum number of iterations where the exact source derivative is used before it is switched to the linearization', '100')
CALL prms%CreateIntOption(    'MaxIterNewton'          ,'Maximum number of iterations in the Newton solver', '10000')
CALL prms%CreateRealOption(   'EpsNonLinear'           ,'Abort residual of the Newton solver', '1.0E-6')
CALL prms%CreateIntOption(    'PrecondType'            ,'Preconditioner type\n 0: no preconditioner\n 1: Side-block SPD preconditioner matrix (already Cholesky decomposed)\n 2: Inverse of diagonal preconditioned')
CALL prms%CreateRealOption(   'epsCG'                  ,'Abort residual of the CG solver', '1.0E-6')
CALL prms%CreateIntOption(    'OutIterCG'              ,'Number of iteration steps between output of CG solver info to std out', '1')
CALL prms%CreateLogicalOption('useRelativeAbortCrit'   ,'Switch between relative and absolute abort criterion', '.FALSE.')
CALL prms%CreateIntOption(    'MaxIterCG'              ,'Maximum number of iterations in the CG solver', '500')
CALL prms%CreateLogicalOption('ExactLambda'            ,'Initially set lambda on all sides (volume and boundaries) via a pre-defined function (ExactFunc)', '.FALSE.')
CALL prms%CreateIntOption(    'HDGSkip'                ,'Number of time step iterations until the HDG solver is called (i.e. all intermediate calls are skipped)', '0')
CALL prms%CreateIntOption(    'HDGSkipInit'            ,'Number of time step iterations until the HDG solver is called (i.e. all intermediate calls are skipped) while time < HDGSkip_t0 (if HDGSkip > 0)', '0')
CALL prms%CreateRealOption(   'HDGSkip_t0'             ,'Time during which HDGSkipInit is used instead of HDGSkip (if HDGSkip > 0)', '0.')
CALL prms%CreateLogicalOption('HDGDisplayConvergence'  ,'Display divergence criteria: Iterations, RunTime and Residual', '.FALSE.')
CALL prms%CreateRealArrayOption( 'EPC-Resistance'      , 'Vector (length corresponds to the number of EPC boundaries) with the resistance for each EPC in Ohm', no=0)
CALL prms%CreateLogicalOption('HDGNSideMin'            ,'Use the minimum polynomial degree at the sides for the HDG solver', '.FALSE.')
#if defined(PARTICLES)
CALL prms%CreateLogicalOption(  'UseBiasVoltage'              , 'Activate usage of bias voltage adjustment (for specific boundaries only)', '.FALSE.')
CALL prms%CreateIntOption(      'BiasVoltage-NPartBoundaries' , 'Number of particle boundaries where the total ion excess is to be calculated for bias voltage model')
CALL prms%CreateIntArrayOption( 'Biasvoltage-PartBoundaries'  , 'Particle boundary index of boundaries where the total ion excess is to be calculated for bias voltage model', no=0)
CALL prms%CreateRealOption(     'BiasVoltage-Frequency'       , 'Frequency of the sinusoidal field boundary where the bias voltage is applied (a value of 0.0 corresponds to a DC potential BC). The total particle electric current over one cycle is required to converge to zero.')
CALL prms%CreateRealOption(     'BiasVoltage-Delta'           , 'Bias voltage difference used for adjusting the DC voltage of the corresponding BC')
#endif /*defined(PARTICLES)*/

! --- BR electron fluid
#if defined(PARTICLES)
CALL DefineParametersBR()
#endif /*defined(PARTICLES)*/

END SUBROUTINE DefineParametersHDG


!===================================================================================================================================
!> Initialize variables of the HDG module
!===================================================================================================================================
SUBROUTINE InitHDG()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_Basis                 ,ONLY: PolynomialDerivativeMatrix
USE MOD_Interpolation_Vars    ,ONLY: N_Inter,NMax
USE MOD_ChangeBasis           ,ONLY: ChangeBasis2D
USE MOD_Elem_Mat              ,ONLY: Elem_Mat,BuildPrecond
USE MOD_ReadInTools           ,ONLY: GETLOGICAL,GETREAL,GETINT
USE MOD_Mesh_Vars             ,ONLY: nBCSides,N_SurfMesh
USE MOD_Mesh_Vars             ,ONLY: BoundaryType,nSides,BC
USE MOD_Mesh_Vars             ,ONLY: nGlobalMortarSides,nMortarMPISides,N_VolMesh
USE MOD_Mesh_Vars             ,ONLY: offSetElem,ElemToSide
USE MOD_Basis                 ,ONLY: InitializeVandermonde,LegendreGaussNodesAndWeights,BarycentricWeights
USE MOD_FillMortar_HDG        ,ONLY: InitMortar_HDG
USE MOD_HDG_Vars              ,ONLY: BRNbrOfRegions,ElemToBRRegion,RegionElectronRef
USE MOD_DG_Vars               ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave
#if defined(PARTICLES)
USE MOD_Part_BR_Elecron_Fluid ,ONLY: UpdateNonlinVolumeFac
USE MOD_Restart_Vars          ,ONLY: DoRestart
#endif /*defined(PARTICLES)*/
#if USE_MPI
USE MOD_MPI_HDG               ,ONLY: StartReceiveMPISurfDataType, StartSendMPISurfDataType, FinishExchangeMPISurfDataType
USE MOD_MPI_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
USE PETSc
USE MOD_HDG_Vars_PETSc
USE MOD_Mesh_Vars             ,ONLY: nMPISides_YOUR
#if USE_MPI
USE MOD_MPI                   ,ONLY: StartReceiveMPIDataInt,StartSendMPIDataInt,FinishExchangeMPIData
#endif /*USE_MPI*/
USE MOD_Elem_Mat              ,ONLY: PETScFillSystemMatrix
USE MOD_HDG_PETSc             ,ONLY: PETScSetSolver
#endif /*USE_PETSC*/
USE MOD_Mesh_Vars             ,ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars             ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_HDG_Init              ,ONLY: InitFPC,InitEPC
#if defined(PARTICLES)
USE MOD_HDG_Init              ,ONLY: InitBV
#endif /*defined(PARTICLES)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,r,iElem,SideID,Nloc,iNeumannBCsides,NSideMin,NSideMax,NSide, iSide
INTEGER           :: BCType,BCState
REAL              :: D(0:Nmax,0:Nmax)
INTEGER           :: nDirichletBCsidesGlobal
#if USE_PETSC
PetscErrorCode    :: ierr
PetscInt          :: major,minor,subminor,release
IS                :: PETScISLocal, PETScISGlobal
INTEGER           :: iUniqueFPCBC
INTEGER             :: iLocalPETScDOF,iDOF
INTEGER             :: OffsetCounter
INTEGER,ALLOCATABLE :: localToGlobalPETScDOF(:)
#if USE_MPI
INTEGER             :: iProc,PETScDOFOffsetsMPI(nProcessors)
#endif
#endif
INTEGER           :: locSide,nMortars
INTEGER           :: MortarSideID,iMortar
REAL              :: StartT,EndT
#if USE_PETSC
CHARACTER(100)    :: hilf
#endif /*USE_PETSC*/
!===================================================================================================================================
IF(HDGInitIsDone)THEN
   LBWRITE(*,*) "InitHDG already called."
   RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT HDG...'
GETTIME(StartT)

! ----------------------------------------------------------------------------------------------------------------------------------
! MAIN STEPS
! ----------------------------------------------------------------------------------------------------------------------------------
!  1. Calculate number of gauss points in sides and volumes for each PP_N
!  2. Allocate SurfExchange for communicating surface information with variable PP_N
!  3. Build SurfElemMin for all sides (including Mortar sides)
!  4. Initialize BR electron fluid model
!  5. Init mortar stuff
!  6. BCs, the first
!  7. Init floating potentials, electric potentials and bias voltage
!  8. BCs the second...
!  9. Initialize interpolation variables for each Polynomial degree (Also fill HDG_Vol_N further)
! 10. Allocate and zero missing HDG_VOL_N and HDG_Surf_N stuff
! 11. Build Preconditioner
! 12. Initialize PETSc
! ----------------------------------------------------------------------------------------------------------------------------------

! (0. Read generic HDG Settings)
HDGDisplayConvergence = GETLOGICAL('HDGDisplayConvergence')

HDGSkip = GETINT('HDGSkip')
IF (HDGSkip.GT.0) THEN
  HDGSkipInit = GETINT('HDGSkipInit')
  HDGSkip_t0  = GETREAL('HDGSkip_t0')
ELSE
  HDGSkip=0
END IF

! Read in CG parameters (also used for PETSc)
#if USE_PETSC
PetscCallA(PetscGetVersionNumber(major,minor,subminor,release,ierr))
#ifdef PETSC_HAVE_HYPRE
hilf = '(built with Hypre and'
#else
hilf = '(built without Hypre and'
#endif /*PETSC_HAVE_HYPRE*/
#ifdef PETSC_HAVE_MUMPS
hilf = TRIM(hilf)//' with Mumps)'
#else
hilf = TRIM(hilf)//' without Mumps)'
#endif /*PETSC_HAVE_MUMPS*/
LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A,I0,A)') ' | Method for HDG solver: PETSc ',major,'.',minor,'.',subminor,' '//TRIM(hilf)
PrecondType          = GETINT('PrecondType','1')
#else /*without PETSC*/
LBWRITE(UNIT_stdOut,'(A)') ' | Method for HDG solver: CG '
PrecondType          = GETINT('PrecondType','2')
#endif /*USE_PETSC*/
epsCG                = GETREAL('epsCG')
OutIterCG            = GETINT('OutIterCG')
useRelativeAbortCrit = GETLOGICAL('useRelativeAbortCrit')
MaxIterCG            = GETINT('MaxIterCG')
ExactLambda          = GETLOGICAL('ExactLambda')
UseNSideMin          = GETLOGICAL('HDGNSideMin')

! Initialize element containers
ALLOCATE(HDG_Vol_N(1:PP_nElems))
ALLOCATE(HDG_Surf_N(1:nSides))

ALLOCATE(MaskedSide(1:nSides))
MaskedSide=0

!mappings (Only used in ELEM_MAT!)
sideDir(  XI_MINUS)=1
sideDir(   XI_PLUS)=1
sideDir( ETA_MINUS)=2
sideDir(  ETA_PLUS)=2
sideDir(ZETA_MINUS)=3
sideDir( ZETA_PLUS)=3
pm(  XI_MINUS)=1
pm(   XI_PLUS)=2
pm( ETA_MINUS)=1
pm(  ETA_PLUS)=2
pm(ZETA_MINUS)=1
pm( ZETA_PLUS)=2

dirPm2iSide(1,1) = XI_MINUS
dirPm2iSide(2,1) = XI_PLUS
dirPm2iSide(1,2) = ETA_MINUS
dirPm2iSide(2,2) = ETA_PLUS
dirPm2iSide(1,3) = ZETA_MINUS
dirPm2iSide(2,3) = ZETA_PLUS

! -------------------------------------------------------------------------------------------------
! Fill NSide for each side.
! For Mortars, DG_Elems_slave(iSide) = -1 for the large mortar side.
! -> Loop over all large mortar sides and fill NSide with the min/max of all shared sides.
DO iSide = 1, nSides
  IF(UseNSideMin) THEN
    NSideMin = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
    N_SurfMesh(iSide)%NSide = NSideMin
  ElSE
    NSideMax = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
    N_SurfMesh(iSide)%NSide = NSideMax
  END IF
END DO

! Mortars: Get minimum/maximum of all sides the Mortar interface
DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
  nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide  = MortarType(2,MortarSideID)
  NSide = N_SurfMesh(MortarSideID)%NSide
  DO iMortar = 1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,locSide) !small SideID
    IF(UseNSideMin) THEN
      NSide = MIN(NSide,N_SurfMesh(SideID)%NSide)
    ELSE
      NSide = MAX(NSide,N_SurfMesh(SideID)%NSide)
    END IF
  END DO !iMortar

  N_SurfMesh(MortarSideID)%NSide = NSide
  DO iMortar = 1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,locSide) !small SideID
    N_SurfMesh(SideID)%NSide = NSide
  END DO !iMortar
END DO !MortarSideID


! -------------------------------------------------------------------------------------------------
! 1. Calculate number of gauss points in sides and volumes for each PP_N
ALLOCATE(nGP_vol(1:NMax))
ALLOCATE(nGP_face(1:NMax))
DO Nloc = 1, NMax
  nGP_vol(Nloc)  = (Nloc+1)**3
  nGP_face(Nloc) = (Nloc+1)**2
END DO ! Nloc = 1, NMax

! -------------------------------------------------------------------------------------------------
! 2. Allocate SurfExchange for communicating surface information with variable PP_N
#if USE_MPI
ALLOCATE(SurfExchange(nNbProcs))
DO iNbProc=1,nNbProcs
  ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv(MAXVAL(DataSizeSurfRecMax(iNbProc,:))))
  ALLOCATE(SurfExchange(iNbProc)%SurfDataSend(MAXVAL(DataSizeSurfSendMax(iNbProc,:))))
END DO !iProc=1,nNBProcs
DO iNbProc=1,nNbProcs
  DEALLOCATE(SurfExchange(iNbProc)%SurfDataRecv)
  DEALLOCATE(SurfExchange(iNbProc)%SurfDataSend)
  IF(UseNSideMin)THEN
    ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv(MAXVAL(DataSizeSurfRecMin( iNbProc,:))))
    ALLOCATE(SurfExchange(iNbProc)%SurfDataSend(MAXVAL(DataSizeSurfSendMin(iNbProc,:))))
    ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv2(MAXVAL(DataSizeSurfRecMin( iNbProc,:))**2))
    ALLOCATE(SurfExchange(iNbProc)%SurfDataSend2(MAXVAL(DataSizeSurfSendMin(iNbProc,:))**2))
  ELSE
    ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv(MAXVAL(DataSizeSurfRecMax( iNbProc,:))))
    ALLOCATE(SurfExchange(iNbProc)%SurfDataSend(MAXVAL(DataSizeSurfSendMax(iNbProc,:))))
    ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv2(MAXVAL(DataSizeSurfRecMax( iNbProc,:))**2))
    ALLOCATE(SurfExchange(iNbProc)%SurfDataSend2(MAXVAL(DataSizeSurfSendMax(iNbProc,:))**2))
  END IF
END DO !iProc=1,nNBProcs
#endif /*USE_MPI*/

! -------------------------------------------------------------------------------------------------


! 4. Initialize BR electron fluid model
HDGNonLinSolver = -1 ! init
#if defined(PARTICLES)
! BR electron fluid model
IF (BRNbrOfRegions .GT. 0) THEN !Regions only used for Boltzmann Electrons so far -> non-linear HDG-sources!
#if USE_PETSC
  CALL CollectiveStop(__STAMP__,' HDG with BR electron fluid (non-linear HDG solver) is not implemented with PETSc')
#endif /*USE_PETSC*/
  HDGNonLinSolver=GETINT('HDGNonLinSolver')

  IF (HDGNonLinSolver.EQ.1) THEN
    NewtonExactSourceDeriv  = GETLOGICAL('NewtonExactSourceDeriv')
    AdaptIterNewton         = GETINT('AdaptIterNewton')
    AdaptIterNewtonOld      = AdaptIterNewton
    NewtonAdaptStartValue   = GETLOGICAL('NewtonAdaptStartValue')
    AdaptIterNewtonToLinear = GETINT('AdaptIterNewtonToLinear')
    IF (NewtonExactSourceDeriv) NewtonAdaptStartValue=.TRUE.
    IF (DoRestart) NewtonAdaptStartValue=.FALSE.
    !ALLOCATE(NonlinVolumeFac(nGP_vol,PP_nElems))
    DO iElem=1,PP_nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      ALLOCATE(HDG_Vol_N(iElem)%NonlinVolumeFac(nGP_vol(Nloc)))
    END DO

    ! Set NonlinVolumeFac = RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0) for each element.
    ! Set zero if
    !  a) NewtonAdaptStartValue is activated
    !  b) fully kinetic currently active and variable electron temperature activated (skips calculation when starting from
    !     a simulation that was fully kinetic and will be switched)
    CALL UpdateNonlinVolumeFac(NewtonAdaptStartValue.OR.((.NOT.UseBRElectronFluid).AND.CalcBRVariableElectronTemp))

  END IF

  MaxIterNewton = GETINT('MaxIterNewton')
  EpsNonLinear  = GETREAL('EpsNonLinear')
END IF
#endif /*defined(PARTICLES)*/

! 5. Init mortar stuff
IF(nGlobalMortarSides.GT.0)THEN !mortar mesh
  IF(nMortarMPISides.GT.0) CALL abort(__STAMP__,&
  "nMortarMPISides >0: HDG mortar MPI implementation relies on big sides having always only master sides (=> nMortarMPISides=0 )")

CALL InitMortar_HDG()
END IF !mortarMesh

! 6. BCs, the first
!boundary conditions
nDirichletBCsides=0
nNeumannBCsides  =0
nConductorBCsides=0
DO SideID=1,nBCSides
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(HDGDIRICHLETBCSIDEIDS) ! Dirichlet
    nDirichletBCsides=nDirichletBCsides+1
  CASE(10,11,12) ! Neumann
    nNeumannBCsides=nNeumannBCsides+1
  CASE(20) ! Conductor: Floating Boundary Condition (FPC)
    nConductorBCsides=nConductorBCsides+1
  CASE DEFAULT ! unknown BCType
    CALL CollectiveStop(__STAMP__,' unknown BC Type in hdg.f90!',IntInfo=BCType)
  END SELECT ! BCType
END DO

! 7. Init floating potentials, electric potentials and bias voltage
! Conductor: Initialize floating boundary condition
CALL InitFPC()

! Conductor: Initialize electric potential condition (resistive decharging of surface and electric potential calculation)
CALL InitEPC()

#if defined(PARTICLES)
! Bias Voltage: Initialize containers and sub-communicator
! BCType: 50,X for bias voltage + DC boundary condition
! BCType: 51,X for bias voltage + cos(wt) function boundary condition
! BCType: 52,X for bias voltage + cos(wt) function + coupled power adjustment (for AC and not DC in this case)
CALL InitBV()
#endif /*defined(PARTICLES)*/

! 8. BCs the second...
! Get the global number of Dirichlet boundaries. If there are none, the potential of a single DOF must be set.
#if USE_MPI
  CALL MPI_ALLREDUCE(nDirichletBCsides , nDirichletBCsidesGlobal , 1 , MPI_INTEGER , MPI_MAX , MPI_COMM_PICLAS , IERROR)
#else
  nDirichletBCsidesGlobal = nDirichletBCsides
#endif /*USE_MPI*/

ZeroPotentialSide = -1
IF(mpiRoot.AND.nDirichletBCsidesGlobal==0) ZeroPotentialSide = ElemToSide(E2S_SIDE_ID,1,1)

IF(nDirichletBCsides.GT.0)ALLOCATE(DirichletBC(nDirichletBCsides))
IF(nNeumannBCsides  .GT.0)THEN
  ALLOCATE(NeumannBC(nNeumannBCsides))
END IF
IF(nConductorBCsides.GT.0)ALLOCATE(ConductorBC(nConductorBCsides))
#if (PP_nVar!=1)
  IF(nDirichletBCsides.GT.0)ALLOCATE(qn_face_MagStat(PP_nVar, nGP_face(PP_N),nDirichletBCsides))
#endif
nDirichletBCsides=0
nNeumannBCsides  =0
nConductorBCsides=0
DO SideID=1,nBCSides
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(HDGDIRICHLETBCSIDEIDS) ! Dirichlet
    nDirichletBCsides=nDirichletBCsides+1
    DirichletBC(nDirichletBCsides)=SideID
    MaskedSide(SideID)=1
  CASE(10,11,12) !Neumann,
    nNeumannBCsides=nNeumannBCsides+1
    NeumannBC(nNeumannBCsides)=SideID
  CASE(20) ! Conductor: Floating Boundary Condition (FPC)
    nConductorBCsides=nConductorBCsides+1
    ConductorBC(nConductorBCsides)=SideID
    MaskedSide(SideID)=2
  CASE DEFAULT ! unknown BCType
    CALL CollectiveStop(__STAMP__,' unknown BC Type in hdg.f90!',IntInfo=BCType)
  END SELECT ! BCType
END DO

IF(nNeumannBCsides.GT.0)THEN
  DO iNeumannBCsides = 1, nNeumannBCsides
    SideID = NeumannBC(iNeumannBCsides)
    Nloc = N_SurfMesh(SideID)%NSide
    ALLOCATE(HDG_Surf_N(SideID)%qn_face(PP_nVar, nGP_face(Nloc)))
  END DO ! iNeumannBCsides = 1, nNeumannBCsides
  END IF

! 9. Initialize interpolation variables for each Polynomial degree (Also fill HDG_Vol_N further)
! Initialize interpolation variables
DO Nloc = 1, NMax
  ALLOCATE(N_Inter(Nloc)%LL_minus(0:Nloc,0:Nloc))
  ALLOCATE(N_Inter(Nloc)%LL_plus( 0:Nloc,0:Nloc))
  DO j=0,Nloc
    DO i=0,Nloc
      N_Inter(Nloc)%LL_minus(i,j) = N_Inter(Nloc)%L_minus(i)*N_Inter(Nloc)%L_minus(j)
      N_Inter(Nloc)%LL_plus(i,j)  = N_Inter(Nloc)%L_plus(i)*N_Inter(Nloc)%L_plus(j)
    END DO
  END DO

  ALLOCATE(N_Inter(Nloc)%Lomega_m(0:Nloc))
  ALLOCATE(N_Inter(Nloc)%Lomega_p(0:Nloc))
! Compute a lifting matrix scaled by the Gaussian weights
  N_Inter(Nloc)%Lomega_m = - N_Inter(Nloc)%L_minus/N_Inter(Nloc)%wGP
  N_Inter(Nloc)%Lomega_p = + N_Inter(Nloc)%L_plus/N_Inter(Nloc)%wGP
  ALLOCATE(N_Inter(Nloc)%Domega(0:Nloc,0:Nloc))
! Compute Differentiation matrix D for given Gausspoints (1D)
  CALL PolynomialDerivativeMatrix(Nloc,N_Inter(Nloc)%xGP,D(0:Nloc,0:Nloc))
  ! Compute a Differentiation matrix scaled by the Gaussian weights
  DO j=0,Nloc
    DO i=0,Nloc
      N_Inter(Nloc)%Domega(i,j) = N_Inter(Nloc)%wGP(i)/N_Inter(Nloc)%wGP(j)*D(i,j)
  END DO !r
END DO !s

  ALLOCATE(N_Inter(Nloc)%wGP_vol(nGP_vol(Nloc)))
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
        N_Inter(Nloc)%wGP_vol(r)=N_Inter(Nloc)%wGP(i)*N_Inter(Nloc)%wGP(j)*N_Inter(Nloc)%wGP(k)
      END DO
    END DO
  END DO !i,j,k
END DO ! Nloc = 1, NMax

DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(HDG_Vol_N(iElem)%InvDhat(nGP_vol(Nloc),nGP_vol(Nloc)))
  ALLOCATE(HDG_Vol_N(iElem)%JwGP_vol(nGP_vol(Nloc)))
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
    HDG_Vol_N(iElem)%JwGP_vol(r)=N_Inter(Nloc)%wGP_vol(r)/N_VolMesh(iElem)%sJ(i,j,k) !omega*J
  END DO; END DO; END DO !i,j,k

  ALLOCATE(HDG_Vol_N(iElem)%Ehat(nGP_face(Nloc),nGP_vol(Nloc),6))
!side matrices
  ALLOCATE(HDG_Vol_N(iElem)%Smat(nGP_face(Nloc),nGP_face(Nloc),6,6))
END DO !iElem

!stabilization parameter
ALLOCATE(Tau(PP_nElems))
DO iElem=1,PP_nElems
  Tau(iElem)=2./((SUM(HDG_Vol_N(iElem)%JwGP_vol(:)))**(1./3.))  !1/h ~ 1/vol^(1/3) (volref=8)
END DO !iElem

CALL Elem_Mat(0_8) ! takes iter=0 (kind=8)

! 10. Allocate and zero missing HDG_VOL_N and HDG_Surf_N stuff
DO iElem = 1, PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(HDG_Vol_N(iElem)%RHS_vol(PP_nVar, nGP_vol(Nloc)))
  HDG_Vol_N(iElem)%RHS_vol=0.
END DO ! iElem = 1, PP_nElems

DO SideID = 1, nSides
  NSide = N_SurfMesh(SideID)%NSide
  ALLOCATE(HDG_Surf_N(SideID)%lambda(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%lambda=0.
  ALLOCATE(HDG_Surf_N(SideID)%RHS_face(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%RHS_face=0.
  ALLOCATE(HDG_Surf_N(SideID)%mv(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%mv=0.
  ALLOCATE(HDG_Surf_N(SideID)%R(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%R=0.
  ALLOCATE(HDG_Surf_N(SideID)%V(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%V=0.
  ALLOCATE(HDG_Surf_N(SideID)%Z(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%Z=0.
#if USE_MPI
  ALLOCATE(HDG_Surf_N(SideID)%buf(PP_nVar,nGP_face(NSide)))
  HDG_Surf_N(SideID)%buf=0.
#if USE_PETSC
#else
  IF(PrecondType.EQ.1)THEN
    ALLOCATE(HDG_Surf_N(SideID)%buf2(nGP_face(NSide),nGP_face(NSide)))
    HDG_Surf_N(SideID)%buf2=0.
  END IF ! PrecondType.EQ.1
#endif /*USE_PETSC*/
#endif
END DO ! SideID = 1, nSides

! 11. Build Preconditioner
! Requires HDG_Surf_N(SideID)%buf
#if !USE_PETSC
CALL BuildPrecond()
#endif

#if USE_PETSC
! -------------------------------------------------------------------------------------------------
! 12. Initialize PETSc
! -------------------------------------------------------------------------------------------------
! Steps:
! 3.1) Create PETSc mappings to build the global system
! 3.2) Initialize PETSc objects

! 3.1) Create PETSc mappings to build the global system
! 3.1.1) Calculate nLocalPETScDOFs without nMPISides_YOUR to compute nGlobalPETScDOFs
nLocalPETScDOFs = 0
DO SideID=1,nSides-nMPISides_YOUR
  IF(MaskedSide(SideID).NE.0) CYCLE ! Skip Dirichlet + small mortar sides
  nLocalPETScDOFs = nLocalPETScDOFs + nGP_face(N_SurfMesh(SideID)%NSide)
END DO

! 3.1.2) Calculate nGlobalPETScDOFs
! This is the total number of PETScDOFs
#if USE_MPI
CALL MPI_ALLGATHER(nLocalPETScDOFs,1,MPI_INTEGER,PETScDOFOffsetsMPI,1,MPI_INTEGER,MPI_COMM_PICLAS,IERROR)
nGlobalPETScDOFs = SUM(PETScDOFOffsetsMPI)
#else
nGlobalPETScDOFs = nLocalPETScDOFs
#endif

! 3.1.3) Calculate OffsetGlobalPETScDOF(SideID)
! This is the GlobalPETScDOF of the first DOF of the side with given SideID
OffsetCounter = 0
#if USE_MPI
DO iProc=1, myrank
  OffsetCounter = OffsetCounter + PETScDOFOffsetsMPI(iProc)
END DO
#endif
ALLOCATE(OffsetGlobalPETScDOF(nSides))
DO SideID=1,nSides-nMPISides_YOUR
  IF(MaskedSide(SideID).NE.0) CYCLE ! Skip Dirichlet + Small mortar sides
  OffsetGlobalPETScDOF(SideID) = OffsetCounter
  OffsetCounter = OffsetCounter + nGP_face(N_SurfMesh(SideID)%NSide)
END DO

! 3.1.3.1 Mortars: Small mortar sides have the same DOFs as the big side
DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
  nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide  = MortarType(2,MortarSideID)
  DO iMortar = 1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,locSide)
    OffsetGlobalPETScDOF(SideID) = OffsetGlobalPETScDOF(MortarSideID)
  END DO !iMortar
END DO !MortarSideID

#if USE_MPI
CALL StartReceiveMPIDataInt(1,OffsetGlobalPETScDOF,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
CALL StartSendMPIDataInt(   1,OffsetGlobalPETScDOF,1,nSides,SendRequest_U,SendID=1) ! Send MINE
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif

! 4.2.4.3) ZeroPotential
ZeroPotentialDOF = -1
IF(nDirichletBCsidesGlobal==0) THEN
  IF(mpiRoot) ZeroPotentialDOF = OffsetGlobalPETScDOF(ZeroPotentialSide)
#if USE_MPI
  CALL MPI_BCAST(ZeroPotentialDOF,1,MPI_INTEGER,0,MPI_COMM_PICLAS,IERROR)
#endif
END IF

! 3.1.3.5) Add All Small Mortar Sides to nLocalPETScDOFs
DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
  nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  nLocalPETScDOFs = nLocalPETScDOFs + nGP_face(N_SurfMesh(MortarSideID)%NSide) * nMortars
END DO !MortarSideID

! 3.1.4) Sum up YOUR sides for nLocalPETScDOFs
! The full nLocalPETScDOFs is used to compute the Scatter context
#if USE_MPI
! 4. Add the nMPISides_YOUR to nLocalPETScDOFs
DO SideID=nSides-nMPISides_YOUR+1,nSides
  IF(MaskedSide(SideID).GT.0) CYCLE
  nLocalPETScDOFs = nLocalPETScDOFs + nGP_face(N_SurfMesh(SideID)%NSide)
END DO
#endif

! 3.1.5) Create localToGlobalPETScDOF(iLocalPETScDOF) mapping
! The mapping is used to create the PETSc Scatter context to extract the local solution from the global solution vector
ALLOCATE(localToGlobalPETScDOF(nLocalPETScDOFs+FPC%nUniqueFPCBounds))
iLocalPETScDOF = 0
DO SideID=1,nSides
  IF(MaskedSide(SideID).GT.0) CYCLE
  DO iDOF=1,nGP_face(N_SurfMesh(SideID)%NSide)
    iLocalPETScDOF = iLocalPETScDOF + 1
    LocalToGlobalPETScDOF(iLocalPETScDOF) = OffsetGlobalPETScDOF(SideID) + iDOF - 1
  END DO
END DO

! 3.1.6) Add each FPC to the DOFs
IF(UseFPC)THEN
  DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
    LocalToGlobalPETScDOF(nLocalPETScDOFs+iUniqueFPCBC) = nGlobalPETScDOFs + iUniqueFPCBC - 1
  END DO
  nLocalPETScDOFs = nLocalPETScDOFs + FPC%nUniqueFPCBounds
  nGlobalPETScDOFs = nGlobalPETScDOFs + FPC%nUniqueFPCBounds
END IF

! ------------------------------------------------------
! 3.2) Initialize PETSc Objects
PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))

! 3.2.1) Set up and fill System matrix
PetscCallA(MatCreate(PETSC_COMM_WORLD,PETScSystemMatrix,ierr))
PetscCallA(MatSetSizes(PETScSystemMatrix,PETSC_DECIDE,PETSC_DECIDE,nGlobalPETScDOFs,nGlobalPETScDOFs,ierr))
PetscCallA(MatSetType(PETScSystemMatrix,MATSBAIJ,ierr)) ! Symmetric sparse matrix
! Conservative guess for the number of nonzeros: With mortars at most 12 sides with Nmax.
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 21)
PetscCallA(MatSEQSBAIJSetPreallocation(PETScSystemMatrix,1,22 * nGP_face(NMax),PETSC_NULL_INTEGER_ARRAY,ierr))
PetscCallA(MatMPISBAIJSetPreallocation(PETScSystemMatrix,1,22 * nGP_face(NMax),PETSC_NULL_INTEGER_ARRAY,22 * nGP_face(NMax),PETSC_NULL_INTEGER_ARRAY,ierr))
#else
PetscCallA(MatSEQSBAIJSetPreallocation(PETScSystemMatrix,1,22 * nGP_face(NMax),PETSC_NULL_INTEGER,ierr))
PetscCallA(MatMPISBAIJSetPreallocation(PETScSystemMatrix,1,22 * nGP_face(NMax),PETSC_NULL_INTEGER,22 * nGP_face(NMax),PETSC_NULL_INTEGER,ierr))
#endif
PetscCallA(MatZeroEntries(PETScSystemMatrix,ierr))
PetscCallA(MatSetOption(PETScSystemMatrix,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)) ! Column oriented for more convenient set up

CALL PETScFillSystemMatrix()

! 3.2.2) Set up Solver
CALL PETScSetSolver()

! 3.2.3) Set up RHS and solution vectors
PetscCallA(VecCreate(PETSC_COMM_WORLD,PETScSolution,ierr))
PetscCallA(VecSetSizes(PETScSolution,PETSC_DECIDE,nGlobalPETScDOFs,ierr))
PetscCallA(VecSetType(PETScSolution,VECSTANDARD,ierr))
PetscCallA(VecSetUp(PETScSolution,ierr))
PetscCallA(VecDuplicate(PETScSolution,PETScRHS,ierr))

! 3.2.4) Set up Scatter stuff
PetscCallA(VecCreateSeq(PETSC_COMM_SELF,nLocalPETScDOFs,PETScSolutionLocal,ierr))
! Create a PETSc Vector 0:(nLocalPETScDOFs-1)
PetscCallA(ISCreateStride(PETSC_COMM_SELF,nLocalPETScDOFs,0,1,PETScISLocal,ierr))
! Create a PETSc Vector of the Global DOF IDs
PetscCallA(ISCreateGeneral(PETSC_COMM_WORLD,nLocalPETScDOFs,localToGlobalPETScDOF,PETSC_COPY_VALUES,PETScISGlobal,ierr))
! Create a scatter context to extract the local dofs
PetscCallA(VecScatterCreate(PETScSolution,PETScISGlobal,PETScSolutionLocal,PETScISLocal,PETScScatter,ierr))

! (Delete local allocated vectors)
DEALLOCATE(localToGlobalPETScDOF)
! Clean-up local PETSc objects
PetscCallA(ISDestroy(PETScISLocal,ierr))
PetscCallA(ISDestroy(PETScISGlobal,ierr))
#endif /*USE_PETSC*/

GETTIME(EndT)
HDGInitIsDone = .TRUE.
LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' INIT HDG'
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!')
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitHDG


!===================================================================================================================================
!> HDG solver for linear or non-linear systems
!===================================================================================================================================
#if defined(PARTICLES)
SUBROUTINE HDG(t,iter,ForceCGSolverIteration_opt,RecomputeLambda_opt)
#else
SUBROUTINE HDG(t,iter,RecomputeLambda_opt)
#endif /*defined(PARTICLES)*/
! MODULES
#if USE_MPI
USE mpi_f08
USE MOD_HDG_Readin      ,ONLY: synchronizevoltageonepc
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars   ,ONLY: iStage
#endif
#if (USE_HDG && (PP_nVar==1))
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars   ,ONLY: iStage,nRKStages
#endif
#if defined(PARTICLES)
USE MOD_TimeDisc_Vars   ,ONLY: dt
#endif /*defined(PARTICLES)*/
USE MOD_Analyze_Vars    ,ONLY: CalcElectricTimeDerivative
#endif /*(USE_HDG && (PP_nVar==1))*/
#if defined(PARTICLES)
USE MOD_HDG_Vars        ,ONLY: UseEPC,EPC
USE MOD_HDG_NonLinear   ,ONLY: HDGNewton
#endif /*defined(PARTICLES)*/
USE MOD_HDG_Linear      ,ONLY: HDGLinear
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: t !time
INTEGER(KIND=8),INTENT(IN)  :: iter
#if defined(PARTICLES)
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only required for BR electron fluid)
#endif /*defined(PARTICLES)*/
LOGICAL,INTENT(IN),OPTIONAL :: RecomputeLambda_opt        ! Is set true during restart and is used to skip the calculation of dD/dt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(PARTICLES)
LOGICAL :: ForceCGSolverIteration_loc
INTEGER :: iUniqueEPCBC
#endif /*defined(PARTICLES)*/
LOGICAL :: RecomputeLambda_loc
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(4))
#endif /*EXTRAE*/

IF (PRESENT(RecomputeLambda_opt)) THEN; RecomputeLambda_loc = RecomputeLambda_opt
ELSE;                                   RecomputeLambda_loc = .FALSE.
END IF

#if (USE_HDG && (PP_nVar==1))
! Calculate temporal derivate of D in last iteration before Analyze_dt is reached: Store E^n here
IF(CalcElectricTimeDerivative.AND.(.NOT.RecomputeLambda_loc)) CALL CalculateElectricTimeDerivative(iter,1)
#endif /*(USE_HDG && (PP_nVar==1))*/

! Check whether the solver should be skipped in this iteration
IF ((iter.GT.0).AND.( HDGSkip.NE.0)) THEN
  IF (t.LT.HDGSkip_t0) THEN
    IF (MOD(iter,INT(HDGSkipInit,8)).NE.0) RETURN
  ELSE
    IF (MOD(iter,INT(HDGSkip,8)).NE.0) RETURN
  END IF
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  IF (iStage.GT.1) RETURN
#endif
END IF

  ! Electric potential condition (EPC)
#if (USE_HDG && (PP_nVar==1) && defined(PARTICLES))
  IF(UseEPC)THEN
#if USE_MPI
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
    IF(iStage.EQ.1)&
#endif
      EPC%Charge = 0.

    ! Communicate the accumulated charged on each BC to all processors on the communicator
    DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
      IF(EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
        IF(MPIRoot)THEN
          CALL MPI_REDUCE(MPI_IN_PLACE, EPC%ChargeProc(iUniqueEPCBC), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, EPC%COMM(iUniqueEPCBC)%UNICATOR, IERROR)
        ELSE
          CALL MPI_REDUCE(EPC%ChargeProc(iUniqueEPCBC), 0           , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, EPC%COMM(iUniqueEPCBC)%UNICATOR, IERROR)
        END IF ! MPIRoot
        EPC%Charge(iUniqueEPCBC) = EPC%Charge(iUniqueEPCBC) + EPC%ChargeProc(iUniqueEPCBC)
      END IF ! EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL
    END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
#endif /*USE_MPI*/
    IF(MPIRoot) EPC%Charge(:) = EPC%Charge(:) + EPC%ChargeProc(:)
    EPC%ChargeProc = 0.
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
    IF(iStage.EQ.nRKStages) THEN
#endif
      ! Only root calculates the voltage and then broadcasts it
      IF(MPIRoot)THEN
        DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
          ! U = R*I = R*(-dQ/dt)
          EPC%Voltage(iUniqueEPCBC) = EPC%Resistance(iUniqueEPCBC)*(- EPC%Charge(iUniqueEPCBC) / dt)
        END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
      END IF ! MPIRoot
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
    END IF !iStage.EQ.nRKStages
#endif
#if USE_MPI
  ! 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
  CALL SynchronizeVoltageOnEPC()
#endif /*USE_MPI*/
  END IF ! UseEPC
#endif /*(USE_HDG && (PP_nVar==1)&& defined(PARTICLES))*/

! Run the appropriate HDG solver: Newton or Linear
#if defined(PARTICLES)
IF(UseBRElectronFluid) THEN
  IF (HDGNonLinSolver.EQ.1) THEN

    IF (PRESENT(ForceCGSolverIteration_opt)) THEN; ForceCGSolverIteration_loc = ForceCGSolverIteration_opt
    ELSE;                                          ForceCGSolverIteration_loc = .FALSE.
    END IF

    !CALL HDGNewton(t, U_out, iter, ForceCGSolverIteration_loc)
    CALL HDGNewton(t, iter, ForceCGSolverIteration_loc)
  ELSE
    CALL CollectiveStop(__STAMP__,'Defined HDGNonLinSolver not implemented (HDGFixPoint has been removed!) HDGNonLinSolver = ',&
    IntInfo=HDGNonLinSolver)
  END IF
ELSE
#endif /*defined(PARTICLES)*/
  CALL HDGLinear(t)
#if defined(PARTICLES)
END IF
#endif /*defined(PARTICLES)*/

#if (USE_HDG && (PP_nVar==1))
! Calculate temporal derivate of D in last iteration before Analyze_dt is reached: Use E^n+1 here and calculate the derivative dD/dt
IF(CalcElectricTimeDerivative.AND.(.NOT.RecomputeLambda_loc)) CALL CalculateElectricTimeDerivative(iter,2)
#endif /*(USE_HDG && (PP_nVar==1))*/

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
END SUBROUTINE HDG


#if (USE_HDG && (PP_nVar==1))
!===================================================================================================================================
!> Calculates the temporal derivate of the electric displacement field strength
!===================================================================================================================================
SUBROUTINE CalculateElectricTimeDerivative(iter,mode)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Globals_Vars           ,ONLY: eps0
USE MOD_TimeDisc_Vars          ,ONLY: dt,dt_Min
USE MOD_Analyze_Vars           ,ONLY: FieldAnalyzeStep
USE MOD_DG_Vars                ,ONLY: U_N
USE MOD_Dielectric_vars        ,ONLY: DoDielectric,isDielectricElem,ElemToDielectric,DielectricVol
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars          ,ONLY: iStage,nRKStages
#endif /*(PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)*/
USE MOD_Particle_Boundary_Vars ,ONLY: DoVirtualDielectricLayer
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)  :: iter
INTEGER,INTENT(IN) :: mode !< 1: store E^n at the beginning of the time step
                           !< 2: store E^n+1 at the end of the time step and subtract E^n to calculate the difference
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iDir,iElem
!===================================================================================================================================
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
IF((iStage.NE.1).AND.(iStage.NE.nRKStages)) RETURN
#endif

! iter is incremented after this function and then checked in analyze routine with iter+1
IF( ( ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR. & ! Analysis dt
      ALMOSTEQUAL(dt,dt_Min(DT_END)).OR.     & ! tEnd is reached
      (MOD(iter+1,FieldAnalyzeStep).EQ.0)    & ! Field analysis is reached
    ) .OR. DoVirtualDielectricLayer)THEN       ! Exception: VDL requires analysis in every time step (integration of ODE)
  IF(mode.EQ.1)THEN
    ! Store E^n at the beginning of the time step
    DO iElem = 1, nElems
      U_N(iElem)%Dt(:,:,:,:) = U_N(iElem)%E(:,:,:,:)
    END DO ! iElem = 1, nElems
  ELSE
    ! Store E^n+1 at the end of the time step and subtract E^n to calculate the difference
    IF(DoDielectric)THEN
      DO iElem=1,PP_nElems
        IF(isDielectricElem(iElem)) THEN
          DO iDir = 1, 3
            U_N(iElem)%Dt(iDir,:,:,:) = DielectricVol(ElemToDielectric(iElem))%DielectricEps(:,:,:)&
                *eps0*(U_N(iElem)%E(iDir,:,:,:)-U_N(iElem)%Dt(iDir,:,:,:)) / dt
          END DO ! iDir = 1, 3
        ELSE
          U_N(iElem)%Dt(:,:,:,:) = eps0*(U_N(iElem)%E(:,:,:,:)-U_N(iElem)%Dt(:,:,:,:)) / dt
        END IF ! isDielectricElem(iElem)
      END DO ! iElem=1,PP_nElems
    ELSE
      DO iElem=1,PP_nElems
        U_N(iElem)%Dt(:,:,:,:) = eps0*(U_N(iElem)%E(:,:,:,:)-U_N(iElem)%Dt(:,:,:,:)) / dt
      END DO ! iElem=1,PP_nElems
    END IF ! DoDielectric
#if defined(PARTICLES)
    ! Calculate the electric VDL surface potential from the particle and electric displacement current
    IF(DoVirtualDielectricLayer) CALL CalculatePhiAndEFieldFromCurrentsVDL(.TRUE.)
#endif /*defined(PARTICLES)*/
  END IF ! mode.EQ.1
END IF

END SUBROUTINE CalculateElectricTimeDerivative


#if defined(PARTICLES)
!===================================================================================================================================
!> description
!===================================================================================================================================
SUBROUTINE CalculatePhiAndEFieldFromCurrentsVDL(UpdatePhiF)
! MODULES
USE MOD_Globals                ,ONLY: VECNORM
USE MOD_Globals_Vars           ,ONLY: eps0
USE MOD_TimeDisc_Vars          ,ONLY: dt
USE MOD_Mesh_Vars              ,ONLY: N_SurfMesh,SideToElem,nBCSides,N_SurfMesh,offSetElem,BC
USE MOD_DG_Vars                ,ONLY: U_N,N_DG_Mapping
USE MOD_PICDepo_Vars           ,ONLY: PS_N
USE MOD_Particle_Boundary_Vars ,ONLY: N_SurfVDL,PartBound
USE MOD_ProlongToFace          ,ONLY: ProlongToFace_Side
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ElemID,SideID,ilocSide,Nloc,iPartBound,p,q
REAL,ALLOCATABLE :: jface(:,:,:)
REAL             :: coeff
LOGICAL          :: UpdatePhiF
!===================================================================================================================================
! 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!     interpolate the vector field E = (/Ex, Ey, Ez/) to the boundary face
DO SideID=1,nBCSides
  ! Get the local element index
    ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  ! Get local polynomial degree of the element
  Nloc   = N_DG_Mapping(2,ElemID+offSetElem)
  ! Get particle boundary index
  iPartBound = PartBound%MapToPartBC(BC(SideID))

  ! Skip sides that are not a VDL boundary (these sides are still in the list of sides)
  IF(PartBound%ThicknessVDL(iPartBound).GT.0.0)THEN

    ! Allocate jface depending on the local polynomial degree
    ALLOCATE(jface(1:3,0:Nloc,0:Nloc))
    ! Get local side index
    ilocSide = SideToElem(S2E_LOC_SIDE_ID,SideID)

    IF(UpdatePhiF)THEN
      ! Calculate coefficient
      coeff = (PartBound%ThicknessVDL(iPartBound)*dt)/(eps0*PartBound%PermittivityVDL(iPartBound))

      ! Update PhiF:  PhiF_From_Currents: PhiF^n - PhiF^n-1 = -(d*dt)/(eps0*epsR)*(jp+jD)
      U_N(ElemID)%PhiF(1:3,:,:,:) = U_N(ElemID)%PhiF(1:3,:,:,:) &
                                  - coeff * ( PS_N(ElemID)%PartSource(1:3,:,:,:) + U_N(ElemID)%Dt(1:3,:,:,:) )
    END IF ! UpdatePhiF
    !WRITE (*,*) "XXXX: jp + jD =", MAXVAL(ABS(PS_N(ElemID)%PartSource(1:3,:,:,:))), MAXVAL(ABS(U_N(ElemID)%Dt(1:3,:,:,:)))

    ! Prolong-to-face depending on orientation in reference element
    CALL ProlongToFace_Side(3, Nloc, ilocSide, 0, U_N(ElemID)%PhiF(1:3,0:Nloc,0:Nloc,0:Nloc), jface(1:3,0:Nloc,0:Nloc))

    !IF(UpdatePhiF) WRITE (*,*) "time   , N_SurfVDL(SideID)%U(9,0,0) =", time, N_SurfVDL(SideID)%U(9,0,0)

    ! 2.) Apply the normal vector to get the flux over the boundary face
    DO q=0,Nloc
      DO p=0,Nloc
        ASSOCIATE( normal => N_SurfMesh(SideID)%NormVec(1:3,p,q), j => jface(1:3,p,q), PhiF => N_SurfVDL(SideID)%U(9,p,q))
          ! j_normal = <j,normal> with normal inverted
          !jface(1,p,q) = DOT_PRODUCT(j,-normal)
          PhiF = DOT_PRODUCT(j,-normal)
        END ASSOCIATE
      END DO ! p
    END DO ! q

    ! 3.) Get E from Phi_F
    DO q=0,Nloc
      DO p=0,Nloc
        ASSOCIATE(      E => N_SurfVDL(SideID)%U(10:12,p,q)     ,&
                   normal => N_SurfMesh(SideID)%NormVec(1:3,p,q),&
                        j => jface(1,p,q)                       ,&
                     PhiF => N_SurfVDL(SideID)%U(9,p,q)         )
          ! Reconstruct Phi_F from the current density and the (uncorrected) electric displacement fields in each element
          ! PhiF_From_Currents: PhiF^n - PhiF^n-1 = -(d*dt)/(eps0*epsR)*(jp+jD)
          !PhiF = j

          ! Reconstruct E from Phi_Max via E = Phi/d
          E =  PhiF/PartBound%ThicknessVDL(iPartBound)*normal

        END ASSOCIATE
      END DO ! p
    END DO ! q
    !WRITE (*,*) "N_SurfVDL(SideID)%U(9,:,:) =", N_SurfVDL(SideID)%U(9,:,:)
    !IF(UpdatePhiF)THEN WRITE (*,*) "time+dt, N_SurfVDL(SideID)%U(9,0,0) =", time+dt, N_SurfVDL(SideID)%U(9,0,0)
    DEALLOCATE(jface)
  END IF ! PartBound%ThicknessVDL(iPartBound).GT.0.0
END DO ! SideID=1,nBCSides

END SUBROUTINE CalculatePhiAndEFieldFromCurrentsVDL
#endif /*defined(PARTICLES)*/
#endif /*(USE_HDG && (PP_nVar==1))*/


!===================================================================================================================================
!> During restart, recalculate the gradient of the HDG solution
!===================================================================================================================================
SUBROUTINE RecomputeEFieldHDG()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_Elem_Mat      ,ONLY: PostProcessGradientHDG
USE MOD_Basis         ,ONLY: getSPDInverse
#if USE_MPI
USE MOD_MPI_Vars
#endif /*USE_MPI*/
#if (PP_nVar==1)
!USE MOD_Equation_Vars ,ONLY: E
#elif (PP_nVar==3)
USE MOD_Equation_Vars ,ONLY: B
#else
USE MOD_Equation_Vars ,ONLY: B, E
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol(PP_N),PP_nElems)
#endif
!===================================================================================================================================

#if (PP_nVar==1)
  CALL PostProcessGradientHDG()
#elif (PP_nVar==3)
  ! NOT IMPLEMENTED
  CALL abort(__STAMP__,"ERROR: Functionality not implemented!")
  ! DO iVar=1, PP_nVar
  !   CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  ! END DO
  ! DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !   r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
  !   B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
  !   B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
  !   B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  ! END DO; END DO; END DO !i,j,k
#else
  ! NOT IMPLEMENTED
  CALL abort(__STAMP__,"ERROR: Functionality not implemented!")
  ! DO iVar=1, 3
  !   CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  ! END DO
  ! DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !   r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
  !   B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
  !   B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
  !   B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  ! END DO; END DO; END DO !i,j,k
  ! CALL PostProcessGradient(U_out(4,:,:),lambda(4,:,:),E)
#endif
END SUBROUTINE RecomputeEFieldHDG


!===================================================================================================================================
!> Finalizes variables necessary for HDG subroutines
!===================================================================================================================================
SUBROUTINE FinalizeHDG()
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars
USE MOD_Interpolation_Vars ,ONLY: NMax
#if USE_PETSC
USE petsc
USE MOD_HDG_Vars_PETSc
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem,nSides,SideToNonUniqueGlobalSide,N_SurfMesh
USE MOD_Mesh_Tools         ,ONLY: LambdaSideToMaster,GetMasteriLocSides
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY: SurfExchange
#endif /*USE_MPI*/
USE MOD_Interpolation_Vars ,ONLY: N_Inter
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_PETSC
PetscErrorCode       :: ierr
#endif
#if USE_LOADBALANCE
INTEGER             :: NonUniqueGlobalSideID
INTEGER             :: iSide
#endif /*USE_LOADBALANCE*/
#if USE_MPI
INTEGER             :: iBC
#endif /*USE_MPI*/
INTEGER           :: Nloc
!===================================================================================================================================
HDGInitIsDone = .FALSE.
#if USE_PETSC
PetscCallA(KSPDestroy(PETScSolver,ierr))
PetscCallA(MatDestroy(PETScSystemMatrix,ierr))
PetscCallA(VecDestroy(PETScSolution,ierr))
PetscCallA(VecDestroy(PETScSolutionLocal,ierr))
PetscCallA(VecScatterDestroy(PETScScatter,ierr))
PetscCallA(VecDestroy(PETScRHS,ierr))
PetscCallA(PetscFinalize(ierr))
SDEALLOCATE(SmallMortarType)
SDEALLOCATE(OffsetGlobalPETScDOF)
#endif
!SDEALLOCATE(NonlinVolumeFac)
SDEALLOCATE(DirichletBC)
SDEALLOCATE(NeumannBC)
SDEALLOCATE(HDG_Vol_N)
SDEALLOCATE(qn_face_MagStat)
!SDEALLOCATE(delta)
!SDEALLOCATE(LL_minus)
!SDEALLOCATE(LL_plus)
!SDEALLOCATE(Lomega_m)
!SDEALLOCATE(Lomega_p)
!SDEALLOCATE(Domega)
!SDEALLOCATE(InvDhat)
!SDEALLOCATE(wGP_vol)
!SDEALLOCATE(JwGP_vol)
!SDEALLOCATE(Ehat)
!SDEALLOCATE(Smat)
SDEALLOCATE(Tau)
!SDEALLOCATE(RHS_vol)
!SDEALLOCATE(Precond)
!SDEALLOCATE(InvPrecondDiag)
SDEALLOCATE(MaskedSide)
SDEALLOCATE(SmallMortarInfo)
!SDEALLOCATE(IntMatMortar)

#if USE_LOADBALANCE
! MPIRoot does not deallocate during load balance because this process sends the info to the other processors via the respective
! sub-communicators after the new domain decomposition is performed
IF(MPIRoot)THEN
  IF(.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))THEN
#endif /*USE_LOADBALANCE*/
    SDEALLOCATE(FPC%Charge)
    SDEALLOCATE(FPC%Voltage)
#if USE_LOADBALANCE
  END IF ! .NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))
ELSE
  SDEALLOCATE(FPC%Charge)
  SDEALLOCATE(FPC%Voltage)
END IF ! MPIRoot
#endif /*USE_LOADBALANCE*/
SDEALLOCATE(ConductorBC)
SDEALLOCATE(FPC%Group)
SDEALLOCATE(FPC%BCState)
SDEALLOCATE(FPC%VoltageProc)
SDEALLOCATE(FPC%ChargeProc)
#if USE_MPI
DO iBC = 1, FPC%nUniqueFPCBounds
  IF(FPC%COMM(iBC)%UNICATOR.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(FPC%COMM(iBC)%UNICATOR,iERROR)
END DO
SDEALLOCATE(FPC%COMM)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
! MPIRoot does not deallocate during load balance because this process sends the info to the other processors via the respective
! sub-communicators after the new domain decomposition is performed
IF(MPIRoot)THEN
  IF(.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))THEN
#endif /*USE_LOADBALANCE*/
    SDEALLOCATE(EPC%Voltage)
#if USE_LOADBALANCE
  END IF ! .NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))
ELSE
  SDEALLOCATE(EPC%Voltage)
END IF ! MPIRoot
#endif /*USE_LOADBALANCE*/
SDEALLOCATE(EPC%Charge)
SDEALLOCATE(EPC%Resistance)
SDEALLOCATE(EPC%Group)
SDEALLOCATE(EPC%BCState)
SDEALLOCATE(EPC%VoltageProc)
SDEALLOCATE(EPC%ChargeProc)
#if USE_MPI
DO iBC = 1, EPC%nUniqueEPCBounds
  IF(EPC%COMM(iBC)%UNICATOR.NE.MPI_COMM_NULL)   CALL MPI_COMM_FREE(EPC%COMM(iBC)%UNICATOR,iERROR)
END DO
SDEALLOCATE(EPC%COMM)
#if defined(PARTICLES)
IF(BiasVoltage%COMM%UNICATOR.NE.MPI_COMM_NULL)  CALL MPI_COMM_FREE(BiasVoltage%COMM%UNICATOR,iERROR)
IF(CPPCOMM%UNICATOR.NE.MPI_COMM_NULL)           CALL MPI_COMM_FREE(CPPCOMM%UNICATOR,iERROR)
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI*/

#if USE_LOADBALANCE
IF(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))THEN
  ! Store lambda solution on global non-unique array for MPI communication
  ASSOCIATE( firstSide => ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1) + 1       ,&
             lastSide  => ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElem    + nElems) )
    ALLOCATE(lambdaLB(PP_nVar,nGP_face(NMax)+1,firstSide:lastSide)) ! +1 comes from the NSideMin info that is sent additionally
    ! TODO NSideMin - What?
    lambdaLB=0.
  END ASSOCIATE
  IF(nProcessors.GT.1) CALL GetMasteriLocSides()
  DO iSide = 1, nSides
    NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,iSide)

    CALL LambdaSideToMaster(1,iSide,lambdaLB(:,:,NonUniqueGlobalSideID),N_SurfMesh(iSide)%NSide)
    ! Check if the same global unique side is encountered twice and store both global non-unique side IDs in the array
    ! SideToNonUniqueGlobalSide(1:2,iSide)
    IF(SideToNonUniqueGlobalSide(2,iSide).NE.-1)THEN
      NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(2,iSide)
      CALL LambdaSideToMaster(1,iSide,lambdaLB(:,:,NonUniqueGlobalSideID),N_SurfMesh(iSide)%NSide)
    END IF ! SideToNonUniqueGlobalSide(1,iSide).NE.-1

  END DO ! iSide = 1, nSides
  IF(nProcessors.GT.1) DEALLOCATE(iLocSides)

END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/

!SDEALLOCATE(lambda)
SDEALLOCATE(HDG_Surf_N)
SDEALLOCATE(nGP_vol)
SDEALLOCATE(nGP_face)
#if USE_MPI
SDEALLOCATE(SurfExchange)
#endif /*USE_MPI*/
DO Nloc = 1, NMax
  DEALLOCATE(N_Inter(Nloc)%LL_minus)
  DEALLOCATE(N_Inter(Nloc)%LL_plus)
  DEALLOCATE(N_Inter(Nloc)%Lomega_m)
  DEALLOCATE(N_Inter(Nloc)%Lomega_p)
  DEALLOCATE(N_Inter(Nloc)%Domega)
  DEALLOCATE(N_Inter(Nloc)%wGP_vol)
END DO

END SUBROUTINE FinalizeHDG
#endif /*USE_HDG*/


END MODULE MOD_HDG