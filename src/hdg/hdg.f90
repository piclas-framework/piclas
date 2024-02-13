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
INTERFACE InitHDG
  MODULE PROCEDURE InitHDG
END INTERFACE

!INTERFACE HDG
!  MODULE PROCEDURE HDG
!END INTERFACE

INTERFACE FinalizeHDG
  MODULE PROCEDURE FinalizeHDG
END INTERFACE

PUBLIC :: InitHDG,FinalizeHDG
PUBLIC :: HDG, RestartHDG
PUBLIC :: DefineParametersHDG
#if USE_MPI
PUBLIC :: SynchronizeChargeOnFPC,SynchronizeVoltageOnEPC
#if defined(PARTICLES)
 PUBLIC :: SynchronizeBV
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI */
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> Define parameters for HDG (Hubridized Discontinous Galerkin)
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
CALL prms%CreateIntOption(    'PrecondType'            ,'Preconditioner type\n 0: no preconditioner\n 1: Side-block SPD preconditioner matrix (already Cholesky decomposed)\n 2: Inverse of diagonal preconditioned', '2')
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
USE MOD_Interpolation_Vars    ,ONLY: xGP,wGP,L_minus,L_plus
USE MOD_Basis                 ,ONLY: PolynomialDerivativeMatrix
USE MOD_Interpolation_Vars    ,ONLY: wGP
USE MOD_Elem_Mat              ,ONLY: Elem_Mat,BuildPrecond
USE MOD_ReadInTools           ,ONLY: GETLOGICAL,GETREAL,GETINT
USE MOD_Mesh_Vars             ,ONLY: sJ,nBCSides
USE MOD_Mesh_Vars             ,ONLY: BoundaryType,nSides,BC
USE MOD_Mesh_Vars             ,ONLY: nGlobalMortarSides,nMortarMPISides
USE MOD_Mesh_Vars             ,ONLY: DoSwapMesh
USE MOD_ChangeBasis           ,ONLY: ChangeBasis2D
USE MOD_Basis                 ,ONLY: InitializeVandermonde,LegendreGaussNodesAndWeights,BarycentricWeights
USE MOD_FillMortar_HDG        ,ONLY: InitMortar_HDG
USE MOD_HDG_Vars              ,ONLY: BRNbrOfRegions,ElemToBRRegion,RegionElectronRef
#if defined(PARTICLES)
USE MOD_Part_BR_Elecron_Fluid ,ONLY: UpdateNonlinVolumeFac
USE MOD_Restart_Vars          ,ONLY: DoRestart
#endif /*defined(PARTICLES)*/
#if USE_PETSC
USE PETSc
USE MOD_Mesh_Vars             ,ONLY: nMPISides_YOUR
#if USE_MPI
USE MOD_MPI                   ,ONLY: StartReceiveMPIDataInt,StartSendMPIDataInt,FinishExchangeMPIData
USE MOD_MPI_Vars
#endif /*USE_MPI*/
USE MOD_Mesh_Vars             ,ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars             ,ONLY: firstMortarInnerSide,lastMortarInnerSide
#endif /*USE_PETSC*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,r,iElem,SideID
INTEGER           :: BCType,BCState
REAL              :: D(0:PP_N,0:PP_N)
INTEGER           :: nDirichletBCsidesGlobal
#if USE_PETSC
PetscErrorCode    :: ierr
INTEGER           :: iProc
INTEGER           :: OffsetPETScSideMPI(nProcessors)
INTEGER           :: OffsetPETScSide
INTEGER           :: PETScLocalID
INTEGER           :: MortarSideID,iMortar
INTEGER           :: locSide,nMortarMasterSides,nMortars
!INTEGER           :: nAffectedBlockSides
INTEGER,ALLOCATABLE :: indx(:)
#endif
!===================================================================================================================================
IF(HDGInitIsDone)THEN
   LBWRITE(*,*) "InitHDG already called."
   RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT HDG...'

HDGDisplayConvergence = GETLOGICAL('HDGDisplayConvergence')

nGP_vol  = (PP_N+1)**3
nGP_face = (PP_N+1)**2

HDGSkip = GETINT('HDGSkip')
IF (HDGSkip.GT.0) THEN
  HDGSkipInit = GETINT('HDGSkipInit')
  HDGSkip_t0  = GETREAL('HDGSkip_t0')
ELSE
  HDGSkip=0
END IF

HDGNonLinSolver = -1 ! init

#if USE_PETSC
! initialize PETSc stuff!
PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
#endif

#if defined(PARTICLES)
! BR electron fluid model
IF (BRNbrOfRegions .GT. 0) THEN !Regions only used for Boltzmann Electrons so far -> non-linear HDG-sources!
  HDGNonLinSolver=GETINT('HDGNonLinSolver')

  IF (HDGNonLinSolver.EQ.1) THEN
    NewtonExactSourceDeriv  = GETLOGICAL('NewtonExactSourceDeriv')
    AdaptIterNewton         = GETINT('AdaptIterNewton')
    AdaptIterNewtonOld      = AdaptIterNewton
    NewtonAdaptStartValue   = GETLOGICAL('NewtonAdaptStartValue')
    AdaptIterNewtonToLinear = GETINT('AdaptIterNewtonToLinear')
    IF (NewtonExactSourceDeriv) NewtonAdaptStartValue=.TRUE.
    IF (DoRestart) NewtonAdaptStartValue=.FALSE.
    ALLOCATE(NonlinVolumeFac(nGP_vol,PP_nElems))

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

!CG parameters
#if USE_PETSC
LBWRITE(UNIT_stdOut,'(A)') ' Method for HDG solver: PETSc '
#else
LBWRITE(UNIT_stdOut,'(A)') ' Method for HDG solver: CG '
#endif /*USE_PETSC*/
PrecondType          = GETINT('PrecondType')
epsCG                = GETREAL('epsCG')
OutIterCG            = GETINT('OutIterCG')
useRelativeAbortCrit = GETLOGICAL('useRelativeAbortCrit')
MaxIterCG            = GETINT('MaxIterCG')

ExactLambda          = GETLOGICAL('ExactLambda')

ALLOCATE(MaskedSide(1:nSides))
MaskedSide=0

IF(nGlobalMortarSides.GT.0)THEN !mortar mesh
  IF(nMortarMPISides.GT.0) CALL abort(__STAMP__,&
  "nMortarMPISides >0: HDG mortar MPI implementation relies on big sides having always only master sides (=> nMortarMPISides=0 )")
END IF !mortarMesh

CALL InitMortar_HDG()

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
  CASE(10,11) ! Neumann
    nNeumannBCsides=nNeumannBCsides+1
  CASE(20) ! Conductor: Floating Boundary Condition (FPC)
    nConductorBCsides=nConductorBCsides+1
  CASE DEFAULT ! unknown BCType
    CALL CollectiveStop(__STAMP__,' unknown BC Type in hdg.f90!',IntInfo=BCType)
  END SELECT ! BCType
END DO

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

! Get the global number of Dirichlet boundaries. If there are none, the potential of a single DOF must be set.
#if USE_MPI
  CALL MPI_ALLREDUCE(nDirichletBCsides , nDirichletBCsidesGlobal , 1 , MPI_INTEGER , MPI_MAX , MPI_COMM_PICLAS , IERROR)
#else
  nDirichletBCsidesGlobal = nDirichletBCsides
#endif /*USE_MPI*/
#if USE_PETSC
IF(nDirichletBCsidesGlobal.EQ.0) THEN
#else
IF(MPIroot .AND. (nDirichletBCsidesGlobal.EQ.0)) THEN
#endif
  SetZeroPotentialDOF = .TRUE.
ELSE
  SetZeroPotentialDOF = .FALSE.
END IF

IF(nDirichletBCsides.GT.0)ALLOCATE(DirichletBC(nDirichletBCsides))
IF(nNeumannBCsides  .GT.0)THEN
  ALLOCATE(NeumannBC(nNeumannBCsides))
  ALLOCATE(qn_face(PP_nVar, nGP_face,nNeumannBCsides))
END IF
IF(nConductorBCsides.GT.0)ALLOCATE(ConductorBC(nConductorBCsides))
#if (PP_nVar!=1)
  IF(nDirichletBCsides.GT.0)ALLOCATE(qn_face_MagStat(PP_nVar, nGP_face,nDirichletBCsides))
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
  CASE(10,11) !Neumann,
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

#if USE_PETSC
! Create PETSc Mappings
OffsetPETScSide=0
#if USE_MPI
! Count all Mortar slave sides and remove them from PETSc vector
! TODO How to compute those
nMortarMasterSides = 0
DO SideID=1,nSides
  IF(SmallMortarInfo(SideID).EQ.1) THEN
    nMortarMasterSides = nMortarMasterSides + 1
  END IF
END DO
nPETScUniqueSides = nSides-nDirichletBCSides-nMPISides_YOUR-nMortarMasterSides-nConductorBCsides
CALL MPI_ALLGATHER(nPETScUniqueSides,1,MPI_INTEGER,OffsetPETScSideMPI,1,MPI_INTEGER,MPI_COMM_PICLAS,IERROR)
DO iProc=1, myrank
  OffsetPETScSide = OffsetPETScSide + OffsetPETScSideMPI(iProc)
END DO
nPETScUniqueSidesGlobal = SUM(OffsetPETScSideMPI) + FPC%nUniqueFPCBounds
#endif

ALLOCATE(PETScGlobal(nSides))
ALLOCATE(PETScLocalToSideID(nPETScUniqueSides+nMPISides_YOUR))
PETScGlobal=-1
PETScLocalToSideID=-1
PETScLocalID=0 ! = nSides-nDirichletBCSides
DO SideID=1,nSides!-nMPISides_YOUR
  IF(MaskedSide(SideID).GT.0) CYCLE
  PETScLocalID=PETScLocalID+1
  PETScLocalToSideID(PETScLocalID)=SideID
  PETScGlobal(SideID)=PETScLocalID+OffsetPETScSide-1 ! PETSc arrays start at 0!
END DO
! Set the Global PETSc Sides of small mortar sides equal to the big mortar side
DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide) !small SideID
    PETScGlobal(SideID)=PETScGlobal(MortarSideID)
  END DO !iMortar
END DO
#if USE_MPI
CALL StartReceiveMPIDataInt(1,PETScGlobal,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
CALL StartSendMPIDataInt(   1,PETScGlobal,1,nSides,SendRequest_U,SendID=1) ! Send MINE
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif
#endif

!mappings
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

ALLOCATE(delta(0:PP_N,0:PP_N))
delta=0.
DO i=0,PP_N
  delta(i,i)=1.
END DO !i

ALLOCATE(LL_minus(0:PP_N,0:PP_N))
ALLOCATE(LL_plus( 0:PP_N,0:PP_N))
DO j=0,PP_N
  DO i=0,PP_N
    LL_minus(i,j) = L_minus(i)*L_minus(j)
    LL_plus(i,j)  = L_plus(i)*L_plus(j)
  END DO
END DO

ALLOCATE(Lomega_m(0:PP_N))
ALLOCATE(Lomega_p(0:PP_N))
! Compute a lifting matrix scaled by the Gaussian weights
Lomega_m = - L_minus/wGP
Lomega_p = + L_plus/wGP
ALLOCATE(Domega(0:PP_N,0:PP_N))
! Compute Differentiation matrix D for given Gausspoints (1D)
CALL PolynomialDerivativeMatrix(PP_N,xGP,D)
! Compute a Differentiation mtarix scaled by the Gaussian weigths
DO j=0,PP_N
  DO i=0,PP_N
    Domega(i,j) = wGP(i)/wGP(j)*D(i,j)
  END DO !r
END DO !s

ALLOCATE(InvDhat(nGP_vol,nGP_vol,PP_nElems))
ALLOCATE(wGP_vol(nGP_vol))
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
  wGP_vol(r)=wGP(i)*wGP(j)*wGP(k)
END DO; END DO; END DO !i,j,k

ALLOCATE(JwGP_vol(nGP_vol,PP_nElems))
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    JwGP_vol(r,iElem)=wGP_vol(r)/sJ(i,j,k,iElem) !omega*J
  END DO; END DO; END DO !i,j,k
END DO !iElem


ALLOCATE(Ehat(nGP_face,nGP_vol,6,PP_nElems))
!side matrices
ALLOCATE(Smat(nGP_face,nGP_face,6,6,PP_nElems))

#if USE_PETSC
ALLOCATE(Smat_BC(nGP_face,nGP_face,6,nDirichletBCSides))
Smat_BC = 0.

PetscCallA(MatCreate(PETSC_COMM_WORLD,Smat_petsc,ierr))
PetscCallA(MatSetBlockSize(Smat_petsc,nGP_face,ierr))
PetscCallA(MatSetSizes(Smat_petsc,PETSC_DECIDE,PETSC_DECIDE,nPETScUniqueSidesGlobal*nGP_Face,nPETScUniqueSidesGlobal*nGP_Face,ierr))
PetscCallA(MatSetType(Smat_petsc,MATSBAIJ,ierr)) ! Symmetric sparse (mpi) matrix
!! TODO Set preallocation row wise
!! 1 Big mortar side is affected by 6 + 4*4 = 22 other sides...
!! TODO Does this require communication over all procs? Global number of sides associated with the i-th FPC
!IF(FPC%nFPCBounds.GT.0)THEN
!  ALLOCATE(FPC%GroupGlobal(1:FPC%nFPCBounds))
!  FPC%GroupGlobal(1:FPC%nFPCBounds) = FPC%Group(1:FPC%nFPCBounds,3)
!  ! TODO is this allreduce required?
!  !CALL MPI_ALLREDUCE(FPC%Group(1:FPC%nFPCBounds,3),FPC%GroupGlobal(1:FPC%nFPCBounds), FPC%nFPCBounds, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_PICLAS, IERROR)
!  nAffectedBlockSides = MAXVAL(FPC%GroupGlobal(:))
!  DEALLOCATE(FPC%GroupGlobal)
!  nAffectedBlockSides = MAX(22,nAffectedBlockSides*6)
!ELSE
!  nAffectedBlockSides = 22
!END IF ! FPC%nFPCBounds
!!IPWRITE(UNIT_StdOut,*) "nAffectedBlockSides =", nAffectedBlockSides
PetscCallA(MatSEQSBAIJSetPreallocation(Smat_petsc,nGP_face,22,PETSC_NULL_INTEGER,ierr))
PetscCallA(MatMPISBAIJSetPreallocation(Smat_petsc,nGP_face,22,PETSC_NULL_INTEGER,22-1,PETSC_NULL_INTEGER,ierr))
PetscCallA(MatZeroEntries(Smat_petsc,ierr))
#endif

!stabilization parameter
ALLOCATE(Tau(PP_nElems))
DO iElem=1,PP_nElems
  Tau(iElem)=2./((SUM(JwGP_vol(:,iElem)))**(1./3.))  !1/h ~ 1/vol^(1/3) (volref=8)
END DO !iElem

IF(.NOT.DoSwapMesh)THEN ! can take very long, not needed for swap mesh run as only the state file is converted
  CALL Elem_Mat(0_8) ! takes iter=0 (kind=8)
END IF

#if USE_PETSC
PetscCallA(KSPCreate(PETSC_COMM_WORLD,ksp,ierr))
PetscCallA(KSPSetOperators(ksp,Smat_petsc,Smat_petsc,ierr))

IF(PrecondType.GE.10) THEN
  PetscCallA(KSPSetType(ksp,KSPPREONLY,ierr)) ! Exact solver
ELSE
  PetscCallA(KSPSetType(ksp,KSPCG,ierr)) ! CG solver for sparse symmetric positive definite matrix

  PetscCallA(KSPSetInitialGuessNonzero(ksp,PETSC_TRUE, ierr))

  PetscCallA(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED, ierr))
  PetscCallA(KSPSetTolerances(ksp,1.E-20,epsCG,PETSC_DEFAULT_REAL,MaxIterCG,ierr))
END IF
#endif

CALL BuildPrecond()

ALLOCATE(lambda(PP_nVar,nGP_face,nSides))
lambda=0.
ALLOCATE(RHS_vol(PP_nVar, nGP_vol,PP_nElems))
RHS_vol=0.

#if USE_PETSC
! allocate RHS & lambda vectors
PetscCallA(VecCreate(PETSC_COMM_WORLD,lambda_petsc,ierr))
PetscCallA(VecSetBlockSize(lambda_petsc,nGP_face,ierr))
PetscCallA(VecSetSizes(lambda_petsc,PETSC_DECIDE,nPETScUniqueSidesGlobal*nGP_Face,ierr))
PetscCallA(VecSetType(lambda_petsc,VECSTANDARD,ierr))
PetscCallA(VecSetUp(lambda_petsc,ierr))
PetscCallA(VecDuplicate(lambda_petsc,RHS_petsc,ierr))

! Create scatter context to access local values from global petsc vector
PetscCallA(VecCreateSeq(PETSC_COMM_SELF,nPETScUniqueSides*nGP_face,lambda_local_petsc,ierr))
PetscCallA(ISCreateStride(PETSC_COMM_SELF,nPETScUniqueSides*nGP_face,0,1,idx_local_petsc,ierr))
PetscCallA(ISCreateBlock(PETSC_COMM_WORLD,nGP_face,nPETScUniqueSides,PETScGlobal(PETScLocalToSideID(1:nPETScUniqueSides)),PETSC_COPY_VALUES,idx_global_petsc,ierr))
PetscCallA(VecScatterCreate(lambda_petsc,idx_global_petsc,lambda_local_petsc,idx_local_petsc,scatter_petsc,ierr))

IF(UseFPC)THEN
  PetscCallA(VecCreateSeq(PETSC_COMM_SELF,nGP_face*FPC%nUniqueFPCBounds,lambda_local_conductors_petsc,ierr))
  PetscCallA(ISCreateStride(PETSC_COMM_SELF,nGP_face*FPC%nUniqueFPCBounds,0,1,idx_local_conductors_petsc,ierr))
  ALLOCATE(indx(FPC%nUniqueFPCBounds))
  DO i=1,FPC%nUniqueFPCBounds
    indx(i) = nPETScUniqueSidesGlobal-FPC%nUniqueFPCBounds+i-1
  END DO
  PetscCallA(ISCreateBlock(PETSC_COMM_WORLD,nGP_face,FPC%nUniqueFPCBounds,indx,PETSC_COPY_VALUES,idx_global_conductors_petsc,ierr))
  DEALLOCATE(indx)
  PetscCallA(VecScatterCreate(lambda_petsc,idx_global_conductors_petsc,lambda_local_conductors_petsc,idx_local_conductors_petsc,scatter_conductors_petsc,ierr))
END IF
#endif

HDGInitIsDone = .TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT HDG DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitHDG


!===================================================================================================================================
!> Create containers and communicators for each floating boundary condition where impacting charges are accumulated.
!>
!> 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
!>     contribute to the total floating boundary condition. If yes, then this processor is part of the communicator
!> 2.) Create Mapping from floating boundary condition BC index to field BC index
!> 3.) Create Mapping from field BC index to floating boundary condition BC index
!> 4.0) Check if field BC is on current proc (or MPI root)
!> 4.1.) Each processor loops over all of his elements
!> 4.2.) Loop over all compute-node elements (every processor loops over all of these elements)
!> 5.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitFPC()
! MODULES
USE MOD_Globals  ! ,ONLY: MPIRoot,iError,myrank,UNIT_stdOut,MPI_COMM_PICLAS
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nBCs,BoundaryType
USE MOD_Analyze_Vars       ,ONLY: DoFieldAnalyze
USE MOD_HDG_Vars           ,ONLY: UseFPC,FPC
USE MOD_Mesh_Vars          ,ONLY: nBCSides,BC
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI && defined(PARTICLES)
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
USE MOD_Globals            ,ONLY: ElementOnProc
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,BoundsOfElem_Shared,SideInfo_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Vars  ,ONLY: halo_eps
USE MOD_Mesh_Vars          ,ONLY: nElems, offsetElem
#endif /*USE_MPI && defined(PARTICLES)*/
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER  :: BCTypeFPC = 20
INTEGER             :: BCType,BCState,iUniqueFPCBC
INTEGER             :: SideID,iBC
#if USE_MPI
INTEGER             :: color,WithSides
LOGICAL,ALLOCATABLE :: BConProc(:)
#if defined(PARTICLES)
INTEGER             :: iElem,iCNElem
REAL                :: iElemCenter(1:3),iGlobElemCenter(1:3)
REAL                :: iElemRadius,iGlobElemRadius
INTEGER             :: iGlobElem,BCIndex,iSide
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI*/
CHARACTER(5)        :: hilf,hilf2
!===================================================================================================================================

! Get global number of FPC boundaries in [1:nBCs], they might belong to the same group (will be reduced to "nUniqueFPCBounds" below)
! FPC boundaries with the same BCState will be in the same group (electrically connected)
UseFPC = .FALSE.
FPC%nFPCBounds = 0
FPC%nUniqueFPCBounds = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
  FPC%nFPCBounds=FPC%nFPCBounds+1
  IF(BCState.LE.0) CALL CollectiveStop(__STAMP__,' BCState for FPC must be >0! BCState=',IntInfo=BCState)
END DO

IF(FPC%nFPCBounds.EQ.0) RETURN ! Already determined in HDG initialization

UseFPC = .TRUE.

ALLOCATE(FPC%Group(1:FPC%nFPCBounds,3))
FPC%Group = 0 ! Initialize

! 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
! contribute to the total floating boundary condition. If yes, then this processor is part of the communicator
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
  WRITE(UNIT=hilf,FMT='(I0)') BCState
  WRITE(UNIT=hilf2,FMT='(I0)') FPC%nFPCBounds
  IF(BCState.GT.FPC%nFPCBounds) CALL CollectiveStop(__STAMP__,&
      'BCState='//TRIM(hilf)//' must be smaller or equal than the total number of '//TRIM(hilf2)//' FPC boundaries!')
  FPC%Group(BCState,1) = FPC%Group(BCState,1) + 1
  IF(FPC%Group(BCState,1).EQ.1) THEN
    FPC%nUniqueFPCBounds = FPC%nUniqueFPCBounds +1 ! Only count once
    FPC%Group(BCState,2) = FPC%nUniqueFPCBounds
  END IF
END DO

! Automatically activate surface model analyze flag
DoFieldAnalyze = .TRUE.

! 2.) Create Mapping from floating boundary condition BC index to BCState
ALLOCATE(FPC%BCState(FPC%nUniqueFPCBounds))
FPC%BCState = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE
  BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
  iUniqueFPCBC = FPC%Group(BCState,2)
  FPC%BCState(iUniqueFPCBC) = BCState
END DO

! Allocate the containers
! This container is not deallocated for MPIRoot when performing load balance (only root needs this info to write it to .csv)
IF(.NOT.ALLOCATED(FPC%Voltage))THEN
  ALLOCATE(FPC%Voltage(1:FPC%nUniqueFPCBounds))
  FPC%Voltage = 0.
END IF ! .NOT.ALLOCATED(FPC%Voltage)
ALLOCATE(FPC%VoltageProc(1:FPC%nUniqueFPCBounds))
FPC%VoltageProc = 0.
! This container is not deallocated for MPIRoot when performing load balance as this process updates the information on the new
! sub-communicator processes during load balance
IF(.NOT.ALLOCATED(FPC%Charge))THEN
  ALLOCATE(FPC%Charge(1:FPC%nUniqueFPCBounds))
  FPC%Charge = 0.
END IF ! .NOT.ALLOCATED(FPC%Charge)
ALLOCATE(FPC%ChargeProc(1:FPC%nUniqueFPCBounds))
FPC%ChargeProc = 0.

! Set initial value depending on IniExactFunc
SELECT CASE (IniExactFunc)
CASE(800,900) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite
  ! Set initial value
  FPC%Charge(1)  = 1.0e-2*(GEO%ymaxglob-GEO%yminglob)*(GEO%zmaxglob-GEO%zminglob) ! C/m2
  FPC%Voltage(1) = 1.1293922903231239 ! V
CASE(801) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite: 2D case
  ! Set initial value
  FPC%Charge(1)  = 5e-11*(1.0 - 0.05) ! C/m2
  FPC%Charge(2)  = 5e-11*(1.0 - 0.15) ! C/m2
  FPC%Charge(3)  = 5e-11*(1.0 - 0.25) ! C/m2
  FPC%Charge(4)  = 5e-11*(1.0 - 0.35) ! C/m2
  FPC%Charge(5)  = 5e-11*(1.0 - 0.45) ! C/m2
  FPC%Charge(6)  = 5e-11*(1.0 - 0.55) ! C/m2
  FPC%Charge(7)  = 5e-11*(1.0 - 0.65) ! C/m2
  FPC%Charge(8)  = 5e-11*(1.0 - 0.75) ! C/m2
  FPC%Charge(9)  = 5e-11*(1.0 - 0.85) ! C/m2
  FPC%Charge(10) = 5e-11*(1.0 - 0.95) ! C/m2
END SELECT

!! 3.) Create Mapping from field BC index to floating boundary condition BC index
!ALLOCATE(FPC%BCIDToFPCBCID(nBCs))
!FPC%BCIDToFPCBCID = -1
!DO iFPCBC = 1, FPC%NBoundaries
!  iBC = EDC%FieldBoundaries(iEDCBC)
!  EDC%BCIDToEDCBCID(iBC) = iEDCBC
!END DO ! iEDCBC = 1, EDC%NBoundaries

! Get processor-local number of FPC sides associated with each i-th FPC boundary
! Check local sides
DO SideID=1,nBCSides
  iBC    = BC(SideID)
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iFPC
  FPC%Group(BCState,3) = FPC%Group(BCState,3) + 1
END DO ! SideID=1,nBCSides

#if USE_MPI
! 4.0) Check if field BC is on current proc (or MPI root)
ALLOCATE(BConProc(FPC%nUniqueFPCBounds))
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE

  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
    BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iFPC
    iUniqueFPCBC = FPC%Group(BCState,2)
    BConProc(iUniqueFPCBC) = .TRUE.
  END DO ! SideID=1,nBCSides

#if defined(PARTICLES)
  ! Check if all FPCs have already been found
  IF(.NOT.(ALL(BConProc)))THEN
    ! Particles might impact the FPC on another proc/node. Therefore check if a particle can travel from a local element to an
    ! element that has at least one side, which is an FPC
    ! 4.1.) Each processor loops over all of his elements
    iElemLoop: DO iElem = 1+offsetElem, nElems+offsetElem

      iElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,2,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
      iElemRadius = VECNORM ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem),&
                                BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem),&
                                BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

      ! 4.2.) Loop over all compute-node elements (every processor loops over all of these elements)
      ! Loop ALL compute-node elements (use global element index)
      iCNElemLoop: DO iCNElem = 1,nComputeNodeTotalElems
        iGlobElem = GetGlobalElemID(iCNElem)

        ! Skip my own elements as they have already been tested when the local sides are checked
        IF(ElementOnProc(iGlobElem)) CYCLE iCNElemLoop

        ! Check if one of the six sides of the compute-node element is a FPC
        ! Note that iSide is in the range of 1:nNonUniqueGlobalSides
        DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
          ! Get BC index of the global side index
          BCIndex = SideInfo_Shared(SIDE_BCID,iSide)
          ! Only check BC sides with BC index > 0
          IF(BCIndex.GT.0)THEN
            ! Get boundary type
            BCType = BoundaryType(BCIndex,BC_TYPE)
            ! Check if FPC has been found
            IF(BCType.EQ.BCTypeFPC)THEN

              ! Check if the BC can be reached
              iGlobElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,2,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,3,iGlobElem)) /) / 2.
              iGlobElemRadius = VECNORM ((/ BoundsOfElem_Shared(2,1,iGlobElem)-BoundsOfElem_Shared(1,1,iGlobElem),&
                                            BoundsOfElem_Shared(2,2,iGlobElem)-BoundsOfElem_Shared(1,2,iGlobElem),&
                                            BoundsOfElem_Shared(2,3,iGlobElem)-BoundsOfElem_Shared(1,3,iGlobElem) /) / 2.)

              ! check if compute-node element "iGlobElem" is within halo_eps of processor-local element "iElem"
              IF (VECNORM( iElemCenter(1:3) - iGlobElemCenter(1:3) ) .LE. ( halo_eps + iElemRadius + iGlobElemRadius ) )THEN
                BCState = BoundaryType(BCIndex,BC_STATE) ! BCState corresponds to iFPC
                IF(BCState.LT.1) CALL abort(__STAMP__,'BCState cannot be <1',IntInfoOpt=BCState)
                iUniqueFPCBC = FPC%Group(BCState,2)
                ! Flag the i-th FPC
                BConProc(iUniqueFPCBC) = .TRUE.
                ! Check if all FPCs have been found -> exit complete loop
                IF(ALL(BConProc)) EXIT iElemLoop
                ! Go to next element
                CYCLE iCNElemLoop
              END IF ! VECNORM( ...
            END IF ! BCType.EQ.BCTypeFPC
          END IF ! BCIndex.GT.0
        END DO ! iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
      END DO iCNElemLoop ! iCNElem = 1,nComputeNodeTotalElems
    END DO iElemLoop ! iElem = 1, nElems
  END IF ! .NOT.(ALL(BConProc))
#endif /*defined(PARTICLES)*/

END IF ! MPIRoot

! 5.) Create MPI sub-communicators
ALLOCATE(FPC%COMM(FPC%nUniqueFPCBounds))
DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  ! create new communicator
  color = MERGE(iUniqueFPCBC, MPI_UNDEFINED, BConProc(iUniqueFPCBC))

  ! set communicator id
  FPC%COMM(iUniqueFPCBC)%ID=iUniqueFPCBC

  ! create new emission communicator for floating boundary condition communication. Pass MPI_INFO_NULL as rank to follow the original ordering
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, MPI_INFO_NULL, FPC%COMM(iUniqueFPCBC)%UNICATOR, iError)

  ! Find my rank on the shared communicator, comm size and proc name
  IF(BConProc(iUniqueFPCBC))THEN
    CALL MPI_COMM_RANK(FPC%COMM(iUniqueFPCBC)%UNICATOR, FPC%COMM(iUniqueFPCBC)%MyRank, iError)
    CALL MPI_COMM_SIZE(FPC%COMM(iUniqueFPCBC)%UNICATOR, FPC%COMM(iUniqueFPCBC)%nProcs, iError)

    ! inform about size of emission communicator
    IF (FPC%COMM(iUniqueFPCBC)%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
      IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
          WRITE(UNIT_StdOut,'(A,I0,A,I0,A,I0)') ' Floating boundary condition: Emission-Communicator ',iUniqueFPCBC,' on ',&
              FPC%COMM(iUniqueFPCBC)%nProcs,' procs for BCState ',FPC%BCState(iUniqueFPCBC)
    END IF
  END IF ! BConProc(iUniqueFPCBC)
END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
DEALLOCATE(BConProc)

! Get the number of procs that actually have a local BC side that is an FPC (required for voltage output to .csv)
! Procs might have zero FPC sides but are in the group because 1.) MPIRoot or 2.) the FPC is in the halo region
! Because only the MPI root process writes the .csv data, the information regarding the voltage on each FPC must be
! communicated with this process even though it might not be connected to each FPC boundary
DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR, nProcsWithSides => FPC%COMM(iUniqueFPCBC)%nProcsWithSides )
    IF(FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
      ! Check if the current processor is actually connected to the FPC via a BC side
      IF(FPC%Group(FPC%BCState(iUniqueFPCBC),3).EQ.0)THEN
        WithSides = 0
      ELSE
        WithSides = 1
      END IF ! FPC%Group(FPC%BCState(iUniqueFPCBC),3).EQ.0
      ! Calculate the sum across the sub-communicator. Only the MPI root process needs this information
      IF(MPIRoot)THEN
        CALL MPI_REDUCE(WithSides, nProcsWithSides, 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, iError)
        ! Sanity check
        IF(nProcsWithSides.EQ.0) CALL abort(__STAMP__,'Found FPC with no processors connected to it')
      ELSE
        CALL MPI_REDUCE(WithSides, 0              , 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, IError)
      END IF ! MPIRoot
    END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
  END ASSOCIATE
END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
#endif /*USE_MPI*/

! When restarting, load the deposited charge on each FPC from the .h5 state file
CALL ReadFPCDataFromH5()

END SUBROUTINE InitFPC


!===================================================================================================================================
!> Create containers and communicators for each electric potential boundary condition where impacting charges are removed and
!> subsequently an electric potential is created.
!>
!> 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
!>     contribute to the total electric potential boundary condition. If yes, then this processor is part of the communicator
!> 2.) Create Mapping from electric potential boundary condition BC index to field BC index
!> 3.) Check if field BC is on current proc (or MPI root)
!> 3.1) Each processor loops over all of his elements
!> 3.2) Loop over all compute-node elements (every processor loops over all of these elements)
!> 4.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitEPC()
! MODULES
USE MOD_Globals  ! ,ONLY: MPIRoot,iError,myrank,UNIT_stdOut,MPI_COMM_PICLAS
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nBCs,BoundaryType
USE MOD_Analyze_Vars       ,ONLY: DoFieldAnalyze
USE MOD_HDG_Vars           ,ONLY: UseEPC,EPC
USE MOD_Mesh_Vars          ,ONLY: nBCSides,BC
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI && defined(PARTICLES)
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
USE MOD_Globals            ,ONLY: ElementOnProc
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,BoundsOfElem_Shared,SideInfo_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Vars  ,ONLY: halo_eps
USE MOD_Mesh_Vars          ,ONLY: nElems, offsetElem
#endif /*USE_MPI && defined(PARTICLES)*/
USE MOD_ReadInTools        ,ONLY: GETREALARRAY
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER  :: BCTypeEPC = 8
INTEGER             :: BCType,BCState,iUniqueEPCBC
INTEGER             :: SideID,iBC
#if USE_MPI
INTEGER             :: color,WithSides
LOGICAL,ALLOCATABLE :: BConProc(:)
#if defined(PARTICLES)
INTEGER             :: iElem,iCNElem
REAL                :: iElemCenter(1:3),iGlobElemCenter(1:3)
REAL                :: iElemRadius,iGlobElemRadius
INTEGER             :: iGlobElem,BCIndex,iSide
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI*/
CHARACTER(5)        :: hilf,hilf2
!===================================================================================================================================

! Get global number of EPC boundaries in [1:nBCs], they might belong to the same group (will be reduced to "nUniqueEPCBounds" below)
! EPC boundaries with the same BCState will be in the same group (electrically connected)
UseEPC = .FALSE.
EPC%nEPCBounds = 0
EPC%nUniqueEPCBounds = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
  EPC%nEPCBounds=EPC%nEPCBounds+1
  IF(BCState.LE.0) CALL CollectiveStop(__STAMP__,' BCState for EPC must be >0! BCState=',IntInfo=BCState)
END DO

IF(EPC%nEPCBounds.EQ.0) RETURN ! Already determined in HDG initialization

UseEPC = .TRUE.

ALLOCATE(EPC%Group(1:EPC%nEPCBounds,3))
EPC%Group = 0 ! Initialize

! 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
! contribute to the total electric potential boundary condition. If yes, then this processor is part of the communicator
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
  WRITE(UNIT=hilf,FMT='(I0)') BCState
  WRITE(UNIT=hilf2,FMT='(I0)') EPC%nEPCBounds
  IF(BCState.GT.EPC%nEPCBounds) CALL CollectiveStop(__STAMP__,&
      'BCState='//TRIM(hilf)//' must be smaller or equal than the total number of '//TRIM(hilf2)//' EPC boundaries!')
  EPC%Group(BCState,1) = EPC%Group(BCState,1) + 1
  IF(EPC%Group(BCState,1).EQ.1) THEN
    EPC%nUniqueEPCBounds = EPC%nUniqueEPCBounds +1 ! Only count once
    EPC%Group(BCState,2) = EPC%nUniqueEPCBounds
  END IF
END DO

! Automatically activate surface model analyze flag
DoFieldAnalyze = .TRUE.

! Read resistances for each unique EPC
EPC%Resistance  = GETREALARRAY('EPC-Resistance',EPC%nUniqueEPCBounds)

! 2.) Create Mapping from electric potential boundary condition BC index to BCState
ALLOCATE(EPC%BCState(EPC%nUniqueEPCBounds))
EPC%BCState = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE
  BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
  iUniqueEPCBC = EPC%Group(BCState,2)
  EPC%BCState(iUniqueEPCBC) = BCState
END DO

! Allocate the containers
! This container is not deallocated for MPIRoot when performing load balance (only root needs this info to write it to .csv)
IF(.NOT.ALLOCATED(EPC%Voltage))THEN
  ALLOCATE(EPC%Voltage(1:EPC%nUniqueEPCBounds))
  EPC%Voltage = 0.
END IF ! .NOT.ALLOCATED(EPC%Voltage)
ALLOCATE(EPC%VoltageProc(1:EPC%nUniqueEPCBounds))
EPC%VoltageProc = 0.
! This container is not deallocated for MPIRoot when performing load balance as this process updates the information on the new
! sub-communicator processes during load balance
IF(.NOT.ALLOCATED(EPC%Charge))THEN
  ALLOCATE(EPC%Charge(1:EPC%nUniqueEPCBounds))
  EPC%Charge = 0.
END IF ! .NOT.ALLOCATED(EPC%Charge)
ALLOCATE(EPC%ChargeProc(1:EPC%nUniqueEPCBounds))
EPC%ChargeProc = 0.

! Get processor-local number of EPC sides associated with each i-th EPC boundary
! Check local sides
DO SideID=1,nBCSides
  iBC    = BC(SideID)
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iEPC
  EPC%Group(BCState,3) = EPC%Group(BCState,3) + 1
END DO ! SideID=1,nBCSides

#if USE_MPI
! 3) Check if field BC is on current proc (or MPI root)
ALLOCATE(BConProc(EPC%nUniqueEPCBounds))
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE

  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
    BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iEPC
    iUniqueEPCBC = EPC%Group(BCState,2)
    BConProc(iUniqueEPCBC) = .TRUE.
  END DO ! SideID=1,nBCSides

#if defined(PARTICLES)
  ! Check if all FPCs have already been found
  IF(.NOT.(ALL(BConProc)))THEN
    ! Particles might impact the EPC on another proc/node. Therefore check if a particle can travel from a local element to an
    ! element that has at least one side, which is an EPC
    ! 3.1) Each processor loops over all of his elements
    iElemLoop: DO iElem = 1+offsetElem, nElems+offsetElem

      iElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,2,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
      iElemRadius = VECNORM ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem),&
                                BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem),&
                                BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

      ! 3.2) Loop over all compute-node elements (every processor loops over all of these elements)
      ! Loop ALL compute-node elements (use global element index)
      iCNElemLoop: DO iCNElem = 1,nComputeNodeTotalElems
        iGlobElem = GetGlobalElemID(iCNElem)

        ! Skip my own elements as they have already been tested when the local sides are checked
        IF(ElementOnProc(iGlobElem)) CYCLE iCNElemLoop

        ! Check if one of the six sides of the compute-node element is a EPC
        ! Note that iSide is in the range of 1:nNonUniqueGlobalSides
        DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
          ! Get BC index of the global side index
          BCIndex = SideInfo_Shared(SIDE_BCID,iSide)
          ! Only check BC sides with BC index > 0
          IF(BCIndex.GT.0)THEN
            ! Get boundary type
            BCType = BoundaryType(BCIndex,BC_TYPE)
            ! Check if EPC has been found
            IF(BCType.EQ.BCTypeEPC)THEN

              ! Check if the BC can be reached
              iGlobElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,2,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,3,iGlobElem)) /) / 2.
              iGlobElemRadius = VECNORM ((/ BoundsOfElem_Shared(2,1,iGlobElem)-BoundsOfElem_Shared(1,1,iGlobElem),&
                                            BoundsOfElem_Shared(2,2,iGlobElem)-BoundsOfElem_Shared(1,2,iGlobElem),&
                                            BoundsOfElem_Shared(2,3,iGlobElem)-BoundsOfElem_Shared(1,3,iGlobElem) /) / 2.)

              ! check if compute-node element "iGlobElem" is within halo_eps of processor-local element "iElem"
              IF (VECNORM( iElemCenter(1:3) - iGlobElemCenter(1:3) ) .LE. ( halo_eps + iElemRadius + iGlobElemRadius ) )THEN
                BCState = BoundaryType(BCIndex,BC_STATE) ! BCState corresponds to iEPC
                IF(BCState.LT.1) CALL abort(__STAMP__,'BCState cannot be <1',IntInfoOpt=BCState)
                iUniqueEPCBC = EPC%Group(BCState,2)
                ! Flag the i-th EPC
                BConProc(iUniqueEPCBC) = .TRUE.
                ! Check if all FPCs have been found -> exit complete loop
                IF(ALL(BConProc)) EXIT iElemLoop
                ! Go to next element
                CYCLE iCNElemLoop
              END IF ! VECNORM( ...
            END IF ! BCType.EQ.BCTypeEPC
          END IF ! BCIndex.GT.0
        END DO ! iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
      END DO iCNElemLoop ! iCNElem = 1,nComputeNodeTotalElems
    END DO iElemLoop ! iElem = 1, nElems
  END IF ! .NOT.(ALL(BConProc))
#endif /*defined(PARTICLES)*/

END IF ! MPIRoot

! 4.) Create MPI sub-communicators
ALLOCATE(EPC%COMM(EPC%nUniqueEPCBounds))
DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  ! Create new communicator
  color = MERGE(iUniqueEPCBC, MPI_UNDEFINED, BConProc(iUniqueEPCBC))

  ! Set communicator id
  EPC%COMM(iUniqueEPCBC)%ID=iUniqueEPCBC

  ! Create new emission communicator for Electric potential boundary condition communication.
  ! Pass MPI_INFO_NULL as rank to follow the original ordering
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, MPI_INFO_NULL, EPC%COMM(iUniqueEPCBC)%UNICATOR, iError)

  ! Find my rank on the shared communicator, comm size and proc name
  IF(BConProc(iUniqueEPCBC))THEN
    CALL MPI_COMM_RANK(EPC%COMM(iUniqueEPCBC)%UNICATOR, EPC%COMM(iUniqueEPCBC)%MyRank, iError)
    CALL MPI_COMM_SIZE(EPC%COMM(iUniqueEPCBC)%UNICATOR, EPC%COMM(iUniqueEPCBC)%nProcs, iError)

    ! inform about size of emission communicator
    IF (EPC%COMM(iUniqueEPCBC)%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
      IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
          WRITE(UNIT_StdOut,'(A,I0,A,I0,A,I0)') ' Electric potential boundary condition: Emission-Communicator ',iUniqueEPCBC,' on ',&
              EPC%COMM(iUniqueEPCBC)%nProcs,' procs for BCState ',EPC%BCState(iUniqueEPCBC)
    END IF
  END IF ! BConProc(iUniqueEPCBC)
END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
DEALLOCATE(BConProc)

! Get the number of procs that actually have a local BC side that is an EPC (required for voltage output to .csv)
! Procs might have zero EPC sides but are in the group because 1.) MPIRoot or 2.) the EPC is in the halo region
! Because only the MPI root process writes the .csv data, the information regarding the voltage on each EPC must be
! communicated with this process even though it might not be connected to each EPC boundary
DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  ASSOCIATE( COMM => EPC%COMM(iUniqueEPCBC)%UNICATOR, nProcsWithSides => EPC%COMM(iUniqueEPCBC)%nProcsWithSides )
    IF(EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
      ! Check if the current processor is actually connected to the EPC via a BC side
      IF(EPC%Group(EPC%BCState(iUniqueEPCBC),3).EQ.0)THEN
        WithSides = 0
      ELSE
        WithSides = 1
      END IF ! EPC%Group(EPC%BCState(iUniqueEPCBC),3).EQ.0
      ! Calculate the sum across the sub-communicator. Only the MPI root process needs this information
      IF(MPIRoot)THEN
        CALL MPI_REDUCE(WithSides, nProcsWithSides, 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, iError)
        ! Sanity check
        IF(nProcsWithSides.EQ.0) CALL abort(__STAMP__,'Found EPC with no processors connected to it')
      ELSE
        CALL MPI_REDUCE(WithSides, 0              , 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, IError)
      END IF ! MPIRoot
    END IF ! EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL
  END ASSOCIATE
END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
#endif /*USE_MPI*/

! When restarting, load the deposited charge on each EPC from the .h5 state file
CALL ReadEPCDataFromH5()

END SUBROUTINE InitEPC


#if defined(PARTICLES)
!===================================================================================================================================
!> Create containers and communicators for each bias-voltage boundary condition where impacting charges are removed and
!> subsequently an electric potential is created (the particle communication is part of the BPO analysis and required here).
!>
!> 1.) Activate bias voltage and check number of boundaries
!> 2.) Get bias voltage parameters
!> 3.) Check if actual bias voltage BC is on current process (or MPI root)
!> 4.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitBV()
! MODULES
USE MOD_Globals                   ,ONLY: CollectiveStop,UNIT_StdOut
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETREAL,GETINT,GETINTARRAY
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
USE MOD_Mesh_Vars                 ,ONLY: nBCs,BoundaryType,BoundaryName
USE MOD_HDG_Vars                  ,ONLY: UseBiasVoltage,BiasVoltage
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcBoundaryParticleOutput,BPO
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_Globals                   ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION,MPI_COMM_PICLAS,MPI_INFO_NULL,MPI_UNDEFINED,MPIRoot
USE MOD_Mesh_Vars                 ,ONLY: nBCSides,BC
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER, PARAMETER :: BCTypeBV(1:3) = (/50,51,52/) ! BCType which allows bias voltage control
!                                                  ! 50: pure DC potential
!                                                  ! 51: cos(wt) function with DC bias
!                                                  ! 52: cos(wt) function with DC bias + coupled power for AC potential adjustment
INTEGER             :: BCType,BVBoundaries,BCState,iBoundary
INTEGER             :: iBC,iPBC
#if USE_MPI
INTEGER             :: color,SideID
LOGICAL             :: BConProc
#endif /*USE_MPI*/
!===================================================================================================================================

!> 1.) Activate bias voltage and check number of boundaries
! Activate model
UseBiasVoltage = GETLOGICAL('UseBiasVoltage')

! Count the number of boundaries that allow bias voltage
BVBoundaries = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(.NOT.ANY(BCType.EQ.BCTypeBV)) CYCLE ! Skip other boundaries
  BVBoundaries = BVBoundaries + 1
END DO

! Skip the following if bias voltage is not active
IF(.NOT.UseBiasVoltage)THEN
  ! Sanity check before returning: bias voltage BCs cannot be used without activating the bias voltage model
  IF(BVBoundaries.GT.0) CALL CollectiveStop(__STAMP__,' Bias voltage BCs require UseBiasVoltage=T!')

  ! Exit this subroutine
  RETURN
END IF

! CalcBoundaryParticleOutput=T and boundaries must be set correctly
IF(.NOT.CalcBoundaryParticleOutput) CALL CollectiveStop(__STAMP__,' UseBiasVoltage=T requires CalcBoundaryParticleOutput=T!')

! Check the number of boundaries that allow bias voltage: Must be exactly 1
IF(BVBoundaries.NE.1) CALL CollectiveStop(__STAMP__,' UseBiasVoltage=T requires exactly one boundary with this feature!')

!> 2.) Get bias voltage parameters
BiasVoltage%NPartBoundaries = GETINT('BiasVoltage-NPartBoundaries')
BiasVoltage%PartBoundaries  = GETINTARRAY('Biasvoltage-PartBoundaries',biasvoltage%npartboundaries)
BiasVoltage%Frequency       = GETReal('BiasVoltage-Frequency')
BiasVoltage%Delta           = GETReal('BiasVoltage-Delta')
#if USE_LOADBALANCE
! Do not nullify during load balance in order to keep the old value on the MPIRoot
IF((.NOT.PerformLoadBalance).OR.(.NOT.MPIRoot))THEN
#endif /*USE_LOADBALANCE*/
  BiasVoltage%BVData = 0.
  ! Update time
  IF(BiasVoltage%Frequency.GT.0.0) BiasVoltage%BVData(3) = 1.0/BiasVoltage%Frequency
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

IF(BiasVoltage%NPartBoundaries.LT.1) CALL CollectiveStop(__STAMP__,' UseBiasVoltage=T requires one or more particle boundaries!')

DO iBoundary=1,BiasVoltage%NPartBoundaries
  iPBC = BiasVoltage%PartBoundaries(iBoundary)
  IF(.NOT.ANY(iPBC.EQ.BPO%PartBoundaries(:))) CALL CollectiveStop(__STAMP__,&
      'One of Biasvoltage-PartBoundaries not defined in any BPO-PartBoundaries')
  iBC = PartBound%MapToFieldBC(iPBC)
  IF(iBC.GT.SIZE(BoundaryName)) CALL CollectiveStop(__STAMP__,'BiasVoltage-PartBoundaries BC index maps to wrong field BCID= ',&
    IntInfo=iBC)
  BCType  = BoundaryType(iBC,BC_TYPE)
  BCState = BoundaryType(iBC,BC_STATE)
  LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)') ' Activated bias voltage by collecting currents from ['//TRIM(BoundaryName(iBC))&
      //'] with BCType [',BCType,'] and BCState [',BCState,']'
END DO

#if USE_MPI
!> 3.) Check if actual bias voltage BC is on current process (or MPI root)
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE
  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(.NOT.ANY(BCType.EQ.BCTypeBV)) CYCLE ! Skip other boundaries
    BConProc = .TRUE.
  END DO ! SideID=1,nBCSides
END IF ! MPIRoot

! 4.) Create MPI sub-communicators
! Create new communicator
color = MERGE(BVBoundaries, MPI_UNDEFINED, BConProc)

! Set communicator id
BiasVoltage%COMM%ID = BVBoundaries

! Create new emission communicator for electric potential boundary condition communication. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, MPI_INFO_NULL, BiasVoltage%COMM%UNICATOR, iError)

! Find my rank on the shared communicator, comm size and process name
IF(BConProc)THEN
  CALL MPI_COMM_RANK(BiasVoltage%COMM%UNICATOR, BiasVoltage%COMM%MyRank, iError)
  CALL MPI_COMM_SIZE(BiasVoltage%COMM%UNICATOR, BiasVoltage%COMM%nProcs, iError)

  ! Inform about size of emission communicator
  IF (BiasVoltage%COMM%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
    IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
        WRITE(UNIT_StdOut,'(A,I0,A,I0)') ' Bias voltage communicator on ',BiasVoltage%COMM%nProcs,' procs for BCState ',BCState
  END IF
END IF ! BConProc
#endif /*USE_MPI*/

! When restarting, load the deposited charge on each EPC from the .h5 state file
CALL ReadBVDataFromH5()

END SUBROUTINE InitBV


!===================================================================================================================================
!> Read the bias voltage (BV) data from a .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes connected to the BV boundary.
!===================================================================================================================================
SUBROUTINE ReadBVDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals          ,ONLY: UNIT_stdOut,MPIRoot,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5          ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars     ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input       ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars         ,ONLY: BiasVoltage,BVDataLength
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(255) :: ContainerName
LOGICAL        :: BVExists
REAL           :: BVDataHDF5(1:BVDataLength)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='BiasVoltage'
  CALL DatasetExists(File_ID,TRIM(ContainerName),BVExists)
  ! Check for new parameter name
  IF(BVExists)THEN
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(BVDataLength,IK)/) , 0_IK , 1 , RealArray=BVDataHDF5)
    WRITE(UNIT_stdOut,'(3(A,ES10.2E3))') " Read bias voltage from restart file ["//TRIM(RestartFile)//&
        "] Bias voltage[V]: ",BVDataHDF5(1),", Ion excess[C]: ",BVDataHDF5(2),", next adjustment time[s]: ",BVDataHDF5(3)
    BiasVoltage%BVData = BVDataHDF5
  END IF ! BVExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
CALL SynchronizeBV()
#endif /*USE_MPI*/

END SUBROUTINE ReadBVDataFromH5


#if USE_MPI
!===================================================================================================================================
!> Communicate the bias voltage values from MPIRoot to sub-communicator processes
!===================================================================================================================================
SUBROUTINE SynchronizeBV()
! MODULES
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
USE MOD_HDG_Vars ,ONLY: BiasVoltage,BVDataLength
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(BiasVoltage%COMM%UNICATOR.NE.MPI_COMM_NULL)THEN
  ! Broadcast from root to other processors on the sub-communicator
  CALL MPI_BCAST(BiasVoltage%BVData, BVDataLength, MPI_DOUBLE_PRECISION, 0, BiasVoltage%COMM%UNICATOR, IERROR)
END IF
END SUBROUTINE SynchronizeBV
#endif /*USE_MPI*/
#endif /*defined(PARTICLES)*/


!===================================================================================================================================
!> Read the electric charge that resides on each FPC boundary from .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes for each FPC
!===================================================================================================================================
SUBROUTINE ReadFPCDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals            ,ONLY: UNIT_stdOut,MPIRoot,IERROR,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5            ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input         ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars           ,ONLY: FPC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueFPCBC,nFPCHDF5,nDimsFPC
CHARACTER(255)     :: ContainerName
CHARACTER(1000)    :: TmpStr
LOGICAL            :: FPCExists
REAL,ALLOCATABLE   :: FPCDataHDF5(:)
INTEGER(HSIZE_T), POINTER :: HSizeFPC(:)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='FloatingPotentialCharge'
  CALL DatasetExists(File_ID,TRIM(ContainerName),FPCExists)
  ! Check for new parameter name
  IF(FPCExists)THEN
    CALL GetDataSize(File_ID,TRIM(ContainerName),nDimsFPC,HSizeFPC)
    CHECKSAFEINT(HSizeFPC(2),4)
    nFPCHDF5=INT(HSizeFPC(2),4)
    DEALLOCATE(HSizeFPC)
    IF(nFPCHDF5.NE.FPC%nUniqueFPCBounds)THEN
      WRITE(UNIT_StdOut,*) "nFPCHDF5 (restart file) must be equal FPC%nUniqueFPCBounds"
      WRITE(UNIT_StdOut,*) "nFPCHDF5             =", nFPCHDF5
      WRITE(UNIT_StdOut,*) "FPC%nUniqueFPCBounds =", FPC%nUniqueFPCBounds
      CALL abort(__STAMP__,'Restarting with a different number of FPC boundaries, which is not implemented!')
    END IF ! nFPCHDF5.NE.FPC%nUniqueFPCBounds

    ! Allocate the containers
    ALLOCATE(FPCDataHDF5(FPC%nUniqueFPCBounds))
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(FPC%nUniqueFPCBounds,IK)/) , 0_IK , 1 , RealArray=FPCDataHDF5)
    TmpStr=''
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      ! Output in this format: ", 1: 1.312e2 [C]" + ", 2: 3.352e3 [C]" + ...
      WRITE(TmpStr,'(A,I0,A,ES10.3,A)') TRIM(TmpStr)//', ',iUniqueFPCBC,': ',FPCDataHDF5(iUniqueFPCBC),' [C]'
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
    TmpStr(1:1) = ' ' ! Remove first comma
    TmpStr      = ADJUSTL(TmpStr) ! Remove leading white spaces
    WRITE(UNIT_stdOut,'(A)') " Read floating boundary condition charges from restart file ["//TRIM(RestartFile)//"]: "//TRIM(TmpStr)
    FPC%Charge(:) = FPCDataHDF5
    DEALLOCATE(FPCDataHDF5)
  END IF ! FPCExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each FPC
CALL SynchronizeChargeOnFPC()
#endif /*USE_MPI*/

END SUBROUTINE ReadFPCDataFromH5


!===================================================================================================================================
!> Read the electric potential (voltage) that was applied on each EPC boundary from .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
!===================================================================================================================================
SUBROUTINE ReadEPCDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals            ,ONLY: UNIT_stdOut,MPIRoot,IERROR,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5            ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input         ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars           ,ONLY: EPC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueEPCBC,nEPCHDF5,nDimsEPC
CHARACTER(255)     :: ContainerName
CHARACTER(1000)    :: TmpStr
LOGICAL            :: EPCExists
REAL,ALLOCATABLE   :: EPCDataHDF5(:)
INTEGER(HSIZE_T), POINTER :: HSizeEPC(:)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='ElectricPotenitalCondition'
  CALL DatasetExists(File_ID,TRIM(ContainerName),EPCExists)
  ! Check for new parameter name
  IF(EPCExists)THEN
    CALL GetDataSize(File_ID,TRIM(ContainerName),nDimsEPC,HSizeEPC)
    CHECKSAFEINT(HSizeEPC(2),4)
    nEPCHDF5=INT(HSizeEPC(2),4)
    DEALLOCATE(HSizeEPC)
    IF(nEPCHDF5.NE.EPC%nUniqueEPCBounds)THEN
      WRITE(UNIT_StdOut,*) "nEPCHDF5 (restart file) must be equal EPC%nUniqueEPCBounds"
      WRITE(UNIT_StdOut,*) "nEPCHDF5             =", nEPCHDF5
      WRITE(UNIT_StdOut,*) "EPC%nUniqueEPCBounds =", EPC%nUniqueEPCBounds
      CALL abort(__STAMP__,'Restarting with a different number of EPC boundaries, which is not implemented!')
    END IF ! nEPCHDF5.NE.EPC%nUniqueEPCBounds

    ! Allocate the containers
    ALLOCATE(EPCDataHDF5(EPC%nUniqueEPCBounds))
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(EPC%nUniqueEPCBounds,IK)/) , 0_IK , 1 , RealArray=EPCDataHDF5)
    TmpStr=''
    DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
      ! Output in this format: ", 1: 1.312e2 [V]" + ", 2: 3.352e3 [V]" + ...
      WRITE(TmpStr,'(A,I0,A,ES10.3,A)') TRIM(TmpStr)//', ',iUniqueEPCBC,': ',EPCDataHDF5(iUniqueEPCBC),' [V]'
    END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
    TmpStr(1:1) = ' ' ! Remove first comma
    TmpStr      = ADJUSTL(TmpStr) ! Remove leading white spaces
    WRITE(UNIT_stdOut,'(A)') " Read electric potential condition from restart file ["//TRIM(RestartFile)//"]: "//TRIM(TmpStr)
    EPC%Voltage(:) = EPCDataHDF5
    DEALLOCATE(EPCDataHDF5)
  END IF ! EPCExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
CALL SynchronizeVoltageOnEPC()
#endif /*USE_MPI*/

END SUBROUTINE ReadEPCDataFromH5


#if USE_MPI
!===================================================================================================================================
!> The MPI root process distributes the information among the sub-communicator processes for each FPC
!===================================================================================================================================
SUBROUTINE SynchronizeChargeOnFPC()
! MODULES
USE MOD_HDG_Vars ,ONLY: FPC
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueFPCBC
!===================================================================================================================================
DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR )
    IF(COMM.NE.MPI_COMM_NULL)THEN
      ! Broadcast from root to other processors on the sub-communicator
      CALL MPI_BCAST(FPC%Charge(iUniqueFPCBC), 1, MPI_DOUBLE_PRECISION, 0, COMM, IERROR)
    END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
  END ASSOCIATE
END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
END SUBROUTINE SynchronizeChargeOnFPC


!===================================================================================================================================
!> The MPI root process distributes the information among the sub-communicator processes for each EPC
!===================================================================================================================================
SUBROUTINE SynchronizeVoltageOnEPC()
! MODULES
USE MOD_HDG_Vars ,ONLY: EPC
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueEPCBC
!===================================================================================================================================
DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  ASSOCIATE( COMM => EPC%COMM(iUniqueEPCBC)%UNICATOR )
    IF(COMM.NE.MPI_COMM_NULL)THEN
      ! Broadcast from root to other processors on the sub-communicator
      CALL MPI_BCAST(EPC%Voltage(iUniqueEPCBC), 1, MPI_DOUBLE_PRECISION, 0, COMM, IERROR)
    END IF ! EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL
  END ASSOCIATE
END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
END SUBROUTINE SynchronizeVoltageOnEPC
#endif /*USE_MPI*/


!===================================================================================================================================
!> HDG solver for linear or non-linear systems
!===================================================================================================================================
#if defined(PARTICLES)
SUBROUTINE HDG(t,U_out,iter,ForceCGSolverIteration_opt)
#else
SUBROUTINE HDG(t,U_out,iter)
#endif /*defined(PARTICLES)*/
! MODULES
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
USE MOD_TimeDisc_Vars   ,ONLY: dt,dt_Min
USE MOD_Equation_Vars   ,ONLY: E,Et
USE MOD_Globals_Vars    ,ONLY: eps0
USE MOD_Analyze_Vars    ,ONLY: CalcElectricTimeDerivative
USE MOD_Analyze_Vars    ,ONLY: FieldAnalyzeStep
USE MOD_Dielectric_vars ,ONLY: DoDielectric,isDielectricElem,ElemToDielectric,DielectricEps
#endif /*(USE_HDG && (PP_nVar==1))*/
#if defined(PARTICLES)
USE MOD_HDG_Vars        ,ONLY: UseEPC,EPC
#endif /*defined(PARTICLES)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: t !time
INTEGER(KIND=8),INTENT(IN)  :: iter
#if defined(PARTICLES)
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only required for BR electron fluid)
#endif /*defined(PARTICLES)*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_out(PP_nVar,nGP_vol,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(PARTICLES)
LOGICAL :: ForceCGSolverIteration_loc
INTEGER :: iUniqueEPCBC
#endif /*defined(PARTICLES)*/
#if (USE_HDG && (PP_nVar==1))
INTEGER           :: iDir,iElem
#endif /*(USE_HDG && (PP_nVar==1))*/
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(4))
#endif /*EXTRAE*/

! Calculate temporal derivate of D in last iteration before Analyze_dt is reached: Store E^n here
#if (USE_HDG && (PP_nVar==1))
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  IF (iStage.EQ.1) THEN
#endif
  IF(CalcElectricTimeDerivative)THEN
    ! iter is incremented after this function and then checked in analyze routine with iter+1
    IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.ALMOSTEQUAL(dt,dt_Min(DT_END)).OR.(MOD(iter+1,FieldAnalyzeStep).EQ.0))THEN
      ! Store old E-field
      Et(:,:,:,:,:) = E(:,:,:,:,:)
    END IF
  END IF ! CalcElectricTimeDerivative
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
END IF
#endif
#endif /*(USE_HDG && (PP_nVar==1))*/

! Check whether the solver should be skipped in this iteration
IF (iter.GT.0 .AND. HDGSkip.NE.0) THEN
  IF (t.LT.HDGSkip_t0) THEN
    IF (MOD(iter,INT(HDGSkipInit,8)).NE.0) RETURN
  ELSE
    IF (MOD(iter,INT(HDGSkip,8)).NE.0) RETURN
  END IF
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  IF (iStage.GT.1) THEN
    RETURN
  END IF
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
      ASSOCIATE( COMM => EPC%COMM(iUniqueEPCBC)%UNICATOR)
        IF(EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
          IF(MPIRoot)THEN
            CALL MPI_REDUCE(MPI_IN_PLACE, EPC%ChargeProc(iUniqueEPCBC), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IERROR)
          ELSE
            CALL MPI_REDUCE(EPC%ChargeProc(iUniqueEPCBC), 0           , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IERROR)
          END IF ! MPIRoot
          EPC%Charge(iUniqueEPCBC) = EPC%Charge(iUniqueEPCBC) + EPC%ChargeProc(iUniqueEPCBC)
        END IF ! EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL
      END ASSOCIATE
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
#endif /*(USE_HDG && (PP_nVar==1)&&  defined(PARTICLES))*/

! Run the appropriate HDG solver: Newton or Linear
#if defined(PARTICLES)
IF(UseBRElectronFluid) THEN
  IF (HDGNonLinSolver.EQ.1) THEN

    IF (PRESENT(ForceCGSolverIteration_opt)) THEN; ForceCGSolverIteration_loc = ForceCGSolverIteration_opt
    ELSE;                                          ForceCGSolverIteration_loc = .FALSE.
    END IF

    CALL HDGNewton(t, U_out, iter, ForceCGSolverIteration_loc)
  ELSE
    CALL CollectiveStop(__STAMP__,'Defined HDGNonLinSolver not implemented (HDGFixPoint has been removed!) HDGNonLinSolver = ',&
    IntInfo=HDGNonLinSolver)
  END IF
ELSE
#endif /*defined(PARTICLES)*/
  CALL HDGLinear(t,U_out)
#if defined(PARTICLES)
END IF
#endif /*defined(PARTICLES)*/

! Calculate temporal derivate of D in last iteration before Analyze_dt is reached: Use E^n+1 here and calculate the derivative dD/dt
#if (USE_HDG && (PP_nVar==1))
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
IF (iStage.EQ.nRKStages) THEN
#endif
  IF(CalcElectricTimeDerivative)THEN
    ! iter is incremented after this function and then checked in analyze routine with iter+1
    IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.ALMOSTEQUAL(dt,dt_Min(DT_END)).OR.(MOD(iter+1,FieldAnalyzeStep).EQ.0))THEN
      IF(DoDielectric)THEN
        DO iElem=1,PP_nElems
          IF(isDielectricElem(iElem)) THEN
            DO iDir = 1, 3
              Et(iDir,:,:,:,iElem) = DielectricEps(:,:,:,ElemToDielectric(iElem))*eps0*(E(iDir,:,:,:,iElem)-Et(iDir,:,:,:,iElem)) / dt
            END DO ! iDir = 1, 3
          ELSE
            Et(:,:,:,:,iElem) = eps0*(E(:,:,:,:,iElem)-Et(:,:,:,:,iElem)) / dt
          END IF ! isDielectricElem(iElem)
        END DO ! iElem=1,PP_nElems
      ELSE
        Et(:,:,:,:,:) = eps0*(E(:,:,:,:,:)-Et(:,:,:,:,:)) / dt
      END IF ! DoDielectric
    END IF
  END IF ! CalcElectricTimeDerivative
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
END IF
#endif
#endif /*(USE_HDG && (PP_nVar==1))*/

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
END SUBROUTINE HDG


!===================================================================================================================================
!> Linear HDG solver
!===================================================================================================================================
SUBROUTINE HDGLinear(time,U_out)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_Equation           ,ONLY: CalcSourceHDG,ExactFunc
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Equation_Vars      ,ONLY: chitens_face
USE MOD_Mesh_Vars          ,ONLY: Face_xGP,BoundaryType,nSides,BC
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,NormVec,SurfElem
USE MOD_Interpolation_Vars ,ONLY: wGP
USE MOD_Elem_Mat           ,ONLY: PostProcessGradient
USE MOD_FillMortar_HDG     ,ONLY: SmallToBigMortar_HDG
#if (PP_nVar==1)
USE MOD_Equation_Vars      ,ONLY: E
#elif (PP_nVar==3)
USE MOD_Equation_Vars      ,ONLY: B
#else
USE MOD_Equation_Vars      ,ONLY: B, E
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
USE PETSc
USE MOD_Mesh_Vars          ,ONLY: SideToElem
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars
#endif
USE MOD_FillMortar_HDG     ,ONLY: BigToSmallMortar_HDG
#endif
USE MOD_Globals_Vars    ,ONLY: ElementaryCharge,eps0
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: time !time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_out(PP_nVar,nGP_vol,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,r,p,q,iElem, iVar!,iter
INTEGER :: BCsideID,BCType,BCState,SideID,iLocSide
REAL    :: RHS_face(PP_nVar,nGP_face,nSides)
REAL    :: rtmp(nGP_vol)
!LOGICAL :: converged
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol,PP_nElems)
#endif
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
PetscErrorCode       :: ierr
PetscScalar, POINTER :: lambda_pointer(:)
KSPConvergedReason   :: reason
PetscInt             :: iterations
PetscReal            :: petscnorm
INTEGER              :: ElemID,iBCSide,PETScLocalID
INTEGER              :: PETScID_start, PETScID_stop
REAL                 :: timeStartPiclas,timeEndPiclas
REAL                 :: RHS_conductor(nGP_face)
#endif
#if USE_PETSC
INTEGER              :: iUniqueFPCBC
#endif /*USE_PETSC*/
!===================================================================================================================================
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iVar = 1, PP_nVar
  !Dirichlet boundary conditions
#if (PP_nVar!=1)
  IF (iVar.EQ.4) THEN
#endif
  DO BCsideID=1,nDirichletBCSides
    SideID=DirichletBC(BCsideID)
    BCType =BoundaryType(BC(SideID),BC_TYPE)
    BCState=BoundaryType(BC(SideID),BC_STATE)
    SELECT CASE(BCType)
    CASE(2) ! exact BC = Dirichlet BC !! ExactFunc via BCState (time is optional)
      ! Determine the exact BC state
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(BCState,Face_xGP(:,p,q,SideID),lambda(iVar,r:r,SideID),t=time)
      END DO; END DO !p,q
    CASE(4) ! exact BC = Dirichlet BC !! Zero potential
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
       lambda(iVar,r:r,SideID)=0.
      END DO; END DO !p,q
    CASE(5,51,52,60) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (with zero crossing)
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-1,Face_xGP(:,p,q,SideID),lambda(iVar,r:r,SideID),t=time,iRefState=BCState)
      END DO; END DO !p,q
    CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (without zero crossing)
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-2,Face_xGP(:,p,q,SideID),lambda(iVar,r:r,SideID),t=time,iRefState=BCState)
      END DO; END DO !p,q
    CASE(7) ! exact BC = Dirichlet BC !! ExactFunc via LinState (time is optional)for linear potential function
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-3,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,iLinState=BCState)
      END DO; END DO !p,q
    CASE(8) ! exact BC = Dirichlet BC !! ExactFunc via electric potential and decharing of a surface
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-4,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,BCState=BCState)
      END DO; END DO !p,q
    CASE(50) ! exact BC = Dirichlet BC !! ExactFunc via DC bias voltage
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-5,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,BCState=BCState)
      END DO; END DO !p,q
    END SELECT ! BCType
  END DO !BCsideID=1,nDirichletBCSides
#if (PP_nVar!=1)
  END IF
#endif

  !neumann BC
  DO BCsideID=1,nNeumannBCSides
    SideID=NeumannBC(BCsideID)
    BCType =BoundaryType(BC(SideID),BC_TYPE)
    BCState=BoundaryType(BC(SideID),BC_STATE)
    SELECT CASE(BCType)
    CASE(10) !neumann q=0 !! Zero gradient
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        qn_face(iVar,r,BCSideID)= 0.
      END DO; END DO !p,q
    CASE(11) !neumann q*n=1 !test
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        qn_face(iVar,r,BCSideID)=SUM((/1.,1.,1./)  &
                            *MATMUL(chitens_face(:,:,p,q,SideID),NormVec(:,p,q,SideID)))*SurfElem(p,q,SideID)*wGP(p)*wGP(q)
      END DO; END DO !p,q
    END SELECT ! BCType
  END DO !BCsideID=1,nNeumannBCSides

!for magnetostatic only neumann
#if (PP_nVar!=1)
  IF (iVar.LT.4) THEN
    DO BCsideID=1,nDirichletBCSides
!      SideID=DirichletBC(BCsideID)
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        qn_face_MagStat(iVar,r,BCSideID)= 0.
      END DO; END DO !p,q
    END DO !BCsideID=1,nDirichletBCSides
  END IF
#endif

  ! Floating boundary BCs
#if USE_PETSC
  IF(UseFPC)THEN
#if USE_MPI
    ! Communicate the accumulated charged on each BC to all processors on the communicator
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR)
        IF(FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, FPC%ChargeProc(iUniqueFPCBC), 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, IERROR)
          FPC%Charge(iUniqueFPCBC) = FPC%Charge(iUniqueFPCBC) + FPC%ChargeProc(iUniqueFPCBC)
        END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
      END ASSOCIATE
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
#else
    FPC%Charge = FPC%Charge + FPC%ChargeProc
#endif /*USE_MPI*/
    FPC%ChargeProc = 0.
    ! Apply charge to RHS, which this is done below: RHS_conductor(1)=FPC%Charge(iUniqueFPCBC)/eps0

  END IF ! UseFPC
#endif /*USE_PETSC*/

  ! Set potential to zero
  IF(SetZeroPotentialDOF) lambda(iVar,1,1) = 0.
END DO

!volume source (volume RHS of u system)
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    CALL CalcSourceHDG(i,j,k,iElem,RHS_vol(1:PP_nVar,r,iElem))
  END DO; END DO; END DO !i,j,k
  DO iVar = 1, PP_nVar
    RHS_Vol(iVar,:,iElem)=-JwGP_vol(:,iElem)*RHS_vol(iVar,:,iElem)
  END DO
END DO !iElem

!replace lambda with exact function (debugging)
IF(ExactLambda)THEN
  DO SideID=1,nSides
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(IniExactFunc,Face_xGP(:,p,q,SideID),lambda( 1:PP_nVar,r,SideID))
    END DO; END DO !p,q
  END DO
END IF

!prepare RHS_face ( RHS for lamdba system.)
DO iVar = 1, PP_nVar
  RHS_face(iVar,:,:)=0.
  DO iElem=1,PP_nElems
    !rtmp=MATMUL(InvDhat(:,:,iElem),-RHS_loc(:,iElem))
    CALL DSYMV('U',nGP_vol,1., InvDhat(:,:,iElem),nGP_vol, &
                               -RHS_vol(iVar,:,iElem),1,0., &
                               rtmp(:),1)

    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      CALL DGEMV('N',nGP_face,nGP_vol,1., &
                          Ehat(:,:,iLocSide,iElem), nGP_face, &
                          rtmp,1,1.,& !add to RHS_face
                          RHS_face(iVar,:,SideID),1)
    END DO
  END DO !iElem
END DO !ivar

!add Neumann
DO BCsideID=1,nNeumannBCSides
  SideID=NeumannBC(BCsideID)
  RHS_face(:,:,SideID)=RHS_face(:,:,SideID)+qn_face(:,:,BCSideID)
END DO

#if USE_PETSC
! add Dirichlet contribution
DO iBCSide=1,nDirichletBCSides
  BCSideID=DirichletBC(iBCSide)
  ElemID    = SideToElem(S2E_ELEM_ID,BCSideID)
  DO iLocSide=1,6
    SideID = ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
    IF(PETScGlobal(SideID).EQ.-1) CYCLE
    CALL DGEMV('N',nGP_face,nGP_face,-1., &
                          Smat_BC(:,:,iLocSide,iBCSide), nGP_face, &
                          lambda(1,:,BCSideID),1,1.,& !add to RHS_face
                          RHS_face(1,:,SideID),1)
  END DO
END DO
#endif

#if (PP_nVar!=1)
DO iVar = 1, PP_nVar
  IF (iVar.LT.4) THEN
    DO BCsideID=1,nDirichletBCSides
      SideID=DirichletBC(BCsideID)
      RHS_face(iVar,:,SideID)=RHS_face(iVar,:,SideID)+qn_face_MagStat(iVar,:,BCSideID)
    END DO !BCsideID=1,nDirichletBCSides
  END IF
END DO
#endif

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
CALL Mask_MPIsides(PP_nVar,RHS_face)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
CALL SmallToBigMortar_HDG(PP_nVar,RHS_face(1:PP_nVar,1:nGP_Face,1:nSides))
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

! SOLVE
DO iVar=1, PP_nVar

#if USE_PETSC
  ! Fill right hand side
  PetscCallA(VecZeroEntries(RHS_petsc,ierr))
  TimeStartPiclas=PICLASTIME()
  DO PETScLocalID=1,nPETScUniqueSides
    SideID=PETScLocalToSideID(PETScLocalID)
    !VecSetValuesBlockedLocal somehow not working...
    PetscCallA(VecSetValuesBlocked(RHS_petsc,1,PETScGlobal(SideID),RHS_face(1,:,SideID),INSERT_VALUES,ierr))
  END DO
  ! The MPIRoot process has charge and voltage of all FPCs, there, this process sets all conductor RHS information
  IF(MPIRoot)THEN
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      RHS_conductor(:)=0.
      RHS_conductor(1)=FPC%Charge(iUniqueFPCBC)/eps0
      PetscCallA(VecSetValuesBlocked(RHS_petsc,1,nPETScUniqueSidesGlobal-1-FPC%nUniqueFPCBounds+iUniqueFPCBC,RHS_conductor,INSERT_VALUES,ierr))
    END DO !iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  END IF ! MPIRoot

  ! Reset the RHS of the first DOF if ZeroPotential must be set
  IF(MPIroot .AND. SetZeroPotentialDOF) THEN
    PetscCallA(VecSetValue(RHS_petsc,0,0,INSERT_VALUES,ierr))
  END IF

  PetscCallA(VecAssemblyBegin(RHS_petsc,ierr))
  PetscCallA(VecAssemblyEnd(RHS_petsc,ierr))

  ! Calculate lambda
  PetscCallA(KSPSolve(ksp,RHS_petsc,lambda_petsc,ierr))
  TimeEndPiclas=PICLASTIME()
  PetscCallA(KSPGetIterationNumber(ksp,iterations,ierr))
  PetscCallA(KSPGetConvergedReason(ksp,reason,ierr))
  PetscCallA(KSPGetResidualNorm(ksp,petscnorm,ierr))
  ! reason - negative value indicates diverged, positive value converged, see KSPConvergedReason
  !  -2: KSP_DIVERGED_NULL
  !  -3: KSP_DIVERGED_ITS            -> Ran out of iterations before any convergence criteria was reached
  !  -4: KSP_DIVERGED_DTOL           -> norm(r) >= dtol*norm(b)
  !  -5: KSP_DIVERGED_BREAKDOWN      -> Breakdown in the Krylov method: the method could not continue to enlarge the Krylov space
  !  -6: KSP_DIVERGED_BREAKDOWN_BICG -> Breakdown in the Krylov method: the method could not continue to enlarge the Krylov space
  !  -7: KSP_DIVERGED_NONSYMMETRIC   -> It appears the operator or preconditioner is not symmetric and this Krylov method
  !  -8: KSP_DIVERGED_INDEFINITE_PC  -> It appears the preconditioner is indefinite
  !  -9: KSP_DIVERGED_NANORINF
  ! -10: KSP_DIVERGED_INDEFINITE_MAT
  ! -11: KSP_DIVERGED_PC_FAILED      -> It was not possible to build or use the requested preconditioner
  ! -11: KSP_DIVERGED_PCSETUP_FAILED_DEPRECATED
  IF(reason.LT.0)THEN
    SWRITE(*,*) 'Attention: PETSc not converged! Reason: ', reason
  END IF
  IF(MPIroot) CALL DisplayConvergence(TimeEndPiclas-TimeStartPiclas, iterations, petscnorm)

  ! Fill element local lambda for post processing
  PetscCallA(VecScatterBegin(scatter_petsc, lambda_petsc, lambda_local_petsc, INSERT_VALUES, SCATTER_FORWARD,ierr))
  PetscCallA(VecScatterEnd(scatter_petsc, lambda_petsc, lambda_local_petsc, INSERT_VALUES, SCATTER_FORWARD,ierr))
  PetscCallA(VecGetArrayReadF90(lambda_local_petsc,lambda_pointer,ierr))
  DO PETScLocalID=1,nPETScUniqueSides
    SideID=PETScLocalToSideID(PETScLocalID)
    PETScID_start=1+(PETScLocalID-1)*nGP_face
    PETScID_stop=PETScLocalID*nGP_face
    lambda(1,:,SideID) = lambda_pointer(PETScID_start:PETScID_stop)
  END DO
  PetscCallA(VecRestoreArrayReadF90(lambda_local_petsc,lambda_pointer,ierr))

  ! Fill Conductor lambda
  IF(UseFPC)THEN
    PetscCallA(VecScatterBegin(scatter_conductors_petsc, lambda_petsc, lambda_local_conductors_petsc, INSERT_VALUES, SCATTER_FORWARD,ierr))
    PetscCallA(VecScatterEnd(scatter_conductors_petsc, lambda_petsc, lambda_local_conductors_petsc, INSERT_VALUES, SCATTER_FORWARD,ierr))
    PetscCallA(VecGetArrayReadF90(lambda_local_conductors_petsc,lambda_pointer,ierr))
    FPC%VoltageProc = 0. ! nullify just to be safe
    ! TODO multiple conductors
    DO BCsideID=1,nConductorBCsides
      SideID=ConductorBC(BCSideID)
      BCState=BoundaryType(BC(SideID),BC_STATE)
      DO i=1,nGP_face
        lambda(1,:,SideID) = lambda_pointer(1 + (FPC%Group(BCState,2) - 1) * nGP_face)
      END DO
      ! Copy the value to the FPC container for output to .csv (only mpi root does this)
      FPC%VoltageProc(FPC%Group(BCState,2)) = lambda(1,1,SideID)
    END DO
#if USE_MPI
    ! Sum the voltages across each sub-group and divide by the size of the group to get the voltage. This is required if the MPI
    ! root is not connected to every FPC (which is usually the case)
    !
    ! 1. Communicate the accumulated voltage (should be the same on each proc that has such a BC) on each BC to all processors on the
    ! communicator
    FPC%Voltage = 0.
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR )
        IF(FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
          ASSOCIATE( VoltageProc => FPC%VoltageProc(iUniqueFPCBC), Voltage => FPC%Voltage(iUniqueFPCBC) )
            IF(MPIroot)THEN
              CALL MPI_REDUCE(VoltageProc, Voltage, 1 ,MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, iError)
            ELSE
              CALL MPI_REDUCE(VoltageProc, 0      , 1 ,MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IError)
            END IF ! MPIroot
          END ASSOCIATE
        END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
      END ASSOCIATE
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds

    ! 2. Divide by group size (procs that have an actual FPC BC side -> FPC%COMM(iUniqueFPCBC)%nProcsWithSides).
    !    The MPI root process definitely has the size of each group, which is at least 1.
    IF(MPIRoot)THEN
      DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
        FPC%Voltage(iUniqueFPCBC) = FPC%Voltage(iUniqueFPCBC) / REAL(FPC%COMM(iUniqueFPCBC)%nProcsWithSides)
      END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
    END IF ! MPIRoot
#else
    FPC%Voltage = FPC%VoltageProc
#endif /*USE_MPI*/
    PetscCallA(VecRestoreArrayReadF90(lambda_local_conductors_petsc,lambda_pointer,ierr))
  END IF ! UseFPC

  ! PETSc Calculate lambda at small mortars from big mortars
  CALL BigToSmallMortar_HDG(1,lambda)
#if USE_MPI
  CALL StartReceiveMPIData(1,lambda,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
  CALL StartSendMPIData(   1,lambda,1,nSides,SendRequest_U,SendID=1) ! Send MINE
  CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif
#else
  CALL CG_solver(RHS_face(iVar,:,:),lambda(iVar,:,:),iVar)
#endif
  !POST PROCESSING

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  !post processing:
  DO iElem=1,PP_nElems
    ! for post-proc
    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      CALL DGEMV('T',nGP_face,nGP_vol,1., &
                          Ehat(:,:,iLocSide,iElem), nGP_face, &
                          lambda(iVar,:,SideID),1,1.,& !add to RHS_face
                          RHS_vol(iVar,:,iElem),1)
    END DO
    !U_out(:,iElem)=MATMUL(InvDhat(:,:,iElem),-RHS_loc(:,iElem))
    CALL DSYMV('U',nGP_vol,1., InvDhat(:,:,iElem),nGP_vol, &
                               -RHS_vol(iVar,:,iElem),1,0., &
                               U_out(iVar,:,iElem),1)
  END DO !iElem
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
END DO !iVar

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
#if (PP_nVar==1)
  CALL PostProcessGradient(U_out(1,:,:),lambda(1,:,:),E)
#elif (PP_nVar==3)
  DO iVar=1, PP_nVar
    CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  END DO
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
    B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
    B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  END DO; END DO; END DO !i,j,k
#else
  DO iVar=1, 3
    CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  END DO
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
    B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
    B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  END DO; END DO; END DO !i,j,k
  CALL PostProcessGradient(U_out(4,:,:),lambda(4,:,:),E)
#endif
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE HDGLinear


!===================================================================================================================================
!> HDG non-linear solver via Newton's method
!===================================================================================================================================
#if defined(PARTICLES)
SUBROUTINE HDGNewton(time,U_out,td_iter,ForceCGSolverIteration_opt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_Equation               ,ONLY: CalcSourceHDG,ExactFunc
USE MOD_FillMortar_HDG         ,ONLY: SmallToBigMortar_HDG
#if defined(IMPA) || defined(ROS)
USE MOD_LinearSolver_Vars      ,ONLY: DoPrintConvInfo
#else
USE MOD_TimeDisc_Vars          ,ONLY: IterDisplayStep,DoDisplayIter
#endif
USE MOD_Globals_Vars           ,ONLY: eps0
USE MOD_Equation_Vars          ,ONLY: chitens_face
USE MOD_Mesh_Vars              ,ONLY: Face_xGP,BoundaryType,nSides,BC
USE MOD_Mesh_Vars              ,ONLY: ElemToSide,NormVec,SurfElem
USE MOD_Interpolation_Vars     ,ONLY: wGP
USE MOD_Elem_Mat               ,ONLY: PostProcessGradient, Elem_Mat,BuildPrecond
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartTime
#if (PP_nVar==1)
USE MOD_Equation_Vars          ,ONLY: E
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_HDG_Vars               ,ONLY: ElemToBRRegion,RegionElectronRef
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: time !time
INTEGER(KIND=8),INTENT(IN)  :: td_iter
#if defined(PARTICLES)
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only BR electron fluid)
#endif /*defined(PARTICLES)*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_out(PP_nVar,nGP_vol,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,r,p,q,iElem, iter,RegionID
INTEGER :: BCsideID,BCType,BCState,SideID,iLocSide
REAL    :: RHS_face(PP_nVar,nGP_face,nSides)
REAL    :: rtmp(nGP_vol),Norm_r2!,Norm_r2_old
LOGICAL :: converged, beLinear
LOGICAL :: warning_linear
REAL    :: warning_linear_phi
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol,PP_nElems)
#endif
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
#if (PP_nVar!=1)
  CALL CollectiveStop(__STAMP__,'Nonlinear Newton solver only available with EQ-system Poisson!')
#endif
Norm_r2=0.

!Dirichlet boundary conditions
DO BCsideID=1,nDirichletBCSides
  SideID=DirichletBC(BCsideID)
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(2) ! exact BC = Dirichlet BC !! ExactFunc via BCState (time is optional)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(BCState,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),time)
    END DO; END DO !p,q
  CASE(4) ! exact BC = Dirichlet BC !! Zero potential
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      lambda(PP_nVar,r:r,SideID)= 0.
    END DO; END DO !p,q
  CASE(5) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (with zero crossing)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-1,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,iRefState=BCState)
    END DO; END DO !p,q
  CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (without zero crossing)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-2,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,iRefState=BCState)
    END DO; END DO !p,q
  CASE(7) ! exact BC = Dirichlet BC !! ExactFunc via LinState (time is optional)for linear potential function
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-3,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,iLinState=BCState)
    END DO; END DO !p,q
    CASE(8,50,51,52,60) ! exact BC = Dirichlet BC !! ExactFunc via electric potential and decharing of a surface
      CALL abort(__STAMP__,'Dirichlet BC=8,50,51,52,60 model not implemented for HDG Newton!')
  END SELECT ! BCType
END DO !BCsideID=1,nDirichletBCSides

!neumann BC
DO BCsideID=1,nNeumannBCSides
  SideID=NeumannBC(BCsideID)
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(10) !neumann q=0 !! Zero gradient
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      qn_face(PP_nVar,r,BCSideID)= 0.
    END DO; END DO !p,q
  CASE(11) !neumann q*n=1 !test
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      qn_face(PP_nVar,r,BCSideID)=SUM((/1.,1.,1./)  &
                          *MATMUL(chitens_face(:,:,p,q,SideID),NormVec(:,p,q,SideID)))*SurfElem(p,q,SideID)*wGP(p)*wGP(q)
    END DO; END DO !p,q
  END SELECT ! BCType
END DO !BCsideID=1,nNeumannBCSides

warning_linear=.FALSE.
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    CALL CalcSourceHDG(i,j,k,iElem,RHS_vol(1:PP_nVar,r,iElem),U_out(1,r,iElem),warning_linear,warning_linear_phi)
  END DO; END DO; END DO !i,j,k
  RHS_Vol(PP_nVar,:,iElem)=-JwGP_vol(:,iElem)*RHS_vol(PP_nVar,:,iElem)
END DO !iElem
IF (warning_linear) THEN
  SWRITE(*,'(A,ES10.2E3)') 'WARNING: during iteration at least one DOF resulted in a phi > phi_max.\n'//&
    '=> Increase Part-RegionElectronRef#-PhiMax if already steady! Phi-Phi_ref=',warning_linear_phi
END IF

!prepare RHS_face ( RHS for lambda system.)
RHS_vol(PP_nVar,:,:)=RHS_vol(PP_nVar,:,:)+JwGP_vol(:,:)*U_out(PP_nVar,:,:)*NonlinVolumeFac(:,:)

RHS_face(PP_nVar,:,:) =0.
DO iElem=1,PP_nElems
  CALL DSYMV('U',nGP_vol,1., InvDhat(:,:,iElem),nGP_vol, &
                             -RHS_vol(PP_nVar,:,iElem),1,0., &
                             rtmp(:),1)
  DO iLocSide=1,6
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    CALL DGEMV('N',nGP_face,nGP_vol,1., &
                        Ehat(:,:,iLocSide,iElem), nGP_face, &
                        rtmp,1,1.,& !add to RHS_face
                        RHS_face(PP_nVar,:,SideID),1)
  END DO
END DO !iElem


DO BCsideID=1,nNeumannBCSides
  SideID=NeumannBC(BCsideID)
  RHS_face(:,:,SideID)=RHS_face(:,:,SideID)+qn_face(:,:,BCSideID)
END DO

#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
CALL Mask_MPISides(PP_nVar,RHS_Face)
#endif /*USE_MPI*/
CALL SmallToBigMortar_HDG(PP_nVar,RHS_face(1:PP_nVar,1:nGP_Face,1:nSides))

#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/


! SOLVE
CALL CheckNonLinRes(RHS_face(1,:,:),lambda(1,:,:),converged,Norm_r2)
IF(PRESENT(ForceCGSolverIteration_opt))THEN
  IF(ForceCGSolverIteration_opt)THEN
    ! Due to the removal of electrons during restart
    converged=.false.
    SWRITE(UNIT_StdOut,*) "Forcing initial CG solver to iterate by setting converged=F (Norm_r2=",Norm_r2,")"
  END IF ! ForceCGSolverIteration_opt
END IF ! ForceCGSolverIteration_opt
IF (converged) THEN
#if defined(IMPA) || defined(ROS)
  IF(DoPrintConvInfo)THEN
    SWRITE(*,*) 'IMPA || ROS - HDGNewton: Newton Iteration has converged in 0 steps...'
  END IF
#else
  SWRITE(*,*) 'HDGNewton: Newton Iteration has converged in 0 steps...'
#endif
ELSE
  CALL CG_solver(RHS_face(PP_nVar,:,:),lambda(PP_nVar,:,:))

  !post processing:
  DO iElem=1,PP_nElems
    ! for post-proc
    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      CALL DGEMV('T',nGP_face,nGP_vol,1., &
                          Ehat(:,:,iLocSide,iElem), nGP_face, &
                          lambda(PP_nVar,:,SideID),1,1.,& !add to RHS_face
                          RHS_vol(PP_nVar,:,iElem),1)
    END DO
    CALL DSYMV('U',nGP_vol,1., InvDhat(:,:,iElem),nGP_vol, &
                               -RHS_vol(PP_nVar,:,iElem),1,0., &
                               U_out(PP_nVar,:,iElem),1)
  END DO !iElem

  IF(NewtonAdaptStartValue) THEN
    IF ((.NOT.DoRestart.AND.ALMOSTEQUAL(time,0.)).OR.(DoRestart.AND.ALMOSTEQUAL(time,RestartTime))) THEN
      DO iElem=1,PP_nElems
        RegionID=ElemToBRRegion(iElem)
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
          r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
          IF (NewtonExactSourceDeriv) THEN
            NonlinVolumeFac(r,iElem) = RegionElectronRef(1,RegionID)/ (RegionElectronRef(3,RegionID)*eps0) &
                         * EXP( (U_out(1,r,iElem)-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
          ELSE
            NonlinVolumeFac(r,iElem)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
          END IF
        END DO; END DO; END DO !i,j,k
      END DO !iElem
      CALL Elem_Mat(td_iter)
      CALL BuildPrecond()
    END IF
  END IF

  converged =.false.
  beLinear=.false.
  AdaptIterNewton = AdaptIterNewtonOld
  DO iter=1,MaxIterNewton

    IF (.NOT.beLinear) THEN
      IF ((iter.EQ.AdaptIterNewtonToLinear)) THEN !.OR.(iter.GT.3*AdaptIterNewtonOld)) THEN
                                                 !removed second cond. to ensure fast convergence with very small AdaptIterNewton
        IF(MPIroot) WRITE(*,*) 'Info: Solver not converging with exact source derivative, switching to linearization (The linear way, baby)'
        DO iElem=1,PP_nElems
          RegionID=ElemToBRRegion(iElem)
          DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
            r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
            NonlinVolumeFac(r,iElem)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
          END DO; END DO; END DO !i,j,k
        END DO !iElem
        CALL Elem_Mat(td_iter)
        CALL BuildPrecond()
        AdaptIterNewton = 0
        beLinear=.true.
      END IF
    END IF

    IF (AdaptIterNewton.GT.0) THEN
      IF (MOD(iter,AdaptIterNewton).EQ.0) THEN
        DO iElem=1,PP_nElems
          RegionID=ElemToBRRegion(iElem)
          DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
            r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
            IF (NewtonExactSourceDeriv) THEN
              NonlinVolumeFac(r,iElem) = RegionElectronRef(1,RegionID)/ (RegionElectronRef(3,RegionID)*eps0) &
                         * EXP( (U_out(1,r,iElem)-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
            ELSE
              NonlinVolumeFac(r,iElem)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
            END IF
          END DO; END DO; END DO !i,j,k
        END DO !iElem
        CALL Elem_Mat(td_iter)
        CALL BuildPrecond()
      END IF
    END IF
    !volume source (volume RHS of u system)
    !SWRITE(*,*) '!!!!!!!!!!!!!!!!!', iter

    warning_linear=.FALSE.
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
        CALL CalcSourceHDG(i,j,k,iElem,RHS_vol(1:PP_nVar,r,iElem),U_out(1,r,iElem),warning_linear,warning_linear_phi)
      END DO; END DO; END DO !i,j,k
      RHS_Vol(PP_nVar,:,iElem)=-JwGP_vol(:,iElem)*RHS_vol(PP_nVar,:,iElem)
    END DO !iElem
    IF (warning_linear) THEN
      !SWRITE(*,'(A,F5.2,A,F5.2,A)') 'HDGNewton WARNING: during iteration at least one DOF resulted in a phi > phi_max. '//&
        !'=> Increase Part-RegionElectronRef#-PhiMax if already steady! Phi-Phi_ref=',warning_linear_phi,'(Phi_ref=',RegionElectronRef(2,RegionID),')'
      SWRITE(*,'(A,ES10.2E3)') 'HDGNewton WARNING: at least one DOF resulted in phi > phi_max. '//&
        'Increase Part-RegionElectronRef#-PhiMax to shift the ref. point! Phi-Phi_ref=',warning_linear_phi!,' (Phi_ref=',RegionElectronRef(2,RegionID),')'
    END IF

    !prepare RHS_face ( RHS for lamdba system.)
    RHS_vol(PP_nVar,:,:)=RHS_vol(PP_nVar,:,:)+JwGP_vol(:,:)*U_out(PP_nVar,:,:)*NonlinVolumeFac(:,:)

    RHS_face(PP_nVar,:,:) =0.
    DO iElem=1,PP_nElems
      CALL DSYMV('U',nGP_vol,1., InvDhat(:,:,iElem),nGP_vol, &
                                 -RHS_vol(PP_nVar,:,iElem),1,0., &
                                 rtmp(:),1)
      DO iLocSide=1,6
        SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
        CALL DGEMV('N',nGP_face,nGP_vol,1., &
                            Ehat(:,:,iLocSide,iElem), nGP_face, &
                            rtmp,1,1.,& !add to RHS_face
                            RHS_face(PP_nVar,:,SideID),1)
      END DO
    END DO !iElem

    !add Neumann
    DO BCsideID=1,nNeumannBCSides
      SideID=NeumannBC(BCsideID)
      RHS_face(:,:,SideID)=RHS_face(:,:,SideID)+qn_face(:,:,BCSideID)
    END DO


#if USE_MPI
  CALL Mask_MPIsides(PP_nVar,RHS_face)
#endif /*USE_MPI*/
  CALL SmallToBigMortar_HDG(PP_nVar,RHS_face(1:PP_nVar,1:nGP_Face,1:nSides))

    ! SOLVE
    CALL CheckNonLinRes(RHS_face(1,:,:),lambda(1,:,:),converged,Norm_r2)
    IF (converged) THEN
#if defined(IMPA) || defined(ROS)
      IF(DoPrintConvInfo)THEN
        SWRITE(*,*) 'HDGNewton: Newton Iteration has converged in ',iter,' steps...'
      END IF
#else
      IF(DoDisplayIter)THEN
        IF(HDGDisplayConvergence.AND.(MOD(td_iter,IterDisplayStep).EQ.0)) THEN
          SWRITE(*,*) 'HDGNewton: Newton Iteration has converged in ',iter,' steps...'
        END IF
      END IF
#endif
      EXIT
    ELSE IF (iter.EQ.MaxIterNewton) THEN
      IPWRITE(UNIT_StdOut,*) "Norm_r2       =", Norm_r2
      IPWRITE(UNIT_StdOut,*) "iter          =", iter
      IPWRITE(UNIT_StdOut,*) "MaxIterNewton =", MaxIterNewton
      CALL abort(__STAMP__,'HDGNewton: Newton Iteration has NOT converged!')
    END IF

    CALL CG_solver(RHS_face(PP_nVar,:,:),lambda(PP_nVar,:,:))

    !post processing:
    DO iElem=1,PP_nElems
      ! for post-proc
      DO iLocSide=1,6
        SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
        CALL DGEMV('T',nGP_face,nGP_vol,1., &
                            Ehat(:,:,iLocSide,iElem), nGP_face, &
                            lambda(PP_nVar,:,SideID),1,1.,& !add to RHS_vol
                            RHS_vol(PP_nVar,:,iElem),1)
      END DO
      CALL DSYMV('U',nGP_vol,1., InvDhat(:,:,iElem),nGP_vol, &
                                 -RHS_vol(PP_nVar,:,iElem),1,0., &
                                 U_out(PP_nVar,:,iElem),1)
    END DO !iElem
  END DO
END IF

#if (PP_nVar==1)
CALL PostProcessGradient(U_out(PP_nVar,:,:),lambda(PP_nVar,:,:),E)
#endif

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE HDGNewton


!===================================================================================================================================
!> Determine the residual of the HDG solution
!===================================================================================================================================
SUBROUTINE CheckNonLinRes(RHS,lambda,converged,Norm_R2)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars  ,ONLY: nGP_face
USE MOD_HDG_Vars  ,ONLY: EpsNonLinear
USE MOD_Mesh_Vars ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars ,ONLY: nMPISides_YOUR
#endif /*USE_MPI*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars  ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)    :: RHS(nGP_face*nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT)    :: lambda(nGP_face*nSides)
LOGICAL, INTENT(INOUT) :: converged
REAL, INTENT(OUT)      :: Norm_r2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(nGP_face*nSides) :: R
INTEGER                         :: VecSize
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
  VecSize=(nSides-nMPIsides_YOUR)*nGP_face
#else
  VecSize=nSides*nGP_face
#endif /*USE_MPI*/
  CALL EvalResidual(RHS,lambda,R)

  CALL VectorDotProduct(VecSize,R(1:VecSize),R(1:VecSize),Norm_R2) !Z=V
!  print*, Norm_R2

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  IF(MPIroot) converged=(Norm_R2.LT.EpsNonLinear**2)
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)
#else
  converged=(Norm_R2.LT.EpsNonLinear**2)
#endif /*USE_MPI*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3)  = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(3) = MPIW8CountField(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/
END SUBROUTINE CheckNonLinRes
#endif /*defined(PARTICLES)*/


!===================================================================================================================================
!> Continuous Gradient solver
!===================================================================================================================================
SUBROUTINE CG_solver(RHS,lambda,iVar)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars           ,ONLY: nGP_face,HDGDisplayConvergence,iteration
USE MOD_HDG_Vars           ,ONLY: EpsCG,MaxIterCG,PrecondType,useRelativeAbortCrit,OutIterCG
USE MOD_TimeDisc_Vars      ,ONLY: iter,IterDisplayStep
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: nMPISides_YOUR
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)    :: RHS(nGP_face*nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT) :: lambda(nGP_face*nSides)
INTEGER, INTENT(IN),OPTIONAL::iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(nGP_face*nSides) :: V,Z,R
REAL                            :: AbortCrit2
REAL                            :: omega,rr,vz,rz1,rz2,Norm_r2
REAL                            :: timestartCG,timeEndCG
INTEGER                         :: VecSize
LOGICAL                         :: converged
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
IF(HDGDisplayConvergence.AND.(MOD(iter,IterDisplayStep).EQ.0)) THEN
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(*,*)'CG solver start'
END IF
TimeStartCG=PICLASTIME()
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
VecSize=(nSides-nMPIsides_YOUR)*nGP_face
#else
VecSize=nSides*nGP_face
#endif /*USE_MPI*/
IF(PRESENT(iVar)) THEN
  CALL EvalResidual(RHS,lambda,R,iVar)
ELSE
  CALL EvalResidual(RHS,lambda,R)
END IF

CALL VectorDotProduct(VecSize,R(1:VecSize),R(1:VecSize),Norm_R2) !Z=V

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

IF(useRelativeAbortCrit)THEN
#if USE_MPI
  IF(MPIroot) converged=(Norm_R2.LT.1e-16)
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)
#else
  converged=(Norm_R2.LT.1e-16)
#endif /*USE_MPI*/
ELSE
#if USE_MPI
  IF(MPIroot) converged=(Norm_R2.LT.EpsCG**2)
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)
#else
  converged=(Norm_R2.LT.EpsCG**2)
#endif /*USE_MPI*/
END IF

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3)  = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(3) = MPIW8CountField(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

IF(converged) THEN !converged
!  SWRITE(*,*)'CG not needed, residual already =0'
!  SWRITE(UNIT_StdOut,'(132("-"))')
    TimeEndCG=PICLASTIME()
    iteration = 0
    IF(MPIroot) CALL DisplayConvergence(TimeEndCG-TimeStartCG, iteration, SQRT(Norm_R2))
  RETURN
END IF !converged
AbortCrit2=EpsCG**2
IF(useRelativeAbortCrit) AbortCrit2=Norm_R2*EpsCG**2

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
IF(PrecondType.NE.0) THEN
  CALL ApplyPrecond(R,V)
ELSE
  V(:)=R(:)
END IF
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
CALL VectorDotProduct(VecSize,R(1:VecSize),V(1:VecSize),rz1) !Z=V

! Conjugate Gradient
!IF(MPIroot) print*, '!!!!!!!!!!!!!!!!!!!!!!'
!IF(MPIroot) print*, iVar
DO iteration=1,MaxIterCG
  ! matrix vector
  IF(PRESENT(iVar)) THEN
    CALL MatVec(V,Z, iVar)
  ELSE
    CALL MatVec(V,Z)
  END IF

  CALL VectorDotProduct(VecSize,V(1:VecSize),Z(1:VecSize),vz)

  omega=rz1/vz

  lambda=lambda+omega*V
  R=R-omega*Z
  CALL VectorDotProduct(VecSize,R(1:VecSize),R(1:VecSize),rr)
  IF(ISNAN(rr)) CALL abort(__STAMP__,'HDG solver residual rr = NaN for CG iteration =', IntInfoOpt=iteration)
#if USE_MPI
  IF(MPIroot) converged=(rr.LT.AbortCrit2)

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3)  = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(3) = MPIW8CountField(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

#else
  converged=(rr.LT.AbortCrit2)
#endif /*USE_MPI*/
  IF(converged) THEN !converged
    TimeEndCG=PICLASTIME()
    CALL EvalResidual(RHS,lambda,R)
    CALL VectorDotProduct(VecSize,R(1:VecSize),R(1:VecSize),Norm_R2) !Z=V (function contains ALLREDUCE)
    IF(MPIroot) CALL DisplayConvergence(TimeEndCG-TimeStartCG, iteration, SQRT(Norm_R2))
    RETURN
  END IF !converged

  IF (MOD(iteration , MAX(INT(REAL(MaxIterCG)/REAL(OutIterCG)),1) ).EQ.0) THEN
    SWRITE(*,'(2(A,I0),2(A,G0))') 'CG solver reached ',iteration, ' of ',MaxIterCG, ' iterations with res = ',rr, ' > ',AbortCrit2
  END IF
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  IF(PrecondType.NE.0) THEN
    CALL ApplyPrecond(R,Z)
  ELSE
    Z(:)=R(:)
  END IF
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
  CALL VectorDotProduct(VecSize,R(1:VecSize),Z(1:VecSize),rz2)
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  V=Z+(rz2/rz1)*V
  rz1=rz2
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
END DO ! iteration
SWRITE(*,*)'CG solver not converged in ',iteration, 'iterations!!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE CG_solver


!===================================================================================================================================
!> Set the global convergence properties of the HDG (CG) Solver and print then to StdOut)
!===================================================================================================================================
SUBROUTINE DisplayConvergence(ElapsedTime, iteration, Norm)
! MODULES
USE MOD_HDG_Vars      ,ONLY: HDGDisplayConvergence,HDGNorm,RunTime,RunTimePerIteration,iterationTotal,RunTimeTotal
USE MOD_Globals       ,ONLY: UNIT_StdOut
USE MOD_TimeDisc_Vars ,ONLY: iter,IterDisplayStep
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: ElapsedTime
INTEGER,INTENT(IN)  :: iteration
REAL,INTENT(IN)     :: Norm
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
RunTime = ElapsedTime
RunTimeTotal = RunTimeTotal + RunTime
IF(iteration.GT.0)THEN
  iterationTotal = iterationTotal + iteration
  RunTimePerIteration = RunTime/REAL(iteration)
ELSE
  RunTimePerIteration = 0.
END IF ! iteration.GT.0
HDGNorm = Norm

IF(HDGDisplayConvergence.AND.(MOD(iter,IterDisplayStep).EQ.0)) THEN
  WRITE(UNIT_StdOut,'(A,1X,I0,A,I0,A)')                '#iterations          :    ',iteration,' (',iterationTotal,' total)'
  WRITE(UNIT_StdOut,'(A,1X,ES25.14E3,A,ES25.14E3,A)')  'RunTime           [s]:',RunTime,' (',RunTimeTotal,' total)'
  WRITE(UNIT_StdOut,'(A,1X,ES25.14E3)')                'RunTime/iteration [s]:',RunTimePerIteration
  !WRITE(UNIT_StdOut,'(A,1X,ES16.7)')'RunTime/iteration/DOF[s]:',(TimeEndCG-TimeStartCG)/REAL(iteration*PP_nElems*nGP_vol)
  WRITE(UNIT_StdOut,'(A,1X,ES25.14E3)')                'Final Residual       :',HDGNorm
  WRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE DisplayConvergence


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE EvalResidual(RHS,lambda,R,iVar)
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars           ,ONLY: nGP_face,nDirichletBCSides,DirichletBC,SetZeroPotentialDOF
USE MOD_Mesh_Vars          ,ONLY: nSides
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)    :: RHS(nGP_face,nSides)
REAL, INTENT(INOUT) :: lambda(nGP_face,nSides)
INTEGER, INTENT(IN),OPTIONAL::iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)   :: R(nGP_face,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: mv(nGP_face,nSides)
INTEGER             :: BCsideID
!===================================================================================================================================
IF(PRESENT(iVar)) THEN
  CALL MatVec(lambda,mv,iVar)
ELSE
  CALL MatVec(lambda,mv)
END IF
R=RHS-mv

!set mv on Dirichlet BC to zero!
#if (PP_nVar!=1)
IF (iVar.EQ.4) THEN
#endif
! TODO direkt als RHS vorgeben! nicht erst hier für PETSC Problem NICHT mehr symmetrisch!
  ! Dirichlet BCs
  DO BCsideID=1,nDirichletBCSides
    R(:,DirichletBC(BCsideID))=0.
  END DO ! SideID=1,nSides

  ! Set residual to zero
  IF(SetZeroPotentialDOF) R(1,1) = 0.

#if (PP_nVar!=1)
END IF
#endif


END SUBROUTINE EvalResidual


!===================================================================================================================================
!> Performs matrix-vector multiplication for lambda system
!>   Parallel Mortar concept:
!>   1) MORTAR, BigToSmall: interpolate lambda from  big to small (small master sides)
!>   2) send lambda from master MPI sides to slave MPI sides (includes small mortar master sides)
!>   3) compute matrix-vector product locally on each proc, in mv array
!>   4) call mask_MPIsides: send  mv contribution from slave MPI sides to master MPI sides and add to master MPI sides
!>   5) MORTAR, SmallToBig: add contribution of finalized small mortar sides to big mortar, via Transpose of interpolation operator
!===================================================================================================================================
SUBROUTINE MatVec(lambda, mv, iVar)
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars          ,ONLY: Smat,nGP_face,nDirichletBCSides,DirichletBC,SetZeroPotentialDOF
USE MOD_Mesh_Vars         ,ONLY: nSides, SideToElem, ElemToSide, nMPIsides_YOUR
USE MOD_FillMortar_HDG    ,ONLY: BigToSmallMortar_HDG,SmallToBigMortar_HDG
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI               ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_HDG_Vars          ,ONLY: Mask_MPIsides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: lambda(nGP_face, nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: mv(nGP_face, nSides)
INTEGER, INTENT(IN),OPTIONAL::iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: firstSideID, lastSideID
INTEGER :: BCsideID,SideID, ElemID, locSideID
INTEGER :: jLocSide,jSideID(6)
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
#if PP_nVar==1
INTEGER           :: dummy
#endif /*PP_nVar==1*/
!===================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
CALL BigToSmallMortar_HDG(1,lambda)
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if USE_MPI
CALL StartReceiveMPIData(1,lambda,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
CALL StartSendMPIData(   1,lambda,1,nSides,SendRequest_U,SendID=1) ! Send MINE
#endif /*USE_MPI*/


#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
firstSideID = 1
lastSideID = nSides-nMPIsides_YOUR

mv=0.

DO SideID=firstSideID,lastSideID
  !master element
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    ElemID    = SideToElem(S2E_ELEM_ID,SideID)
    jSideID(:) = ElemToSide(E2S_SIDE_ID,:,ElemID)
    DO jLocSide = 1,6
      CALL DGEMV('N',nGP_face,nGP_face,1., &
                        Smat(:,:,jLocSide,locSideID,ElemID), nGP_face, &
                        lambda(:,SideID),1,1.,& !add to mv
                        mv(:,jSideID(jLocSide)),1)
    END DO !jLocSide
  END IF !locSideID.NE.-1
  ! neighbour element
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
    jSideID(:)=ElemToSide(E2S_SIDE_ID,:,ElemID)
    DO jLocSide = 1,6
      CALL DGEMV('N',nGP_face,nGP_face,1., &
                        Smat(:,:,jLocSide,locSideID,ElemID), nGP_face, &
                        lambda(:,SideID),1,1.,& !add to mv
                        mv(:,jSideID(jLocSide)),1)
    END DO !jLocSide
  END IF !locSideID.NE.-1
  !add mass matrix
END DO ! SideID=1,nSides
!SWRITE(*,*)'DEBUG---------------------------------------------------------'
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if USE_MPI
! Finish lambda communication
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
firstSideID=nSides-nMPIsides_YOUR+1
lastSideID =nSides
DO SideID=firstSideID,lastSideID
  !master element
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    ElemID    = SideToElem(S2E_ELEM_ID,SideID)
    jSideID(:) = ElemToSide(E2S_SIDE_ID,:,ElemID)
    DO jLocSide = 1,6
      CALL DGEMV('N',nGP_face,nGP_face,1., &
                        Smat(:,:,jLocSide,locSideID,ElemID), nGP_face, &
                        lambda(:,SideID),1,1.,& !add to mv
                        mv(:,jSideID(jLocSide)),1)
    END DO !jLocSide
  END IF !locSideID.NE.-1
  ! neighbour element
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
    jSideID(:)=ElemToSide(E2S_SIDE_ID,:,ElemID)
    DO jLocSide = 1,6
      CALL DGEMV('N',nGP_face,nGP_face,1., &
                        Smat(:,:,jLocSide,locSideID,ElemID), nGP_face, &
                        lambda(:,SideID),1,1.,& !add to mv
                        mv(:,jSideID(jLocSide)),1)
    END DO !jLocSide
  END IF !locSideID.NE.-1
  !add mass matrix
END DO ! SideID=1,nSides
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
CALL Mask_MPIsides(1,mv)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
CALL SmallToBigMortar_HDG(1,mv)

#if (PP_nVar!=1)
IF (iVar.EQ.4) THEN
#endif

  !set mv on Dirichlet BC to zero!
  DO BCsideID=1,nDirichletBCSides
    mv(:,DirichletBC(BCsideID))=0.
  END DO ! SideID=1,nSides

  ! Set potential to zero
  IF(SetZeroPotentialDOF) mv(1,1) = 0.

#if (PP_nVar!=1)
END IF
#endif

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

! Suppress compiler warning
RETURN
#if PP_nVar==1
dummy=iVar
#endif /*PP_nVar==1*/

END SUBROUTINE MatVec


!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProduct(dim1,A,B,Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: dim1
REAL,INTENT(IN)   :: A(dim1)
REAL,INTENT(IN)   :: B(dim1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
#if USE_MPI
REAL              :: ResuSend
#endif
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
Resu=0.
DO i=1,dim1
  Resu=Resu + A(i)*B(i)
END DO
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4)  = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(4) = MPIW8CountField(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE VectorDotProduct


!===================================================================================================================================
!> Apply the block-diagonal preconditioner for the lambda system
!===================================================================================================================================
SUBROUTINE ApplyPrecond(R, V)
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars  ,ONLY: nGP_face, Precond, PrecondType,InvPrecondDiag
USE MOD_HDG_Vars  ,ONLY: MaskedSide
USE MOD_Mesh_Vars ,ONLY: nSides
USE MOD_Mesh_Vars ,ONLY: nMPIsides_YOUR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: R(nGP_face, nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: V(nGP_face, nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: firstSideID, lastSideID, SideID, igf
!===================================================================================================================================
firstSideID = 1
lastSideID = nSides-nMPIsides_YOUR

SELECT CASE(PrecondType)
CASE(0)
  ! do nothing, should not be called
CASE(1) !apply side-block SPD Preconditioner matrix, already Cholesky decomposed
  DO SideID=firstSideID,lastSideID
    IF(MaskedSide(sideID).GT.0) THEN
      V(:,SideID)=0.
    ELSE
      ! solve the preconditioner linear system
      CALL solveSPD(nGP_face,Precond(:,:,SideID),1,R(:,SideID), V(:,SideID))
    END IF !maskedSide
  END DO ! SideID=1,nSides
CASE(2)
  DO SideID=firstSideID,lastSideID
    IF(MaskedSide(sideID).GT.0) THEN
      V(:,SideID)=0.
    ELSE
      ! apply inverse of diagonal preconditioned
      DO igf = 1, nGP_face
        V(igf, SideID) = InvPrecondDiag(igf,SideID)*R(igf,SideID)
      END DO ! igf
    END IF !maskedSide
  END DO ! SideID=1,nSides
END SELECT ! PrecondType
END SUBROUTINE ApplyPrecond



!===================================================================================================================================
!> Solve a symmetrical positive definite linear system of dimension dims
!===================================================================================================================================
SUBROUTINE solveSPD(dimA,A,nRHS,RHS, X)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: dimA
REAL,INTENT(IN)   :: A(dimA, dimA)
INTEGER,INTENT(IN):: nRHS
REAL,INTENT(IN)   :: RHS(dimA,nRHS)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT):: X(dimA,nRHS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: lapack_info
!===================================================================================================================================
X = RHS
CALL DPOTRS('U',dimA,nRHS,A,dimA,X,dimA,lapack_info)
!IF (lapack_info .NE. 0) THEN
!  STOP 'LAPACK ERROR IN SOLVE CHOLESKY!'
!END IF
END SUBROUTINE solveSPD


!===================================================================================================================================
!> During restart, recalculate the gradient of the HDG solution
!===================================================================================================================================
SUBROUTINE RestartHDG(U_out)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_Elem_Mat      ,ONLY: PostProcessGradient
USE MOD_Basis         ,ONLY: getSPDInverse
#if USE_MPI
USE MOD_MPI_Vars
#endif /*USE_MPI*/
#if (PP_nVar==1)
USE MOD_Equation_Vars ,ONLY: E
#elif (PP_nVar==3)
USE MOD_Equation_Vars ,ONLY: B
#else
USE MOD_Equation_Vars ,ONLY: B, E
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(INOUT)  :: U_out(PP_nVar,nGP_vol,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol,PP_nElems)
#endif
!===================================================================================================================================

#if (PP_nVar==1)
  CALL PostProcessGradient(U_out(1,:,:),lambda(1,:,:),E)
#elif (PP_nVar==3)
  DO iVar=1, PP_nVar
    CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  END DO
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
    B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
    B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  END DO; END DO; END DO !i,j,k
#else
  DO iVar=1, 3
    CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  END DO
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
    B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
    B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  END DO; END DO; END DO !i,j,k
  CALL PostProcessGradient(U_out(4,:,:),lambda(4,:,:),E)
#endif
END SUBROUTINE RestartHDG
#endif /*USE_HDG*/


!===================================================================================================================================
!> Finalizes variables necessary for HDG subroutines
!===================================================================================================================================
SUBROUTINE FinalizeHDG()
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars
#if USE_PETSC
USE petsc
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem,nSides,SideToNonUniqueGlobalSide
USE MOD_Mesh_Tools         ,ONLY: LambdaSideToMaster,GetMasteriLocSides
#endif /*USE_LOADBALANCE*/
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
!===================================================================================================================================
HDGInitIsDone = .FALSE.
#if USE_PETSC
PetscCallA(KSPDestroy(ksp,ierr))
PetscCallA(MatDestroy(Smat_petsc,ierr))
PetscCallA(VecDestroy(lambda_petsc,ierr))
PetscCallA(VecDestroy(RHS_petsc,ierr))
PetscCallA(PetscFinalize(ierr))
SDEALLOCATE(PETScGlobal)
SDEALLOCATE(PETScLocalToSideID)
SDEALLOCATE(Smat_BC)
SDEALLOCATE(SmallMortarType)
#endif
SDEALLOCATE(NonlinVolumeFac)
SDEALLOCATE(DirichletBC)
SDEALLOCATE(NeumannBC)
SDEALLOCATE(qn_face)
SDEALLOCATE(qn_face_MagStat)
SDEALLOCATE(delta)
SDEALLOCATE(LL_minus)
SDEALLOCATE(LL_plus)
SDEALLOCATE(Lomega_m)
SDEALLOCATE(Lomega_p)
SDEALLOCATE(Domega)
SDEALLOCATE(InvDhat)
SDEALLOCATE(wGP_vol)
SDEALLOCATE(JwGP_vol)
SDEALLOCATE(Ehat)
SDEALLOCATE(Smat)
SDEALLOCATE(Tau)
SDEALLOCATE(RHS_vol)
SDEALLOCATE(Precond)
SDEALLOCATE(InvPrecondDiag)
SDEALLOCATE(MaskedSide)
SDEALLOCATE(SmallMortarInfo)
SDEALLOCATE(IntMatMortar)

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
    ALLOCATE(lambdaLB(PP_nVar,nGP_face,firstSide:lastSide))
    lambdaLB=0.
  END ASSOCIATE
  IF(nProcessors.GT.1) CALL GetMasteriLocSides()
  DO iSide = 1, nSides
    NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,iSide)

    CALL LambdaSideToMaster(iSide,lambdaLB(:,:,NonUniqueGlobalSideID))
    ! Check if the same global unique side is encountered twice and store both global non-unique side IDs in the array
    ! SideToNonUniqueGlobalSide(1:2,iSide)
    IF(SideToNonUniqueGlobalSide(2,iSide).NE.-1)THEN
      NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(2,iSide)
      CALL LambdaSideToMaster(iSide,lambdaLB(:,:,NonUniqueGlobalSideID))
    END IF ! SideToNonUniqueGlobalSide(1,iSide).NE.-1

  END DO ! iSide = 1, nSides
  IF(nProcessors.GT.1) DEALLOCATE(iLocSides)

END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/

SDEALLOCATE(lambda)
END SUBROUTINE FinalizeHDG


END MODULE MOD_HDG
