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
USE MOD_Basis                 ,ONLY: PolynomialDerivativeMatrix
USE MOD_Interpolation_Vars    ,ONLY: N_Inter,NMax
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_Elem_Mat              ,ONLY: Elem_Mat,BuildPrecond
USE MOD_ReadInTools           ,ONLY: GETLOGICAL,GETREAL,GETINT
USE MOD_Mesh_Vars             ,ONLY: nBCSides,N_SurfMesh
USE MOD_Mesh_Vars             ,ONLY: BoundaryType,nSides,BC
USE MOD_Mesh_Vars             ,ONLY: nGlobalMortarSides,nMortarMPISides,N_VolMesh
USE MOD_Mesh_Vars             ,ONLY: DoSwapMesh, offSetElem
USE MOD_Basis                 ,ONLY: InitializeVandermonde,LegendreGaussNodesAndWeights,BarycentricWeights
USE MOD_FillMortar_HDG        ,ONLY: InitMortar_HDG
USE MOD_HDG_Vars              ,ONLY: BRNbrOfRegions,ElemToBRRegion,RegionElectronRef
USE MOD_DG_Vars               ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave
#if defined(PARTICLES)
USE MOD_Part_BR_Elecron_Fluid ,ONLY: UpdateNonlinVolumeFac
USE MOD_Restart_Vars          ,ONLY: DoRestart
#endif /*defined(PARTICLES)*/
#if USE_PETSC
USE PETSc
USE MOD_Mesh_Vars             ,ONLY: nMPISides_YOUR
#if USE_MPI
USE MOD_MPI                   ,ONLY: StartReceiveMPIDataInt,StartSendMPIDataInt,FinishExchangeMPIData
#endif /*USE_MPI*/
USE MOD_Mesh_Vars             ,ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars             ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars             ,ONLY: ElemToSide
#endif /*USE_PETSC*/
#if USE_MPI
USE MOD_MPI                   ,ONLY: StartReceiveMPISurfDataType, StartSendMPISurfDataType, FinishExchangeMPISurfDataType
USE MOD_MPI_Vars
#endif
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
INTEGER           :: i,j,k,r,iElem,SideID,Nloc,iNeumannBCsides,NSideMin,NSideMax, iSide
INTEGER           :: BCType,BCState
REAL              :: D(0:Nmax,0:Nmax)
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
INTEGER             :: iLocSide

INTEGER             :: nLocalPETScDOFs
INTEGER             :: iLocalPETScDOF,iDOF
INTEGER             :: OffsetCounter
INTEGER             :: PETScDOFOffsetsMPI(nProcessors)
INTEGER,ALLOCATABLE :: localToGlobalPETScDOF(:)
#endif
!#if USE_MPI
REAL              :: tmp(3,0:Nmax,0:Nmax)
!#endif
!===================================================================================================================================
IF(HDGInitIsDone)THEN
   LBWRITE(*,*) "InitHDG already called."
   RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT HDG...'

HDGDisplayConvergence = GETLOGICAL('HDGDisplayConvergence')

ALLOCATE(nGP_vol(1:NMax))
ALLOCATE(nGP_face(1:NMax))
DO Nloc = 1, NMax
  nGP_vol(Nloc)  = (Nloc+1)**3
  nGP_face(Nloc) = (Nloc+1)**2
END DO ! Nloc = 1, NMax

HDGSkip = GETINT('HDGSkip')
IF (HDGSkip.GT.0) THEN
  HDGSkipInit = GETINT('HDGSkipInit')
  HDGSkip_t0  = GETREAL('HDGSkip_t0')
ELSE
  HDGSkip=0
END IF

HDGNonLinSolver = -1 ! init

#if USE_MPI
ALLOCATE(SurfExchange(nNbProcs))
DO iNbProc=1,nNbProcs
  ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv(MAXVAL(DataSizeSurfRecMax(iNbProc,:))))
  ALLOCATE(SurfExchange(iNbProc)%SurfDataSend(MAXVAL(DataSizeSurfSendMax(iNbProc,:))))
END DO !iProc=1,nNBProcs
CALL StartReceiveMPISurfDataType(RecRequest_Geo, 1, 1)
CALL StartSendMPISurfDataType(SendRequest_Geo,1,1)
CALL FinishExchangeMPISurfDataType(SendRequest_Geo,RecRequest_Geo,1, 1)
DO iNbProc=1,nNbProcs
  DEALLOCATE(SurfExchange(iNbProc)%SurfDataRecv)
  DEALLOCATE(SurfExchange(iNbProc)%SurfDataSend)
  ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv(MAXVAL(DataSizeSurfRecMin( iNbProc,:))))
  ALLOCATE(SurfExchange(iNbProc)%SurfDataSend(MAXVAL(DataSizeSurfSendMin(iNbProc,:))))
  ALLOCATE(SurfExchange(iNbProc)%SurfDataRecv2(MAXVAL(DataSizeSurfRecMin( iNbProc,:))**2))
  ALLOCATE(SurfExchange(iNbProc)%SurfDataSend2(MAXVAL(DataSizeSurfSendMin(iNbProc,:))**2))
END DO !iProc=1,nNBProcs
#endif /*USE_MPI*/

! Build SurfElemMin for all sides (including Mortar sides)
DO iSide = 1, nSides
  ! Get SurfElemMin
  NSideMax = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
  NSideMin = N_SurfMesh(iSide)%NSideMin
  IF(NSideMax.EQ.NSideMin)THEN
    N_SurfMesh(iSide)%SurfElemMin(:,:) = N_SurfMesh(iSide)%SurfElem(:,:)
  ELSE
    ! From high to low
    ! Transform the slave side to the same degree as the master: switch to Legendre basis
    CALL ChangeBasis2D(1, NSideMax, NSideMax, N_Inter(NSideMax)%sVdm_Leg, N_SurfMesh(iSide)%SurfElem(0:NSideMax,0:NSideMax), tmp(1,0:NSideMax,0:NSideMax))
     !Switch back to nodal basis
    CALL ChangeBasis2D(1, NSideMin, NSideMin, N_Inter(NSideMin)%Vdm_Leg , tmp(1,0:NSideMin,0:NSideMin)                     , N_SurfMesh(iSide)%SurfElemMin(0:NSideMin,0:NSideMin))
  END IF ! NSideMax.EQ.NSideMin
END DO ! iSide = 1, nSides

#if USE_PETSC
! initialize PETSc stuff!
PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
#endif

! Initialize element containers
ALLOCATE(HDG_Vol_N(1:PP_nElems))

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

  CALL InitMortar_HDG()
END IF !mortarMesh

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
  CASE(10) ! Neumann
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
  CASE(10) !Neumann,
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

ALLOCATE(HDG_Surf_N(1:nSides))

IF(nNeumannBCsides.GT.0)THEN
  DO iNeumannBCsides = 1, nNeumannBCsides
    SideID = NeumannBC(iNeumannBCsides)
    Nloc = DG_Elems_master(SideID)
    ALLOCATE(HDG_Surf_N(SideID)%qn_face(PP_nVar, nGP_face(Nloc)))
  END DO ! iNeumannBCsides = 1, nNeumannBCsides
END IF

#if USE_PETSC
! Create PETSc Mappings
OffsetPETScSide=0

! TODO PETSC P-Adaption - MORTARS
! Count all Mortar slave sides and remove them from PETSc vector
! TODO How to compute those
nMortarMasterSides = 0
!DO SideID=1,nSides
!  IF(SmallMortarInfo(SideID).EQ.1) THEN
!    nMortarMasterSides = nMortarMasterSides + 1
!  END IF
!END DO
nPETScUniqueSides = nSides-nDirichletBCSides-nMPISides_YOUR-nMortarMasterSides-nConductorBCsides

! TODO PETSC P-Adaption - MPI (Should work?) We probably do not need nPETScUniqueSidesGlobal
#if USE_MPI

CALL MPI_ALLGATHER(nPETScUniqueSides,1,MPI_INTEGER,OffsetPETScSideMPI,1,MPI_INTEGER,MPI_COMM_PICLAS,IERROR)
DO iProc=1, myrank
  OffsetPETScSide = OffsetPETScSide + OffsetPETScSideMPI(iProc)
END DO
nPETScUniqueSidesGlobal = SUM(OffsetPETScSideMPI) + FPC%nUniqueFPCBounds
#endif

! TODO PETSC P-Adaption - Improvement: Delete PETScGlobal
! PETScGlobal stuff
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
! TODO PETSC P-Adaption - MORTARS
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

! PETSc P-Adaption Mappings (with MPI) -----------------------------------------------------

! 1. Calculate nLocalPETScDOFs without nMPISides_YOUR
nLocalPETScDOFs = 0
DO SideID=1,nSides-nMPISides_YOUR
  IF(PETScGlobal(SideID).LT.0) CYCLE ! TODO is this MaskedSide?
  nLocalPETScDOFs = nLocalPETScDOFs + nGP_face(N_SurfMesh(SideID)%NSideMin)
END DO

! 2. Get Offsets
OffsetCounter = 0
#if USE_MPI
! TODO Where to define OffsetCounter and PETScDOFOffsetsMPI
CALL MPI_ALLGATHER(nLocalPETScDOFs,1,MPI_INTEGER,PETScDOFOffsetsMPI,1,MPI_INTEGER,MPI_COMM_PICLAS,IERROR)
DO iProc=1, myrank
  OffsetCounter = OffsetCounter + PETScDOFOffsetsMPI(iProc)
END DO
! 3. Calculate nGlobalPETScDOFs from the offsets
nGlobalPETScDOFs = SUM(PETScDOFOffsetsMPI)
#else
nGlobalPETScDOFs = nLocalPETScDOFs
#endif

#if USE_MPI
! 4. Add the nMPISides_YOUR to nLocalPETScDOFs
DO SideID=nSides-nMPISides_YOUR+1,nSides
  IF(PETScGlobal(SideID).LT.0) CYCLE ! TODO is this MaskedSide?
  nLocalPETScDOFs = nLocalPETScDOFs + nGP_face(N_SurfMesh(SideID)%NSideMin)
END DO
#endif

! 5. Calculate the global offset for each side
ALLOCATE(OffsetGlobalPETScDOF(nSides))
DO SideID=1,nSides-nMPISides_YOUR
  IF(PETScGlobal(SideID).LT.0) CYCLE ! TODO is this MaskedSide?
  OffsetGlobalPETScDOF(SideID) = OffsetCounter
  OffsetCounter = OffsetCounter + nGP_face(N_SurfMesh(SideID)%NSideMin)
END DO

#if USE_MPI
! 6. Communicate OffsetGlobalPETScDOF to fill YOUR sides
CALL StartReceiveMPIDataInt(1,OffsetGlobalPETScDOF,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
CALL StartSendMPIDataInt(   1,OffsetGlobalPETScDOF,1,nSides,SendRequest_U,SendID=1) ! Send MINE
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif

! 7. Create LocalToGlobal mapping (we now know the mapping)
! TODO Which mappings start at 0?
! TODO where to introduce localToGlobalPETScDOF
ALLOCATE(localToGlobalPETScDOF(nLocalPETScDOFs))
iLocalPETScDOF = 0
DO SideID=1,nSides
  IF(PETScGlobal(SideID).LT.0) CYCLE ! TODO is this MaskedSide?
  DO iDOF=1,nGP_face(N_SurfMesh(SideID)%NSideMin)
    iLocalPETScDOF = iLocalPETScDOF + 1
    LocalToGlobalPETScDOF(iLocalPETScDOF) = OffsetGlobalPETScDOF(SideID) + iDOF - 1
  END DO
END DO
! ------------------------------------------------------
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

!ALLOCATE(delta(0:PP_N,0:PP_N))
!delta=0.
!DO i=0,PP_N
  !delta(i,i)=1.
!END DO !i

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



#if USE_PETSC
PetscCallA(MatCreate(PETSC_COMM_WORLD,PETScSystemMatrix,ierr))
PetscCallA(MatSetSizes(PETScSystemMatrix,PETSC_DECIDE,PETSC_DECIDE,nGlobalPETScDOFs,nGlobalPETScDOFs,ierr))
PetscCallA(MatSetType(PETScSystemMatrix,MATAIJ,ierr)) ! Sparse (mpi) matrix, TODO P-Adaption Symmetricity is set later!

! TODO PETSC P-Adaption - FPC Preallocation (also with FPCs)
PetscCallA(MatSetUp(PETScSystemMatrix, ierr))
#endif

!stabilization parameter
ALLOCATE(Tau(PP_nElems))
DO iElem=1,PP_nElems
  Tau(iElem)=2./((SUM(HDG_Vol_N(iElem)%JwGP_vol(:)))**(1./3.))  !1/h ~ 1/vol^(1/3) (volref=8)
END DO !iElem

IF(.NOT.DoSwapMesh)THEN ! can take very long, not needed for swap mesh run as only the state file is converted
  CALL Elem_Mat(0_8) ! takes iter=0 (kind=8)
END IF

#if USE_PETSC
PetscCallA(KSPCreate(PETSC_COMM_WORLD,PETScSolver,ierr))
PetscCallA(KSPSetOperators(PETScSolver,PETScSystemMatrix,PETScSystemMatrix,ierr))

IF(PrecondType.GE.10) THEN
  PetscCallA(KSPSetType(PETScSolver,KSPPREONLY,ierr)) ! Exact solver
ELSE
  PetscCallA(KSPSetType(PETScSolver,KSPCG,ierr)) ! CG solver for sparse symmetric positive definite matrix

  PetscCallA(KSPSetInitialGuessNonzero(PETScSolver,PETSC_TRUE, ierr))

  PetscCallA(KSPSetNormType(PETScSolver, KSP_NORM_UNPRECONDITIONED, ierr))
  PetscCallA(KSPSetTolerances(PETScSolver,1.E-20,epsCG,PETSC_DEFAULT_REAL,MaxIterCG,ierr))
END IF
#endif

DO iElem = 1, PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(HDG_Vol_N(iElem)%RHS_vol(PP_nVar, nGP_vol(Nloc)))
  HDG_Vol_N(iElem)%RHS_vol=0.
END DO ! iElem = 1, PP_nElems

DO SideID = 1, nSides
  NSideMin = N_SurfMesh(SideID)%NSideMin
  ALLOCATE(HDG_Surf_N(SideID)%lambda(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%lambda=0.
  NSideMax = MAX(DG_Elems_master(SideID),DG_Elems_slave(SideID))
  ALLOCATE(HDG_Surf_N(SideID)%lambdaMax(PP_nVar,nGP_face(NSideMax)))
  HDG_Surf_N(SideID)%lambdaMax=0.
  ALLOCATE(HDG_Surf_N(SideID)%RHS_face(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%RHS_face=0.
  ALLOCATE(HDG_Surf_N(SideID)%mv(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%mv=0.
  ALLOCATE(HDG_Surf_N(SideID)%R(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%R=0.
  ALLOCATE(HDG_Surf_N(SideID)%V(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%V=0.
  ALLOCATE(HDG_Surf_N(SideID)%Z(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%Z=0.
#if USE_MPI
  ALLOCATE(HDG_Surf_N(SideID)%buf(PP_nVar,nGP_face(NSideMin)))
  HDG_Surf_N(SideID)%buf=0.
#if USE_PETSC
! TODO PETSC P-Adaption - ?
#else
  IF(PrecondType.EQ.1)THEN
    ALLOCATE(HDG_Surf_N(SideID)%buf2(nGP_face(NSideMin),nGP_face(NSideMin)))
    HDG_Surf_N(SideID)%buf2=0.
  END IF ! PrecondType.EQ.1
#endif /*USE_PETSC*/
#endif
END DO ! SideID = 1, nSides

! Requires HDG_Surf_N(SideID)%buf
CALL BuildPrecond()

#if USE_PETSC
! allocate RHS & lambda vectors
PetscCallA(VecCreate(PETSC_COMM_WORLD,PETScSolution,ierr))
PetscCallA(VecSetSizes(PETScSolution,PETSC_DECIDE,nGlobalPETScDOFs,ierr))
PetscCallA(VecSetType(PETScSolution,VECSTANDARD,ierr))
PetscCallA(VecSetUp(PETScSolution,ierr))
PetscCallA(VecDuplicate(PETScSolution,PETScRHS,ierr))


! Scatter Context
! Create scatter context to access local values from global petsc vector
! Create a local vector for all local DOFs
PetscCallA(VecCreateSeq(PETSC_COMM_SELF,nLocalPETScDOFs,PETScSolutionLocal,ierr))
! Create a PETSc Vector 0:(nLocalPETScDOFs-1)
PetscCallA(ISCreateStride(PETSC_COMM_SELF,nLocalPETScDOFs,0,1,idx_local_petsc,ierr))
! Create a PETSc Vector of the Global DOF IDs
PetscCallA(ISCreateGeneral(PETSC_COMM_WORLD,nLocalPETScDOFs,localToGlobalPETScDOF,PETSC_COPY_VALUES,idx_global_petsc,ierr))
! Create a scatter context to extract the local dofs
PetscCallA(VecScatterCreate(PETScSolution,idx_global_petsc,PETScSolutionLocal,idx_local_petsc,PETScScatter,ierr))

! TODO PETSC P-Adaption - FPC
!IF(UseFPC)THEN
!  PetscCallA(VecCreateSeq(PETSC_COMM_SELF,nGP_face*FPC%nUniqueFPCBounds,lambda_local_conductors_petsc,ierr))
!  PetscCallA(ISCreateStride(PETSC_COMM_SELF,nGP_face*FPC%nUniqueFPCBounds,0,1,idx_local_conductors_petsc,ierr))
!  ALLOCATE(indx(FPC%nUniqueFPCBounds))
!  DO i=1,FPC%nUniqueFPCBounds
!    indx(i) = nPETScUniqueSidesGlobal-FPC%nUniqueFPCBounds+i-1
!  END DO
!  PetscCallA(ISCreateBlock(PETSC_COMM_WORLD,nGP_face,FPC%nUniqueFPCBounds,indx,PETSC_COPY_VALUES,idx_global_conductors_petsc,ierr))
!  DEALLOCATE(indx)
!  PetscCallA(VecScatterCreate(PETScSolution,idx_global_conductors_petsc,lambda_local_conductors_petsc,idx_local_conductors_petsc,scatter_conductors_petsc,ierr))
!END IF

DEALLOCATE(localToGlobalPETScDOF)
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
CASE(800,900,901,1000,1100) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite
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
SUBROUTINE HDG(t,iter,ForceCGSolverIteration_opt)
#else
SUBROUTINE HDG(t,iter)
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
USE MOD_TimeDisc_Vars   ,ONLY: dt
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
REAL,INTENT(IN)     :: t !time
INTEGER(KIND=8),INTENT(IN)  :: iter
#if defined(PARTICLES)
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only required for BR electron fluid)
#endif /*defined(PARTICLES)*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(PARTICLES)
LOGICAL :: ForceCGSolverIteration_loc
INTEGER :: iUniqueEPCBC
#endif /*defined(PARTICLES)*/
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(4))
#endif /*EXTRAE*/

#if (USE_HDG && (PP_nVar==1))
! Calculate temporal derivate of D in last iteration before Analyze_dt is reached: Store E^n here
IF(CalcElectricTimeDerivative) CALL CalculateElectricTimeDerivative(iter,1)
#endif /*(USE_HDG && (PP_nVar==1))*/

! Check whether the solver should be skipped in this iteration
IF (iter.GT.0 .AND. HDGSkip.NE.0) THEN
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
IF(CalcElectricTimeDerivative) CALL CalculateElectricTimeDerivative(iter,2)
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
#if (USE_HDG && (PP_nVar==1))
INTEGER           :: iDir,iElem
#endif /*(USE_HDG && (PP_nVar==1))*/
!===================================================================================================================================

#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
IF((iStage.NE.1).AND.(iStage.NE.nRKStages)) RETURN
#endif

! iter is incremented after this function and then checked in analyze routine with iter+1
IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.ALMOSTEQUAL(dt,dt_Min(DT_END)).OR.(MOD(iter+1,FieldAnalyzeStep).EQ.0))THEN
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

    ! Calculate the electric VDL surface potential from the particle and electric displacement current
    IF(DoVirtualDielectricLayer) CALL CalculatePhiAndEFieldFromCurrentsVDL()
    WRITE (*,*) "iter =", iter, "yes"
  END IF ! mode.EQ.1
ELSE
  WRITE (*,*) "iter =", iter, "no"
END IF

END SUBROUTINE CalculateElectricTimeDerivative


!===================================================================================================================================
!> description
!===================================================================================================================================
SUBROUTINE CalculatePhiAndEFieldFromCurrentsVDL()
! MODULES
USE MOD_Globals                ,ONLY: VECNORM
USE MOD_Globals_Vars           ,ONLY: eps0
USE MOD_TimeDisc_Vars          ,ONLY: dt
USE MOD_Mesh_Vars              ,ONLY: N_SurfMesh,SideToElem,nBCSides,N_SurfMesh,offSetElem,BC
USE MOD_DG_Vars                ,ONLY: U_N,N_DG_Mapping
USE MOD_PICDepo_Vars           ,ONLY: PS_N
USE MOD_Particle_Boundary_Vars ,ONLY: N_SurfVDL,PartBound,ElementThicknessVDL
USE MOD_ProlongToFace          ,ONLY: ProlongToFace_Side
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ElemID,SideID,ilocSide,Nloc,iPartBound,p,q
REAL,ALLOCATABLE :: jface(:,:,:),Dtface(:,:,:)
REAL             :: coeff
!===================================================================================================================================
! 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!     interpolate the vector field E = (/Ex, Ey, Ez/) to the boundary face
DO SideID=1,nBCSides
  ! Get the local element index
  ElemID = SideToElem(S2E_ELEM_ID,SideID)
  ! Get local polynomial degree of the element
  Nloc   = N_DG_Mapping(2,ElemID+offSetElem)
  ! Get particle boundary index
  iPartBound = PartBound%MapToPartBC(BC(SideID))

  ! Skip sides that are not a VDL boundary (these sides are still in the list of sides)
  IF(PartBound%ThicknessVDL(iPartBound).GT.0.0)THEN

    ! Allocate jface and Dtface depending on the local polynomial degree
    ALLOCATE(jface( 1:3,0:Nloc,0:Nloc))
    ALLOCATE(Dtface(1:3,0:Nloc,0:Nloc))
    ! Get local side index
    ilocSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
    ! Prolong-to-face depending on orientation in reference element
    CALL ProlongToFace_Side(3, Nloc, ilocSide, 0, U_N(ElemID)%Dt(:,:,:,:), Dtface)
    CALL ProlongToFace_Side(3, Nloc, ilocSide, 0, PS_N(ElemID)%PartSource(1:3,:,:,:), jface)

    ! 2.) Apply the normal vector to get the normal electric field
    DO q=0,Nloc
      DO p=0,Nloc
        ASSOCIATE( normal => -N_SurfMesh(SideID)%NormVec(1:3,p,q), D => Dtface(1:3,p,q), j => jface(1:3,p,q))
          ! D_normal =  <D,normal>
          Dtface(1,p,q) = DOT_PRODUCT(D,normal)
          ! j_normal =  <j,normal>
          jface(1,p,q)  = DOT_PRODUCT(j,normal)
        END ASSOCIATE
      END DO ! p
    END DO ! q

    ! Calculate coefficient
    coeff = (PartBound%ThicknessVDL(iPartBound)*dt)/(eps0*PartBound%PermittivityVDL(iPartBound))

    ! Calculate the corrected E-field
    !N_SurfVDL(SideID)%U(2:4,:,:) = N_SurfVDL(SideID)%U(2:4,:,:) * (ElementThicknessVDL(ElemID)/PartBound%ThicknessVDL(iPartBound))

    ! Get Phi_F
    DO q=0,Nloc
      DO p=0,Nloc
        ASSOCIATE(      E => N_SurfVDL(SideID)%U(10:12,p,q)     ,&
                   normal => N_SurfMesh(SideID)%NormVec(1:3,p,q),&
                       jp => jface(1,p,q)                       ,&
                       jD => Dtface(1,p,q)                      ,&
                     PhiF => N_SurfVDL(SideID)%U(9,p,q)         )
          ! Normal vector points outwards on BC sides, hence, invert it
          !Edir = DOT_PRODUCT(E,-normal)

          ! Reconstruct Phi_F from the current density and the (uncorrected) electric displacement fields in each element
          ! PhiF_From_Currents: PhiF^n - PhiF^n-1 = -(d*dt)/(eps0*epsR)*(jp+jD)
          PhiF = PhiF - coeff*(jp+jD)

          ! Reconstruct E from Phi_Max via E = Phi/d
          E =  PhiF/PartBound%ThicknessVDL(iPartBound)*normal

        END ASSOCIATE
      END DO ! p
    END DO ! q
    !WRITE (*,*) "N_SurfVDL(SideID)%U(9,:,:) =", N_SurfVDL(SideID)%U(9,:,:)

    DEALLOCATE(jface)
    DEALLOCATE(Dtface)

  END IF ! PartBound%ThicknessVDL(iPartBound).GT.0.0
END DO ! SideID=1,nBCSides

END SUBROUTINE CalculatePhiAndEFieldFromCurrentsVDL
#endif /*(USE_HDG && (PP_nVar==1))*/


!===================================================================================================================================
!> During restart, recalculate the gradient of the HDG solution
!===================================================================================================================================
SUBROUTINE RestartHDG()
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
USE MOD_Interpolation_Vars ,ONLY: NMax
#if USE_PETSC
USE petsc
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
PetscCallA(VecDestroy(PETScRHS,ierr))
PetscCallA(PetscFinalize(ierr))
SDEALLOCATE(PETScGlobal)
SDEALLOCATE(PETScLocalToSideID)
SDEALLOCATE(Smat_BC)
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
    lambdaLB=0.
  END ASSOCIATE
  IF(nProcessors.GT.1) CALL GetMasteriLocSides()
  DO iSide = 1, nSides
    NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,iSide)

    CALL LambdaSideToMaster(1,iSide,lambdaLB(:,:,NonUniqueGlobalSideID),N_SurfMesh(iSide)%NSideMin)
    ! Check if the same global unique side is encountered twice and store both global non-unique side IDs in the array
    ! SideToNonUniqueGlobalSide(1:2,iSide)
    IF(SideToNonUniqueGlobalSide(2,iSide).NE.-1)THEN
      NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(2,iSide)
      CALL LambdaSideToMaster(1,iSide,lambdaLB(:,:,NonUniqueGlobalSideID),N_SurfMesh(iSide)%NSideMin)
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


END MODULE MOD_HDG
