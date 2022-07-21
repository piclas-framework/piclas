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
CALL prms%CreateIntOption(    'PrecondType'            ,'Preconditioner type\n 0: no preconditioner\n 1: Side-block SPD preconditioner'&
                                                      //' matrix (already Cholesky decomposed)\n 2: Inverse of diagonal preconditioned', '2')
CALL prms%CreateRealOption(   'epsCG'                  ,'Abort residual of the CG solver', '1.0E-6')
CALL prms%CreateIntOption(    'OutIterCG'              ,'Number of iteration steps between output of CG solver info to std out', '1')

CALL prms%CreateLogicalOption('useRelativeAbortCrit'   ,'Switch between relative and absolute abort criterion', '.FALSE.')
CALL prms%CreateIntOption(    'MaxIterCG'              ,'Maximum number of iterations in the CG solver', '500')
CALL prms%CreateLogicalOption('ExactLambda'            ,'Initially set lambda on all sides (volume and boundaries) via a pre-defined function (ExactFunc)', '.FALSE.')

CALL prms%CreateIntOption(    'HDGSkip'                ,'Number of time step iterations until the HDG solver is called (i.e. all intermediate calls are skipped)', '0')
CALL prms%CreateIntOption(    'HDGSkipInit'            ,'Number of time step iterations until the HDG solver is called (i.e. all intermediate calls are skipped)'&
                                                      //' while time < HDGSkip_t0 (if HDGSkip > 0)', '0')
CALL prms%CreateRealOption(   'HDGSkip_t0'             ,'Time during which HDGSkipInit is used instead of HDGSkip (if HDGSkip > 0)', '0.')

CALL prms%CreateLogicalOption('HDGDisplayConvergence'  ,'Display divergence criteria: Iterations, RunTime and Residual', '.FALSE.')

CALL prms%CreateIntOption(    'HDGZeroPotentialDir'    ,'Direction in which a Dirichlet condition with phi=0 is superimposed on the boundary conditions'&
                                                      //' (1: x, 2: y, 3: z). The default chooses the direction automatically when no other Dirichlet boundary conditions are defined.'&
                                                       ,'-1')

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
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
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
PrecondType          = GETINT('PrecondType')
epsCG                = GETREAL('epsCG')
OutIterCG            = GETINT('OutIterCG')
useRelativeAbortCrit = GETLOGICAL('useRelativeAbortCrit')
MaxIterCG            = GETINT('MaxIterCG')

ExactLambda          = GETLOGICAL('ExactLambda')

ALLOCATE(MaskedSide(1:nSides))
MaskedSide=.FALSE.

IF(nGlobalMortarSides.GT.0)THEN !mortar mesh
  IF(nMortarMPISides.GT.0) CALL abort( &
  __STAMP__,&
  "nMortarMPISides >0: HDG mortar MPI implementation relies on big sides having always only master sides (=> nMortarMPISides=0 )")
END IF !mortarMesh

CALL InitMortar_HDG()

!boundary conditions
nDirichletBCsides=0
nNeumannBCsides  =0
DO SideID=1,nBCSides
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(2,4,5,6) ! Dirichlet
    nDirichletBCsides=nDirichletBCsides+1
  CASE(10,11) ! Neumann
    nNeumannBCsides=nNeumannBCsides+1
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,' unknown BC Type in hdg.f90!',IntInfoOpt=BCType)
  END SELECT ! BCType
END DO

! Check if zero potential must be set on a boundary (or periodic side)
CALL InitZeroPotential()

IF(nDirichletBCsides.GT.0)ALLOCATE(DirichletBC(nDirichletBCsides))
IF(nNeumannBCsides  .GT.0)THEN
  ALLOCATE(NeumannBC(nNeumannBCsides))
  ALLOCATE(qn_face(PP_nVar, nGP_face,nNeumannBCsides))
END IF
#if (PP_nVar!=1)
  IF(nDirichletBCsides.GT.0)ALLOCATE(qn_face_MagStat(PP_nVar, nGP_face,nDirichletBCsides))
#endif
nDirichletBCsides=0
nNeumannBCsides  =0
DO SideID=1,nBCSides
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(2,4,5,6) ! Dirichlet
    nDirichletBCsides=nDirichletBCsides+1
    DirichletBC(nDirichletBCsides)=SideID
    MaskedSide(SideID)=.TRUE.
  CASE(10,11) !Neumann,
    nNeumannBCsides=nNeumannBCsides+1
    NeumannBC(nNeumannBCsides)=SideID
  END SELECT ! BCType
END DO

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


!stabilization parameter
ALLOCATE(Tau(PP_nElems))
DO iElem=1,PP_nElems
  Tau(iElem)=2./((SUM(JwGP_vol(:,iElem)))**(1./3.))  !1/h ~ 1/vol^(1/3) (volref=8)
END DO !iElem

IF(.NOT.DoSwapMesh)THEN ! can take very long, not needed for swap mesh run as only the state file is converted
  CALL Elem_Mat(0_8) ! takes iter=0 (kind=8)
END IF


CALL BuildPrecond()

ALLOCATE(lambda(PP_nVar,nGP_face,nSides))
lambda=0.
ALLOCATE(RHS_vol(PP_nVar, nGP_vol,PP_nElems))
RHS_vol=0.

HDGInitIsDone = .TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT HDG DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitHDG


!===================================================================================================================================
!> Check if any Dirichlet BCs are present (globally, not only on the local processor).
!> If there are none, an arbitrary potential is set at one of the boundaries to ensure
!> convergence of the HDG solver. This is required for setups where fully periodic and/or Neumann boundaries are solely used.
!> This only works for Cartesian meshes, i.e., that the boundary faces must be perpendicular to two of the three Cartesian axes
!===================================================================================================================================
SUBROUTINE InitZeroPotential()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_HDG_Vars    ,ONLY: ZeroPotentialSideID,HDGZeroPotentialDir
USE MOD_Mesh        ,ONLY: GetMeshMinMaxBoundaries
USE MOD_Mesh_Vars   ,ONLY: nBCs,BoundaryType,nSides,BC,xyzMinMax,NGeo,Face_xGP
USE MOD_ReadInTools ,ONLY: PrintOption,GETINT
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: SideID,BCAlpha,BCType,BCState,iBC,nZeroPotentialSides,nZeroPotentialSidesGlobal,ZeroPotentialDir
INTEGER           :: nZeroPotentialSidesMax,BCSide
REAL,DIMENSION(3) :: x,v1,v2,v3
REAL              :: norm,I(3,3)
!===================================================================================================================================
! Initialize variables
HDGZeroPotentialDir = GETINT('HDGZeroPotentialDir')
ZeroPotentialSideID = -1
I(:,1)              = (/1. , 0. , 0./)
I(:,2)              = (/0. , 1. , 0./)
I(:,3)              = (/0. , 0. , 1./)

! Every processor has to check every BC
DO iBC=1,nBCs
  BCType  = BoundaryType(iBC,BC_TYPE)  ! 1
  BCState = BoundaryType(iBC,BC_STATE) ! 2
  SELECT CASE(BCType)
  CASE(1) ! periodic
    ! do nothing
  CASE(2,4,5,6) ! Dirichlet
    ZeroPotentialSideID = 0 ! no zero potential required
    EXIT ! as soon as one Dirichlet BC is found, no zero potential must be used
  CASE(10,11) ! Neumann
    ! do nothing
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,' unknown BC Type in hdg.f90!',IntInfoOpt=BCType)
  END SELECT ! BCType
END DO

! If a Dirichlet BC is found ZeroPotentialSideID is zero and the following is skipped
IF(ZeroPotentialSideID.EQ.-1)THEN
  ! Check if the user has selected a specific direction
  IF(HDGZeroPotentialDir.EQ.-1)THEN
    ! Select the direction (x, y or z), which has the largest extent (to account for 1D and 2D setups for example)
    CALL GetMeshMinMaxBoundaries()

    ! Calc max extents in each direction for comparison
    x(1) = xyzMinMax(2)-xyzMinMax(1)
    x(2) = xyzMinMax(4)-xyzMinMax(3)
    x(3) = xyzMinMax(6)-xyzMinMax(5)
    ZeroPotentialDir=MAXLOC(x,DIM=1)
  ELSE
    ZeroPotentialDir = HDGZeroPotentialDir
  END IF ! HDGZeroPotentialDir.EQ.-1
  CALL PrintOption('Zero potential side activated in direction (1: x, 2: y, 3: z)','OUTPUT',IntOpt=ZeroPotentialDir)

  nZeroPotentialSides = 0 ! Initialize
  DO SideID=1,nSides ! Periodic sides are not within the 1,nBCSides list !
    IF(MAXVAL(ABS(Face_xGP(:,:,:,SideID))).LE.0.) CYCLE ! slave sides
    BCSide=BC(SideID)
    IF(BCSide.EQ.0) CYCLE ! inner sides
    BCType =BoundaryType(BCSide,BC_TYPE)
    BCState=BoundaryType(BCSide,BC_STATE)
    BCAlpha=BoundaryType(BCSide,BC_ALPHA)
    IF(BCType.EQ.0) CYCLE ! skip inner sides

    ! Check if the normal vector of the face points in the direction (or negative direction) of the ZeroPotentialDir (tolerance 1e-5)
    v1(:) = Face_xGP(1:3 , NGeo , 0    , SideID) - Face_xGP(1:3 , 0 , 0 , SideID)
    v2(:) = Face_xGP(1:3 , 0    , NGeo , SideID) - Face_xGP(1:3 , 0 , 0 , SideID)
    v3(:) = CROSSNORM(v1,v2)
    norm = ABS(DOT_PRODUCT(I(:,ZeroPotentialDir),v3))
    IF(ALMOSTEQUALRELATIVE(norm, 1.0, 1E-5))THEN
      ZeroPotentialSideID = SideID
      nZeroPotentialSides = nZeroPotentialSides + 1
    END IF ! ALMOSTEQUALRELATIVE(norm, 1.0, 1E-5)
  END DO

#if USE_MPI
  ! Combine number of found zero potential sides to make sure that at least one is found
  IF(MPIroot)THEN
    CALL MPI_REDUCE(nZeroPotentialSides , nZeroPotentialSidesGlobal , 1 , MPI_INTEGER , MPI_SUM , 0 , MPI_COMM_WORLD , IERROR)
    CALL MPI_REDUCE(nZeroPotentialSides , nZeroPotentialSidesMax    , 1 , MPI_INTEGER , MPI_MAX , 0 , MPI_COMM_WORLD , IERROR)
  ELSE
    CALL MPI_REDUCE(nZeroPotentialSides , 0                         , 1 , MPI_INTEGER , MPI_SUM , 0 , MPI_COMM_WORLD , IERROR)
    CALL MPI_REDUCE(nZeroPotentialSides , 0                         , 1 , MPI_INTEGER , MPI_MAX , 0 , MPI_COMM_WORLD , IERROR)
  END IF
#else
  nZeroPotentialSidesGlobal = nZeroPotentialSides
#endif /*USE_MPI*/
  LBWRITE(UNIT_StdOut,'(A,I0)') " Found (global) number of zero potential sides: ", nZeroPotentialSidesMax

  ! Sanity checks for root
  IF(MPIroot)THEN
    ! 1) multiples sides found
    IF(nZeroPotentialSidesMax.GT.1)THEN
      LBWRITE(UNIT_StdOut,'(A)') " WARNING: Found more than 1 zero potential side on a proc and currently, only one can be considered."
      LBWRITE(UNIT_StdOut,'(A,I0,A)') " WARNING: nZeroPotentialSidesGlobal: ", nZeroPotentialSidesMax, " (may lead to problems)"
    END IF

    ! 2) no sides found
    IF(nZeroPotentialSidesGlobal.EQ.0)THEN
      WRITE(UNIT_StdOut,*) " Sanity check: this fails when the mesh is not Cartesian. This needs to be implemented if required."
      CALL abort(__STAMP__,'Setup has no Dirichlet BCs and no zero potential sides where found.')
    END IF
  END IF
END IF ! ZeroPotentialSideID.EQ.0

END SUBROUTINE InitZeroPotential


!===================================================================================================================================
!> HDG solver for linear or non-linear systems
!===================================================================================================================================
SUBROUTINE HDG(t,U_out,iter,ForceCGSolverIteration_opt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars ,ONLY: iStage
#endif
#if (USE_HDG && (PP_nVar==1))
USE MOD_TimeDisc_Vars ,ONLY: dt,dt_Min
USE MOD_Equation_Vars ,ONLY: E,Et
USE MOD_Globals_Vars  ,ONLY: eps0
USE MOD_Analyze_Vars  ,ONLY: CalcElectricTimeDerivative
#endif /*(USE_HDG && (PP_nVar==1))*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: t !time
INTEGER(KIND=8),INTENT(IN)  :: iter
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only required for BR electron fluid)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_out(PP_nVar,nGP_vol,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(PARTICLES)
LOGICAL :: ForceCGSolverIteration_loc
#endif /*defined(PARTICLES)*/
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(4))
#endif /*EXTRAE*/

! Calculate temporal derivate of E in last iteration before Analyze_dt is reached: Store E^n here
#if (USE_HDG && (PP_nVar==1))
IF(CalcElectricTimeDerivative)THEN
  IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.ALMOSTEQUAL(dt,dt_Min(DT_END)))THEN
    Et(:,:,:,:,:) = E(:,:,:,:,:)
  END IF
END IF ! CalcElectricTimeDerivative
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

! Run the appropriate HDG solver: Newton or Linear
#if defined(PARTICLES)
IF(UseBRElectronFluid) THEN
  IF (HDGNonLinSolver.EQ.1) THEN
    ForceCGSolverIteration_loc = MERGE(ForceCGSolverIteration_opt, .FALSE., PRESENT(ForceCGSolverIteration_opt))
    CALL HDGNewton(t, U_out, iter, ForceCGSolverIteration_loc)
  ELSE
    CALL abort(__STAMP__,'Defined HDGNonLinSolver not implemented (HDGFixPoint has been removed!) HDGNonLinSolver = ',&
    IntInfoOpt=HDGNonLinSolver)
  END IF
ELSE
#endif /*defined(PARTICLES)*/
  CALL HDGLinear(t,U_out)
#if defined(PARTICLES)
END IF
#endif /*defined(PARTICLES)*/

! Calculate temporal derivate of E in last iteration before Analyze_dt is reached: Store E^n+1 here and calculate the derivative
#if (USE_HDG && (PP_nVar==1))
IF(CalcElectricTimeDerivative)THEN
  IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.ALMOSTEQUAL(dt,dt_Min(DT_END)))THEN
    Et(:,:,:,:,:) = eps0*(E(:,:,:,:,:)-Et(:,:,:,:,:)) / dt
  END IF
END IF ! CalcElectricTimeDerivative
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
USE MOD_Equation               ,ONLY: CalcSourceHDG,ExactFunc
USE MOD_Equation_Vars          ,ONLY: IniExactFunc
USE MOD_Equation_Vars          ,ONLY: chitens_face
USE MOD_Mesh_Vars              ,ONLY: Face_xGP,BoundaryType,nSides,BC
USE MOD_Mesh_Vars              ,ONLY: ElemToSide,NormVec,SurfElem
USE MOD_Interpolation_Vars     ,ONLY: wGP
USE MOD_Elem_Mat               ,ONLY: PostProcessGradient
USE MOD_FillMortar_HDG         ,ONLY: SmallToBigMortar_HDG
#if (PP_nVar==1)
USE MOD_Equation_Vars          ,ONLY: E
#elif (PP_nVar==3)
USE MOD_Equation_Vars          ,ONLY: B
#else
USE MOD_Equation_Vars          ,ONLY: B, E
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
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
    CASE(5) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional)
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-1,Face_xGP(:,p,q,SideID),lambda(iVar,r:r,SideID),t=time,iRefState=BCState)
      END DO; END DO !p,q
    CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional)
      DO q=0,PP_N; DO p=0,PP_N
        r=q*(PP_N+1) + p+1
        CALL ExactFunc(-2,Face_xGP(:,p,q,SideID),lambda(iVar,r:r,SideID),t=time,iRefState=BCState)
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

  ! Check if zero potential sides are present
  IF(ZeroPotentialSideID.GT.0) lambda(iVar,:,ZeroPotentialSideID) = ZeroPotentialValue
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

  CALL CG_solver(RHS_face(iVar,:,:),lambda(iVar,:,:),iVar)
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
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only BR electron fluid)
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
  CALL abort(__STAMP__,'Nonlinear Newton solver only available with EQ-system Poisson!')
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
  CASE(5) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-1,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,iRefState=BCState)
    END DO; END DO !p,q
  CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-2,Face_xGP(:,p,q,SideID),lambda(PP_nVar,r:r,SideID),t=time,iRefState=BCState)
    END DO; END DO !p,q
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
      CALL abort(&
        __STAMP__&
        ,'HDGNewton: Newton Iteration has NOT converged!')
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
USE MOD_Mesh_Vars ,ONLY: nSides,nMPISides_YOUR
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars  ,ONLY: MPIW8TimeField
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
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,iError)
#else
  converged=(Norm_R2.LT.EpsNonLinear**2)
#endif /*USE_MPI*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3) = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
END SUBROUTINE CheckNonLinRes


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
USE MOD_Mesh_Vars          ,ONLY: nSides,nMPISides_YOUR
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField
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
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,iError)
#else
  converged=(Norm_R2.LT.1e-16)
#endif /*USE_MPI*/
ELSE
#if USE_MPI
  IF(MPIroot) converged=(Norm_R2.LT.EpsCG**2)
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,iError)
#else
  converged=(Norm_R2.LT.EpsCG**2)
#endif /*USE_MPI*/
END IF

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3) = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
#endif /*defined(MEASURE_MPI_WAIT)*/

IF(converged) THEN !converged
!  SWRITE(*,*)'CG not needed, residual already =0'
!  SWRITE(UNIT_StdOut,'(132("-"))')
    TimeEndCG=PICLASTIME()
    iteration = 0
    IF(MPIroot) CALL DisplayConvergence(TimeEndCG-TimeStartCG, iteration, Norm_R2)
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

  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,iError)

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3) = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
#endif /*defined(MEASURE_MPI_WAIT)*/

#else
  converged=(rr.LT.AbortCrit2)
#endif /*USE_MPI*/
  IF(converged) THEN !converged
    TimeEndCG=PICLASTIME()
    CALL EvalResidual(RHS,lambda,R)
    CALL VectorDotProduct(VecSize,R(1:VecSize),R(1:VecSize),Norm_R2) !Z=V (function contains ALLREDUCE)
    IF(MPIroot) CALL DisplayConvergence(TimeEndCG-TimeStartCG, iteration, Norm_R2)
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
SUBROUTINE DisplayConvergence(ElapsedTime, iteration, Norm_R2)
! MODULES
USE MOD_HDG_Vars      ,ONLY: HDGDisplayConvergence,HDGNorm,RunTime,RunTimePerIteration,iterationTotal,RunTimeTotal
USE MOD_Globals       ,ONLY: UNIT_StdOut
USE MOD_TimeDisc_Vars ,ONLY: iter,IterDisplayStep
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: ElapsedTime
INTEGER,INTENT(IN)  :: iteration
REAL,INTENT(IN)     :: Norm_R2
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
HDGNorm = SQRT(Norm_R2)

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
USE MOD_HDG_Vars           ,ONLY: nGP_face,nDirichletBCSides,DirichletBC,ZeroPotentialSideID,ZeroPotentialValue
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

  ! Dirichlet BCs
  DO BCsideID=1,nDirichletBCSides
    R(:,DirichletBC(BCsideID))=0.
  END DO ! SideID=1,nSides

  ! Set potential to zero
  IF(ZeroPotentialSideID.GT.0) R(:,ZeroPotentialSideID)= ZeroPotentialValue

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
USE MOD_HDG_Vars          ,ONLY: Smat,nGP_face,nDirichletBCSides,DirichletBC,ZeroPotentialSideID,ZeroPotentialValue
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
  IF(ZeroPotentialSideID.GT.0) mv(:,ZeroPotentialSideID) = ZeroPotentialValue

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
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField
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
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4) = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
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
    IF(MaskedSide(sideID)) THEN
      V(:,SideID)=0.
    ELSE
      ! solve the preconditioner linear system
      CALL solveSPD(nGP_face,Precond(:,:,SideID),1,R(:,SideID), V(:,SideID))
    END IF !maskedSide
  END DO ! SideID=1,nSides
CASE(2)
  DO SideID=firstSideID,lastSideID
    IF(MaskedSide(sideID)) THEN
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
USE MOD_HDG_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
HDGInitIsDone = .FALSE.
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
SDEALLOCATE(lambda)
SDEALLOCATE(RHS_vol)
SDEALLOCATE(Precond)
SDEALLOCATE(InvPrecondDiag)
SDEALLOCATE(MaskedSide)
SDEALLOCATE(SmallMortarInfo)
SDEALLOCATE(IntMatMortar)
END SUBROUTINE FinalizeHDG


END MODULE MOD_HDG
