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
!> Contains global variables used by the HDG modules.
!===================================================================================================================================
MODULE MOD_HDG_Vars
! MODULES
#if USE_MPI
USE MOD_Globals
#endif /*USE_MPI*/
#if USE_PETSC
USE PETSc
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
INTEGER             :: nGP_vol                !< =(PP_N+1)**3
INTEGER             :: nGP_face               !< =(PP_N+1)**2

#if USE_PETSC
Mat                 :: Smat_petsc
Vec                 :: RHS_petsc
Vec                 :: lambda_petsc
KSP                 :: ksp
Vec                 :: lambda_local_petsc
VecScatter          :: scatter_petsc
IS                  :: idx_local_petsc
IS                  :: idx_global_petsc
Vec                 :: lambda_local_conductors_petsc
VecScatter          :: scatter_conductors_petsc
IS                  :: idx_local_conductors_petsc
IS                  :: idx_global_conductors_petsc
INTEGER,ALLOCATABLE :: PETScGlobal(:)         !< PETScGlobal(SideID) maps the local SideID to global PETScSideID
INTEGER,ALLOCATABLE :: PETScLocalToSideID(:)  !< PETScLocalToSideID(PETScLocalSideID) maps the local PETSc side to SideID
REAL,ALLOCATABLE    :: Smat_BC(:,:,:,:)       !< side to side matrix for dirichlet (D) BCs, (ngpface,ngpface,6Sides,DSides)
INTEGER             :: nPETScSides            !< nSides - nDirichletSides
INTEGER             :: nPETScUniqueSides      !< nPETScSides - nMPISides_YOUR
INTEGER             :: nPETScUniqueSidesGlobal
#endif
LOGICAL             :: useHDG=.FALSE.
LOGICAL             :: ExactLambda =.FALSE.   !< Flag to initialize exact function for lambda
REAL,ALLOCATABLE    :: InvDhat(:,:,:)         !< Inverse of Dhat matrix (nGP_vol,nGP_vol,nElems)
REAL,ALLOCATABLE    :: Ehat(:,:,:,:)          !< Ehat matrix (nGP_Face,nGP_vol,6sides,nElems)
REAL,ALLOCATABLE    :: wGP_vol(:)             !< 3D quadrature weights
REAL,ALLOCATABLE    :: JwGP_vol(:,:)          !< 3D quadrature weights*Jacobian for all elements
REAL,ALLOCATABLE    :: lambda(:,:,:)          !< lambda, ((PP_N+1)^2,nSides)
REAL,ALLOCATABLE    :: lambdaLB(:,:,:)        !< lambda, ((PP_N+1)^2,nSides)
REAL,ALLOCATABLE    :: iLocSides(:,:,:)       !< iLocSides, ((PP_N+1)^2,nSides) - used for I/O and ALLGATHERV of lambda
REAL,ALLOCATABLE    :: RHS_vol(:,:,:)         !< Source RHS
REAL,ALLOCATABLE    :: Tau(:)                 !< Stabilization parameter, per element
REAL,ALLOCATABLE    :: Smat(:,:,:,:,:)        !< side to side matrix, (ngpface, ngpface, 6sides, 6sides, nElems)
REAL,ALLOCATABLE    :: Precond(:,:,:)         !< block diagonal preconditioner for lambda(nGP_face, nGP-face, nSides)
REAL,ALLOCATABLE    :: InvPrecondDiag(:,:)    !< 1/diagonal of Precond
REAL,ALLOCATABLE    :: qn_face(:,:,:)         !< for Neumann BC
REAL,ALLOCATABLE    :: qn_face_MagStat(:,:,:) !< for Neumann BC
INTEGER             :: nDirichletBCsides
INTEGER             :: nNeumannBCsides
INTEGER             :: nConductorBCsides      !< Number of processor-local sides that are conductors (FPC) in [1:nBCSides]
LOGICAL             :: SetZeroPotentialDOF    !< Flag to set a single DOF, if only periodic and Neumann boundaries are present
INTEGER,ALLOCATABLE :: ConductorBC(:)
INTEGER,ALLOCATABLE :: DirichletBC(:)
INTEGER,ALLOCATABLE :: NeumannBC(:)
LOGICAL             :: HDGnonlinear           !< Use non-linear sources for HDG? (e.g. Boltzmann electrons)
LOGICAL             :: NewtonExactSourceDeriv
LOGICAL             :: NewtonAdaptStartValue
INTEGER             :: AdaptIterNewton
INTEGER             :: AdaptIterNewtonToLinear
INTEGER             :: AdaptIterNewtonOld
INTEGER             :: HDGNonLinSolver        !< 1 Newton, 2 Fixpoint
REAL,ALLOCATABLE    :: NonlinVolumeFac(:,:)   !< Factor for Volumeintegration necessary for nonlinear sources
!mappings
INTEGER             :: sideDir(6),pm(6),dirPm2iSide(2,3)
REAL,ALLOCATABLE    :: delta(:,:)
REAL,ALLOCATABLE    :: LL_minus(:,:),LL_plus(:,:)
REAL,ALLOCATABLE    :: Domega(:,:)
REAL,ALLOCATABLE    :: Lomega_m(:),Lomega_p(:)
!CG parameters
INTEGER             :: PrecondType=0          !< 0: none 1: block diagonal 2: only diagonal 3:Identity, debug
INTEGER             :: MaxIterCG
INTEGER             :: MaxIterNewton
INTEGER             :: OutIterCG
REAL                :: EpsCG,EpsNonLinear
LOGICAL             :: UseRelativeAbortCrit
LOGICAL             :: HDGInitIsDone=.FALSE.
INTEGER             :: HDGSkip, HDGSkipInit
REAL                :: HDGSkip_t0
INTEGER,ALLOCATABLE :: MaskedSide(:)          !< 1:nSides: all sides which are set to zero in matvec
!mortar variables
REAL,ALLOCATABLE    :: IntMatMortar(:,:,:,:)  !< Interpolation matrix for mortar: (nGP_face,nGP_Face,1:4(iMortar),1:3(MortarType))
INTEGER,ALLOCATABLE :: SmallMortarInfo(:)     !< 1:nSides: info on small Mortar sides:
                                              !< -1: is neighbor small mortar , 0: not a small mortar, 1: small mortar on big side
#if USE_PETSC
INTEGER,ALLOCATABLE :: SmallMortarType(:,:)   !< Type of Mortar side ([1] Type, [2] Side, nSides)
                                              !< [1] Type: mortar type this small side belongs to (1-3)
                                              !< [2] Side: Small side number (1-4)
#endif
LOGICAL             :: HDGDisplayConvergence  !< Display divergence criteria: Iterations, Runtime and Residual
REAL                :: RunTime                !< CG Solver runtime
REAL                :: RunTimePerIteration    !< CG Solver runtime per iteration
REAL                :: HDGNorm                !< Norm
INTEGER             :: iteration              !< number of iterations to achieve the norm
INTEGER             :: iterationTotal         !< number of iterations over the course of a time step (possibly multiple stages)
REAL                :: RunTimeTotal           !< CG Solver runtime sum over the course of a time step (possibly multiple stages)

! --- Boltzmann relation (BR) electron fluid
LOGICAL               :: UseBRElectronFluid            !< Indicates usage of BR electron fluid model
INTEGER               :: BRNbrOfRegions                !< Nbr of regions to be mapped to Elems
LOGICAL               :: CalcBRVariableElectronTemp    !< Use variable ref. electron temperature for BR electron fluid
CHARACTER(255)        :: BRVariableElectronTemp        !< Variable electron reference temperature when using Boltzmann relation
                                                       !< electron model (default is using a constant temperature)
REAL                  :: BRVariableElectronTempValue   !< Final electron temperature
INTEGER, ALLOCATABLE  :: ElemToBRRegion(:)             !< ElemToBRRegion(1:nElems)
REAL, ALLOCATABLE     :: BRRegionBounds(:,:)           !< BRRegionBounds ((xmin,xmax,ymin,...)|1:BRNbrOfRegions)
REAL, ALLOCATABLE     :: RegionElectronRef(:,:)        !< RegionElectronRef((rho0[C/m^3],phi0[V],Te[eV])|1:BRNbrOfRegions)
REAL, ALLOCATABLE     :: RegionElectronRefBackup(:,:)  !< RegionElectronRefBackup(rho0[C/m^3],phi0[V],Te[eV])|1:BRNbrOfRegions) when using variable
                                                       !< reference electron temperature
REAL                  :: BRTimeStepMultiplier          !< Factor that is multiplied with the ManualTimeStep when using BR model
REAL                  :: BRTimeStepBackup              !< Original time step
LOGICAL               :: BRAutomaticElectronRef        !< Automatically obtain the reference parameters (from a fully kinetic
                                                       !< simulation), store them in .h5 state and in .csv
INTEGER               :: nBRAverageElems               !< Processor local number of elements in which the reference values are averaged
INTEGER               :: nBRAverageElemsGlobal         !< Global number of elements in which the reference values are averaged
INTEGER, ALLOCATABLE  :: BRAverageElemToElem(:)        !< Mapping BR average elem to processo-local elem
#if defined(PARTICLES)
! --- Switching between BR and fully kinetic HDG
LOGICAL               :: BRConvertElectronsToFluid     !< User variable for removing all electrons and using BR instead
REAL                  :: BRConvertElectronsToFluidTime !< Time when kinetic electrons should be converted to BR fluid electrons
LOGICAL               :: BRConvertFluidToElectrons     !< User variable for creating particles from BR electron fluid (uses
REAL                  :: BRConvertFluidToElectronsTime !< Time when BR fluid electrons should be converted to kinetic electrons
INTEGER               :: BRConvertMode                 !< Mode used for switching BR->kin->BR OR kin->BR->kin
                                                       !< and ElectronDensityCell ElectronTemperatureCell from .h5 state file)
LOGICAL               :: BRConvertModelRepeatedly      !< Repeat the switch between BR and kinetic multiple times
LOGICAL               :: BRElectronsRemoved            !< True if electrons were removed during restart (only BR electrons)
REAL                  :: DeltaTimeBRWindow             !< Time length when BR is active (possibly multiple times)
LOGICAL               :: BRNullCollisionDefault        !< Flag (backup of read-in parameter) whether null collision method
                                                       !< (determining number of pairs based on maximum relaxation frequency) is used
#endif /*defined(PARTICLES)*/

! --- Sub-communicator groups

#if USE_MPI
TYPE tMPIGROUP
  INTEGER                     :: ID                     !< MPI communicator ID
  INTEGER                     :: UNICATOR=MPI_COMM_NULL !< MPI communicator for floating boundary condition
  INTEGER                     :: nProcs                 !< number of MPI processes part of the FPC group
  INTEGER                     :: nProcsWithSides        !< number of MPI processes part of the FPC group and actual FPC sides
  INTEGER                     :: MyRank                 !< MyRank within communicator
END TYPE
#endif /*USE_MPI*/

!===================================================================================================================================
!-- Floating boundary condition
!===================================================================================================================================

LOGICAL                       :: UseFPC             !< Automatic flag when FPCs are active

TYPE tFPC
  REAL,ALLOCATABLE            :: Voltage(:)         !< Electric potential on floating boundary condition for each (required) BC index over all processors. This is the value that is reduced to the MPI root process
  REAL,ALLOCATABLE            :: VoltageProc(:)     !< Electric potential on floating boundary condition for each (required) BC index for a single processor. This value is non-zero only when the processor has an actual FPC side
  REAL,ALLOCATABLE            :: Charge(:)          !< Accumulated charge on floating boundary condition for each (required) BC index over all processors
  REAL,ALLOCATABLE            :: ChargeProc(:)      !< Accumulated charge on floating boundary condition for each (required) BC index for a single processor
#if USE_MPI
  TYPE(tMPIGROUP),ALLOCATABLE :: COMM(:)            !< communicator and ID for parallel execution
#endif /*USE_MPI*/
  !INTEGER                     :: NBoundaries       !< Total number of boundaries where the floating boundary condition is evaluated
  INTEGER                     :: nFPCBounds         !< Global number of boundaries that are FPC with BCType=20 in [1:nBCs],
!                                                   !< they might belong to the same group (electrically connected)
  INTEGER                     :: nUniqueFPCBounds   !< Global number of independent FPC after grouping certain BC sides together
!                                                   !< (electrically connected) with the same BCState ID
  INTEGER,ALLOCATABLE         :: BCState(:)         !< BCState of the i-th FPC index
  !INTEGER,ALLOCATABLE         :: BCIDToFPCBCID(:)  !< Mapping BCID to FPC BCID (1:nPartBound)
  INTEGER,ALLOCATABLE         :: Group(:,:)         !< FPC%Group(1:FPC%nFPCBounds,3)
                                                    !<   1: BCState
                                                    !<   2: iUniqueFPC (i-th FPC group ID)
                                                    !<   3: number of BCSides for each FPC group
  INTEGER,ALLOCATABLE         :: GroupGlobal(:)     !< Sum of nSides associated with each i-th FPC boundary
END TYPE

TYPE(tFPC)   :: FPC
!===================================================================================================================================
!-- Electric Potential Condition (for decharging)
!===================================================================================================================================

LOGICAL                       :: UseEPC             !< Automatic flag when EPCs are active

TYPE tEPC
  REAL,ALLOCATABLE            :: Voltage(:)         !< Electric potential on floating boundary condition for each (required) BC index over all processors. This is the value that is reduced to the MPI root process
  REAL,ALLOCATABLE            :: VoltageProc(:)     !< Electric potential on floating boundary condition for each (required) BC index for a single processor. This value is non-zero only when the processor has an actual EPC side
  REAL,ALLOCATABLE            :: Charge(:)          !< Accumulated charge on floating boundary condition for each (required) BC index over all processors
  REAL,ALLOCATABLE            :: ChargeProc(:)      !< Accumulated charge on floating boundary condition for each (required) BC index for a single processor
  REAL,ALLOCATABLE            :: Resistance(:)      !< Vector (length corresponds to the number of EPC boundaries) with the resistance for each EPC in Ohm
#if USE_MPI
  TYPE(tMPIGROUP),ALLOCATABLE :: COMM(:)            !< communicator and ID for parallel execution
#endif /*USE_MPI*/
  !INTEGER                     :: NBoundaries       !< Total number of boundaries where the floating boundary condition is evaluated
  INTEGER                     :: nEPCBounds         !< Global number of boundaries that are EPC with BCType=20 in [1:nBCs],
!                                                   !< they might belong to the same group (electrically connected)
  INTEGER                     :: nUniqueEPCBounds   !< Global number of independent EPC after grouping certain BC sides together
!                                                   !< (electrically connected) with the same BCState ID
  INTEGER,ALLOCATABLE         :: BCState(:)         !< BCState of the i-th EPC index
  !INTEGER,ALLOCATABLE         :: BCIDToEPCBCID(:)  !< Mapping BCID to EPC BCID (1:nPartBound)
  INTEGER,ALLOCATABLE         :: Group(:,:)         !< EPC%Group(1:EPC%nEPCBounds,3)
                                                    !<   1: BCState
                                                    !<   2: iUniqueEPC (i-th EPC group ID)
                                                    !<   3: number of BCSides for each EPC group
  INTEGER,ALLOCATABLE         :: GroupGlobal(:)     !< Sum of nSides associated with each i-th EPC boundary
END TYPE

TYPE(tEPC)   :: EPC
#if defined(PARTICLES)
!===================================================================================================================================
!-- Coupled Power Potential (CPP)
!-- Special BC with floating potential that is defined by the absorbed power of the charged particles
!===================================================================================================================================

LOGICAL           :: UseCoupledPowerPotential !< Switch calculation on/off
INTEGER,PARAMETER :: CPPDataLength=6          !< Number of variables in BVData

#if USE_MPI
TYPE(tMPIGROUP) :: CPPCOMM       !< communicator and ID for parallel execution
#endif /*USE_MPI*/

REAL    :: CoupledPowerPotential(CPPDataLength) !< (/min, start, max/) electric potential, e.g., used at all BoundaryType = (/2,2/)
REAL    :: CoupledPowerTarget                   !< Target input power at all BoundaryType = (/2,2/)
REAL    :: CoupledPowerRelaxFac                 !< Relaxation factor for calculation of new electric potential
REAL    :: CoupledPowerFrequency                !< Frequency with which the integrated power is calculated (must be consistent Part-AnalyzeStep, i.e., that one cycle with period T=1/f must be larger than Part-AnalyzeStep * dt)
INTEGER :: CoupledPowerMode                     !< Method for power adjustment with 1: instantaneous power, 2: moving average power, 3: integrated power

!===================================================================================================================================
!-- Bias Voltage (for calculating a BC voltage bias for certain BCs)
!===================================================================================================================================

LOGICAL           :: UseBiasVoltage !< Automatic flag when bias voltage is to be used
INTEGER,PARAMETER :: BVDataLength=3 !< Number of variables in BVData

TYPE tBV
#if USE_MPI
  TYPE(tMPIGROUP)     :: COMM                 !< communicator and ID for parallel execution
#endif /*USE_MPI*/
  INTEGER             :: NPartBoundaries      !< Total number of boundaries where the particles are counted
  INTEGER,ALLOCATABLE :: PartBoundaries(:)    !< Part-boundary number on which the particles are counted
  REAL                :: Frequency            !< Adaption nrequency with which the bias voltage is adjusted (every period T = 1/f the bias voltage is changed)
  REAL                :: Delta                !< Voltage difference used to change the current bias voltage (may also be adjusted over time automatically)
  REAL                :: BVData(BVDataLength) !< 1: bias voltage
!                                             !< 2: Ion excess
!                                             !< 3: sim. time when next adjustment happens
END TYPE

TYPE(tBV)   :: BiasVoltage
#endif /*defined(PARTICLES)*/
!===================================================================================================================================

#if USE_MPI
!no interface for reshape inout vector
!INTERFACE Mask_MPIsides
!  MODULE PROCEDURE Mask_MPIsides
!END INTERFACE

PUBLIC :: Mask_MPIsides
#endif /*USE_MPI*/

CONTAINS

#if USE_MPI
!===================================================================================================================================
!> communicate contribution from MPI slave sides to MPI master sides  and set slaves them to zero afterwards.
!===================================================================================================================================
SUBROUTINE Mask_MPIsides(firstdim,v)
! MODULES
USE MOD_MPI_Vars
USE MOD_Mesh_Vars      ,ONLY: nSides
USE MOD_Mesh_Vars      ,ONLY: nMPIsides_YOUR,nMPIsides,nMPIsides_MINE
USE MOD_MPI            ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN   ) :: firstdim
REAL   ,INTENT(INOUT) :: v(firstdim,nGP_face, nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: vbuf(firstdim,nGP_Face,nMPISides_MINE)
INTEGER :: startbuf,endbuf
!===================================================================================================================================

startbuf=nSides-nMPISides+1
endbuf=nSides-nMPISides+nMPISides_MINE
IF(nMPIsides_MINE.GT.0)vbuf=v(:,:,startbuf:endbuf)
CALL StartReceiveMPIData(firstdim,v,1,nSides,RecRequest_U ,SendID=2)  ! Receive MINE
CALL StartSendMPIData(   firstdim,v,1,nSides,SendRequest_U,SendID=2)  ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) ! Send YOUR - receive MINE
IF(nMPIsides_MINE.GT.0) v(:,:,startbuf:endbuf)=v(:,:,startbuf:endbuf)+vbuf
IF(nMPIsides_YOUR.GT.0) v(:,:,nSides-nMPIsides_YOUR+1:nSides)=0. !set send buffer to zero!

END SUBROUTINE Mask_MPIsides
#endif /*USE_MPI*/


#endif /*USE_HDG*/
END MODULE MOD_HDG_Vars
