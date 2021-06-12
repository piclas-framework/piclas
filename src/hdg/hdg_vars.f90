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
MODULE MOD_HDG_Vars
!===================================================================================================================================
! Contains global variables used by the HDG modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
INTEGER             :: nGP_vol              !=(PP_N+1)**3
INTEGER             :: nGP_face             !=(PP_N+1)**2

LOGICAL             :: useHDG=.FALSE.
LOGICAL             :: ExactLambda =.FALSE. ! Flag to initialize exact function for lambda
REAL,ALLOCATABLE    :: InvDhat(:,:,:)       ! Inverse of Dhat matrix (nGP_vol,nGP_vol,nElems)
REAL,ALLOCATABLE    :: Ehat(:,:,:,:)        ! Ehat matrix (nGP_Face,nGP_vol,6sides,nElems)
REAL,ALLOCATABLE    :: wGP_vol(:)           ! 3D quadrature weights
REAL,ALLOCATABLE    :: JwGP_vol(:,:)        ! 3D quadrature weights*Jacobian for all elements
REAL,ALLOCATABLE    :: lambda(:,:,:)          ! lambda, ((PP_N+1)^2,nSides)
REAL,ALLOCATABLE    :: RHS_vol(:,:,:)         ! Source RHS
REAL,ALLOCATABLE    :: Tau(:)               ! Stabilization parameter, per element
REAL,ALLOCATABLE    :: Smat(:,:,:,:,:)      ! side to side matrix, (ngpface, ngpface, 6sides, 6sides, nElems)
REAL,ALLOCATABLE    :: Precond(:,:,:)       ! block diagonal preconditioner for lambda(nGP_face, nGP-face, nSides)
REAL,ALLOCATABLE    :: InvPrecondDiag(:,:)  ! 1/diagonal of Precond
REAL,ALLOCATABLE    :: qn_face(:,:,:)         ! for Neumann BC
REAL,ALLOCATABLE    :: qn_face_MagStat(:,:,:)         ! for Neumann BC
INTEGER             :: nDirichletBCsides
INTEGER             :: ZeroPotentialSideID   ! SideID, where the solution is set zero to enforce convergence
REAL,PARAMETER      :: ZeroPotentialValue=0. ! This can be set to an arbitrary value (in the range of the potential solution)
INTEGER             :: nNeumannBCsides
INTEGER,ALLOCATABLE :: DirichletBC(:)
INTEGER,ALLOCATABLE :: NeumannBC(:)
LOGICAL             :: HDGnonlinear            ! Use non-linear sources for HDG? (e.g. Boltzmann electrons)
LOGICAL             :: NewtonExactSourceDeriv
LOGICAL             :: NewtonAdaptStartValue
INTEGER             :: AdaptIterNewton
INTEGER             :: AdaptIterNewtonToLinear
INTEGER             :: AdaptIterNewtonOld
INTEGER             :: HDGNonLinSolver  ! 1 Newton, 2 Fixpoint
REAL,ALLOCATABLE    :: NonlinVolumeFac(:,:)      !Factor for Volumeintegration necessary for nonlinear sources
!mappings
INTEGER             :: sideDir(6),pm(6),dirPm2iSide(2,3)
REAL,ALLOCATABLE    :: delta(:,:)
REAL,ALLOCATABLE    :: LL_minus(:,:),LL_plus(:,:)
REAL,ALLOCATABLE    :: Domega(:,:)
REAL,ALLOCATABLE    :: Lomega_m(:),Lomega_p(:)
!CG parameters
INTEGER             :: PrecondType=0  !0: none 1: block diagonal 2: only diagonal 3:Identity, debug
INTEGER             :: MaxIterCG, MaxIterNewton, OutIterCG
REAL                :: EpsCG,EpsNonLinear
LOGICAL             :: UseRelativeAbortCrit
LOGICAL             :: HDGInitIsDone=.FALSE.
INTEGER             :: HDGSkip, HDGSkipInit
REAL                :: HDGSkip_t0
LOGICAL,ALLOCATABLE :: MaskedSide(:)      ! 1:nSides: all sides which are set to zero in matvec
!mortar variables
REAL,ALLOCATABLE    :: IntMatMortar(:,:,:,:) ! Interpolation matrix for mortar: (nGP_face,nGP_Face,1:4(iMortar),1:3(MortarType))
INTEGER,ALLOCATABLE :: SmallMortarInfo(:)      ! 1:nSides: info on small Mortar sides:
                                               ! -1: is neighbor small mortar , 0: not a small mortar, 1: small mortar on big side
LOGICAL             :: HDGDisplayConvergence  !< Display divergence criteria: Iterations, Runtime and Residual
REAL                :: RunTime                !< CG Solver runtime
REAL                :: RunTimePerIteration    !< CG Solver runtime per iteration
REAL                :: HDGNorm                !< Norm
INTEGER             :: iteration              !< number of iterations to achieve the norm 

! --- Boltzmann relation (BR) electron fluid
LOGICAL               :: UseBRElectronFluid            ! Indicates usage of BR electron fluid model
INTEGER               :: BRNbrOfRegions                  ! Nbr of regions to be mapped to Elems
INTEGER, ALLOCATABLE  :: ElemToBRRegion(:)             ! ElemToBRRegion(1:nElems)
REAL, ALLOCATABLE     :: RegionBounds(:,:)             ! RegionBounds ((xmin,xmax,ymin,...)|1:BRNbrOfRegions)
REAL, ALLOCATABLE     :: RegionElectronRef(:,:)        ! RegionElectronRef((rho0,phi0,Te[eV])|1:BRNbrOfRegions)
REAL                  :: BRTimeStepMultiplier          ! Factor that is multiplied with the ManualTimeStep when using BR model
REAL                  :: BRTimeStepBackup              ! Original time step
#if defined(PARTICLES)
! --- Switching between BR and fully kinetic HDG
LOGICAL               :: BRConvertElectronsToFluid     ! User variable for removing all electrons and using BR instead
REAL                  :: BRConvertElectronsToFluidTime ! Time when kinetic electrons should be converted to BR fluid electrons
LOGICAL               :: BRConvertFluidToElectrons     ! User variable for creating particles from BR electron fluid (uses
REAL                  :: BRConvertFluidToElectronsTime ! Time when BR fluid electrons should be converted to kinetic electrons
INTEGER               :: BRConvertMode                 ! Mode used for switching BR->kin->BR OR kin->BR->kin
!                                                      ! and ElectronDensityCell ElectronTemperatureCell from .h5 state file)
LOGICAL               :: BRConvertModelRepeatedly      ! Repeat the switch between BR and kinetic multiple times
LOGICAL               :: BRElectronsRemoved            ! True if electrons were removed during restart (only BR electrons)
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
SUBROUTINE Mask_MPIsides(firstdim,v)
!===================================================================================================================================
! communicate contribution from MPI slave sides to MPI master sides  and set slaves them to zero afterwards.
!===================================================================================================================================
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
