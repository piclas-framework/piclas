#include "boltzplatz.h"

MODULE MOD_LinearSolver_Vars
!===================================================================================================================================
! Contains global variables used by the Timedisc modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                 :: eps_LinearSolver,eps2_LinearSolver,epsTilde_LinearSolver    ! abort tolerance for DG linear solver
REAL,ALLOCATABLE     :: ImplicitSource(:,:,:,:,:)                                   ! temp. storage of source terms
REAL,ALLOCATABLE     :: LinSolverRHS  (:,:,:,:,:)                                   ! RHS for linear solver
REAL,ALLOCATABLE     :: FieldStage(:,:,:,:,:,:)                                     ! FieldStage, don't no of used
REAL,ALLOCATABLE     :: Upredict(:,:,:,:,:,:)                                       ! Upredictor, don't no of used
REAL,ALLOCATABLE     :: Upast(:,:,:,:,:,:)                                          ! history of upast, required for predictor
REAL,ALLOCATABLE     :: tpast(:)                                                    ! history of tpast, required for predictor
REAL,ALLOCATABLE     :: Mass(:,:,:,:,:)                                             ! mass matrix
INTEGER              :: LinSolver                                                   ! selection of linear solver, CGS,BiCGStab,...
INTEGER              :: nKDim,nRestarts                                             ! Number of Subspaces GMRES  and Restarts
INTEGER              :: nDofGlobal, nDofGlobalMPI                                   ! number of DOFs or global DOFs
LOGICAL              :: LinearSolverInitIsDone=.FALSE.                              ! init routine
INTEGER              :: nGP2D,nGP3D                ! (N+1)**2, (N+1)**3             ! number of GP points on face, volume
INTEGER              :: nDOFside,nDOFelem          ! nVar*nGP2D, nVar*nGP3D         ! number of GP*nvar
INTEGER              :: nDOFLine                                                    ! number of 1D DOFs
INTEGER              :: maxIter_LinearSolver                                        ! limit of iter for linear solver
INTEGER              :: totalIterLinearSolver,nInnerIter                            ! global counter of iterations
INTEGER              :: ldim                                                        ! Number of BiCGStab(l) subspaces
#if defined(PARTICLES)
#if defined(IMPA) || (PP_TimeDiscMethod==110)
LOGICAL              :: DoFieldUpdate
INTEGER              :: totalPartIterLinearSolver,nPartInnerIter                    ! Counter for Particle newton
INTEGER              :: nPartNewton                                                 ! some limits or counter
INTEGER              :: nPartNewtonIter                                             ! some limits or counter
INTEGER              :: FreezePartInNewton                                          ! particle is moved after each Newton step
REAL                 :: Eps2PartNewton                                              ! PartNewton abort criterion
LOGICAL              :: EisenstatWalker                                             ! EisenstatWalker for ParticleNewton
REAL                 :: EpsPartLinSolver                                            ! Abort tolerance for linear solver of parts
REAL                 :: PartgammaEW                                                 ! gamma value of PartEisenstatWalker
REAL                 :: rEps0,srEps0                                                ! FD-step-size for PartMV in PartNewton
REAL,ALLOCATABLE     :: PartXK(:,:)                                                 ! ParticlePosition for linearization
REAL,ALLOCATABLE     :: R_PartXK(:,:)                                               ! Part_dt of PartXK
REAL,PARAMETER       :: Part_alpha=0.0001
REAL,PARAMETER       :: Part_sigma(1:2) = (/0.1, 0.5/)
#endif
#endif /*PARTICLES*/
#if (PP_TimeDiscMethod==104)
INTEGER              :: nNewton
INTEGER              :: nNewtonIter
REAL,ALLOCATABLE     :: R_xk(:,:,:,:,:)
REAL,ALLOCATABLE     :: xk(:,:,:,:,:)
REAL                 :: Eps2Newton
LOGICAL              :: EisenstatWalker
REAL                 :: gammaEW
#endif
#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) ||(PP_TimeDiscMethod==122) || (PP_TimeDiscMethod==131)
REAL                 :: PartNewtonRelaxation                                        ! scaling factor for lambda. A value <0
                                                                                    ! disables Armijo rule and uses a fixed value
REAL,ALLOCATABLE     :: ExplicitPartSource(:,:,:,:,:)                               ! temp. storage of source terms 121,122
LOGICAL              :: DoPrintConvInfo =.FALSE.                                    ! flag to print current norm in outer iteration
                                                                                    ! and number of parts in Newton
INTEGER              :: maxFullNewtonIter                                           ! limit of fullnewton iterations
INTEGER              :: TotalFullNewtonIter                                         ! counter for all total full newton iters
REAL                 :: Eps_FullNewton                                              ! abort tolerance for fullnewtoniter
REAL                 :: Eps2_FullNewton                                             ! square of above
INTEGER              :: FullEisenstatWalker                                         ! Switch for outer eisenstat walker
                                                                                    ! 0 - no Eisenstat-Walker
                                                                                    ! 1 - Field Solver
                                                                                    ! 2 - Particle Newton and Field Solver
REAL                 :: PartRelaxationFac                                           ! relaxation factor for particles
REAL                 :: PartRelaxationFac0                                          ! relaxation factor for particles
INTEGER              :: AdaptIterRelaxation0                                        ! iter to adapt relaxation
LOGICAL              :: DoPartRelaxation                                            ! flag for particle relaxation
REAL                 :: FullgammaEW                                                 ! Eisenstat-Walker parameter
INTEGER              :: PartImplicitMethod                                          ! selection for particle implicit method
#ifdef PARTICLES
LOGICAL              :: DoUpdateInStage                                             ! perform updatenextfree position 
                                                                                    ! in each rk stage
LOGICAL              :: DoFullNewton                                                ! use a full Newton instate of iteration 
INTEGER              :: UpdateInIter                                                ! additional update in iteration. required
                                                                                    ! due to overflow of free positions...
                                                                                    ! UNFP each nth iteration
                                                                                    ! scheme
#endif /*PARTICLES*/
#endif
!===================================================================================================================================
END MODULE MOD_LinearSolver_Vars
