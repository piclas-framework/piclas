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
REAL                 :: eps_LinearSolver,eps2_LinearSolver,epsTilde_LinearSolver
REAL,ALLOCATABLE     :: ImplicitSource(:,:,:,:,:)
REAL,ALLOCATABLE     :: LinSolverRHS  (:,:,:,:,:)
REAL,ALLOCATABLE     :: FieldSource(:,:,:,:,:,:)
REAL,ALLOCATABLE     :: Upast(:,:,:,:,:,:)
REAL,ALLOCATABLE     :: Mass(:,:,:,:,:)
INTEGER              :: LinSolver
INTEGER              :: nKDim,nRestarts
INTEGER              :: nDofGlobal, nDofGlobalMPI
LOGICAL              :: LinearSolverInitIsDone=.FALSE.
INTEGER              :: nGP2D,nGP3D                ! (N+1)**2, (N+1)**3
INTEGER              :: nDOFside,nDOFelem          ! nVar*nGP2D, nVar*nGP3D
INTEGER              :: nDOFLine
INTEGER              :: maxIter_LinearSolver
INTEGER              :: totalIterLinearSolver,nInnerIter
INTEGER              :: ldim
#if defined(PARTICLES)
#if defined(IMPA) || (PP_TimeDiscMethod==110)
INTEGER              :: totalPartIterLinearSolver,nPartInnerIter
INTEGER              :: nPartNewton
INTEGER              :: nPartNewtonIter
INTEGER              :: FreezePartInNewton
REAL                 :: Eps2PartNewton
LOGICAL              :: EisenstatWalker
REAL                 :: PartgammaEW
REAL                 :: rEps0,srEps0
REAL,ALLOCATABLE     :: PartXK(:,:)
REAL,ALLOCATABLE     :: R_PartXK(:,:)
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
!===================================================================================================================================
END MODULE MOD_LinearSolver_Vars
