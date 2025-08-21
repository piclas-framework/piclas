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
#if USE_PETSC && USE_HDG
#include "petsc/finclude/petsc.h"
#endif /*USE_PETSC && USE_HDG*/

!===================================================================================================================================
!> Module for the HDG method
!===================================================================================================================================
MODULE MOD_HDG_PETSc
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_PETSC && USE_HDG
PUBLIC :: PETScSetSolver
#endif /*USE_PETSC && USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_PETSC && USE_HDG
SUBROUTINE PETScSetSolver()
!===================================================================================================================================
!> Set the solver and/or preconditioner combination in PETSc
!> Iterative solvers
!>    1: CG + Block Jacobi
!>    2: Pipelined CG + Block Jacobi
!>    3: GMRES + BoomerAMG (with hypre) or Block Jacobi (built-in)
!> Direct solvers
!>    10: CHOLESKY (requires the MUMPS package to support the matrix type)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars
USE PETSc
USE MOD_HDG_Vars_PETSc
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
PetscErrorCode      :: ierr
PC                  :: pc
Mat                 :: F
CHARACTER(LEN=100)  :: ksp_type, pc_type, mat_type
!===================================================================================================================================
PetscCallA(KSPCreate(PETSC_COMM_WORLD,PETScSolver,ierr))
PetscCallA(KSPSetOperators(PETScSolver,PETScSystemMatrix,PETScSystemMatrix,ierr))

PetscCallA(KSPGetPC(PETScSolver,pc,ierr))
! Set the tolerances defaults: rtol=1e-5, atol=1e-50, dtol=1e5, maxits=1e4
! ASSOCIATE( rtol => PETSC_DEFAULT_REAL )
ASSOCIATE( rtol => 1e-16, atol => epsCG )
SELECT CASE(PrecondType)
CASE(0)
  ! ====== Iterative solver: Conjugate Gradient
  PetscCallA(KSPSetType(PETScSolver,KSPCG,ierr))
  PetscCallA(KSPSetInitialGuessNonzero(PETScSolver,PETSC_TRUE, ierr))
  PetscCallA(KSPSetNormType(PETScSolver, KSP_NORM_UNPRECONDITIONED, ierr))
  PetscCallA(KSPSetTolerances(PETScSolver,rtol,atol,PETSC_DEFAULT_REAL,MaxIterCG,ierr))
  ! ===  Preconditioner: None
  PetscCallA(PCSetType(pc,PCNONE,ierr))
CASE(1)
  ! ====== Iterative solver: Conjugate Gradient
  PetscCallA(KSPSetType(PETScSolver,KSPCG,ierr))
  PetscCallA(KSPSetInitialGuessNonzero(PETScSolver,PETSC_TRUE, ierr))
  PetscCallA(KSPSetNormType(PETScSolver, KSP_NORM_UNPRECONDITIONED, ierr))
  PetscCallA(KSPSetTolerances(PETScSolver,rtol,atol,PETSC_DEFAULT_REAL,MaxIterCG,ierr))
  ! ===  Preconditioner: Block Jacobi
  PetscCallA(PCSetType(pc,PCBJACOBI,ierr))
CASE(2)
  ! ====== Iterative solver: Pipelined Conjugate Gradient (only a single non-blocking communication instead of 2 blocking compared to KSPCG)
  PetscCallA(KSPSetType(PETScSolver,KSPPIPECG,ierr))
  PetscCallA(KSPSetInitialGuessNonzero(PETScSolver,PETSC_TRUE, ierr))
  PetscCallA(KSPSetNormType(PETScSolver, KSP_NORM_UNPRECONDITIONED, ierr))
  ! Tolerances defaults: rtol=1e-5, atol=1e-50, dtol=1e5, maxits=1e4
  PetscCallA(KSPSetTolerances(PETScSolver,rtol,atol,PETSC_DEFAULT_REAL,MaxIterCG,ierr))
  ! ===  Preconditioner: Block Jacobi
  PetscCallA(PCSetType(pc,PCBJACOBI,ierr))
CASE(3)
  ! ====== Iterative solver: Flexible Generalized Minimal Residual method
  PetscCallA(KSPSetType(PETScSolver,KSPFGMRES, ierr))
  ! Number of iterations at which the solver restarts [default = 30]: "A larger restart parameter generally leads to faster convergence
  ! of GMRES but the memory usage is higher than with a smaller restart parameter, as is the average time to perform each iteration.
  ! For more ill-conditioned problems a larger restart value may be necessary." https://petsc.org/release/manualpages/KSP/KSPGMRESSetRestart/
  PetscCallA(KSPGMRESSetRestart(PETScSolver, 100, ierr))
  PetscCallA(KSPSetInitialGuessNonzero(PETScSolver,PETSC_TRUE, ierr))
  PetscCallA(KSPSetNormType(PETScSolver, KSP_NORM_UNPRECONDITIONED, ierr))
  PetscCallA(KSPSetTolerances(PETScSolver,rtol,atol,PETSC_DEFAULT_REAL,MaxIterCG,ierr))
#ifdef PETSC_HAVE_HYPRE
  ! ===  Preconditioner: BoomerAMG
  PetscCallA(PCSetType(pc, PCHYPRE, ierr))
  PetscCallA(PCHYPRESetType(pc, "boomeramg", ierr))
  ! BoomerAMG options
  ! Strong threshold for coarsening: greater value means more coarsening; default = 0.25, which is only sufficient for 2D
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_strong_threshold", "0.7", ierr))
  ! Coarsening strategy: HMIS coarsening
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_coarsen_type", "HMIS", ierr))
  ! Maximum number of levels (default: 25)
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_max_levels", "25", ierr))
  ! Number of coarsening levels for "aggressive coarsening"
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_agg_nl", "4", ierr))
  ! Number of pathways within a coarsening level: 1 is most agressive value; balance between the number of levels and paths
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_agg_num_paths", "5", ierr))
  ! Interpolation type
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_interp_type", "ext+i", ierr))
  ! Coarsen during the interpolation
  PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pc_hypre_boomeramg_truncfactor", "0.3", ierr))
  PetscCallA(PCSetFromOptions(pc,ierr))
#else /*NOT PETSC_HAVE_HYPRE*/
  ! ===  Preconditioner: Block Jacobi
  PetscCallA(PCSetType(pc,PCBJACOBI,ierr))
#endif /*PETSC_HAVE_HYPRE*/
#ifdef PETSC_HAVE_MUMPS
CASE(10)
  ! ====== Direct solver: Cholesky
  PetscCallA(KSPSetType(PETScSolver,KSPPREONLY,ierr))
  PetscCallA(PCSetType(pc,PCCHOLESKY,ierr))
  ! PETSc will most likely use MUMPS anyway
  PetscCallA(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr))
  PetscCallA(PCFactorSetUpMatSolverType(pc,ierr))
  ! We need to get the internal matrix to set its options
  PetscCallA(PCFactorGetMatrix(pc,F,ierr))
#if USE_DEBUG
  ! Increase MUMPS diagnostics level: Errors, warnings, and main statistics printed.
  PetscCallA(MatMumpsSetIcntl(F, 4, 2, ierr))
#endif /*USE_DEBUG*/
  ! === Compression
  ! Enable BLR compression with automatic settings: showed better performance for initial factorization and better memory footprint
  PetscCallA(MatMumpsSetIcntl(F, 35, 1, ierr))
  ! ! Enable BLR compression of the contribution blocks, reducing the memory consumption at the cost of some additional operations
  ! ! during factorization
  ! PetscCallA(MatMumpsSetIcntl(F, 37, 1, ierr))
  ! === Parallel ordering: select one of the following or let PETSc decide (recommended)
  ! PetscCallA(MatMumpsSetIcntl(F, 28, 2, ierr))
  ! ! Use PT-SCOTCH for ordering
  ! PetscCallA(MatMumpsSetIcntl(F, 29, 1, ierr))
  ! ! Use ParMetis for parallel ordering
  ! PetscCallA(MatMumpsSetIcntl(F, 29, 2, ierr))

  ! === Memory handling
  ! Workspace allocation: Allow 2x estimated memory (default is at 35%)
  PetscCallA(MatMumpsSetIcntl(F,14,100,ierr))
  ! ! Limit to 2GB per process, or default (=0): each processor will allocate workspace based on the estimates computed during the analysis
  ! PetscCallA(MatMumpsSetIcntl(F,23,2000,ierr))
#endif /*PETSC_HAVE_MUMPS*/
CASE DEFAULT
  SWRITE(*,*) 'PrecondType:', PrecondType
  CALL CollectiveStop(__STAMP__,'ERROR in PETScSetSolver: Unknown option! Direct solver (10) is only available with MUMPS.')
END SELECT
END ASSOCIATE

! Get solver and preconditioner types
PetscCallA(KSPGetType(PETScSolver, ksp_type, ierr))
PetscCallA(PCGetType(pc, pc_type, ierr))

! Reuse preconditioner (might be unneccessary since the system matrix remains the same during the simulation)
PetscCallA(KSPSetReusePreconditioner(PETScSolver, PETSC_TRUE, ierr))

! If using direct solver, print factorization type
IF (TRIM(ksp_type) .EQ. 'preonly') THEN
  ! Print factorization details when using Cholesky/LU
  IF ((TRIM(pc_type) .EQ. 'cholesky') .OR. (TRIM(pc_type) .EQ. 'lu')) then
    PetscCallA(PCFactorGetMatrix(pc, F, ierr))
    PetscCallA(MatGetType(F, mat_type, ierr))
    LBWRITE(UNIT_stdOut,'(A)') ' | Direct solver: '//TRIM(pc_type)//', using factorization type: '//TRIM(mat_type)
  END IF
ELSE
  LBWRITE(UNIT_stdOut,'(A)') ' | Iterative solver: '//TRIM(ksp_type)//' with '//TRIM(pc_type)//' preconditioning'
END IF

END SUBROUTINE PETScSetSolver
#endif /*USE_PETSC && USE_HDG*/

END MODULE MOD_HDG_PETSc
