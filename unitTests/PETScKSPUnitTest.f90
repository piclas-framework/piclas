!==================================================================================================================================
! Copyright (c) 2024 boltzplatz - numerical plasma dynamics GmbH
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
#include "petsc/finclude/petsc.h"
!==================================================================================================================================
!> Unit test 'PETScKSPUnitTest'
!> Test the PETSc installation using the direct solver and the CG solver by solving a simple 3x3 matrix
!==================================================================================================================================

PROGRAM PETScKSPUnitTest
! MODULES
USE MOD_Globals,    ONLY: Abort
USE PETSc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode  :: ierr
Vec             :: x, b, u
Mat             :: A
KSP             :: ksp
PetscInt        :: n, i, j, Ii, Jj
PetscScalar     :: v, one, neg_one
PetscReal       :: norm, tol
!==================================================================================================================================

n = 3  ! Size of the linear system

! Initialize PETSc
PetscCallA(PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-log_view', PETSC_NULL_CHARACTER, ierr));

PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))
PetscCallA(PetscLogDefaultBegin(ierr))

one     = 1.0
neg_one = -1.0

! Create vectors
PetscCallA(VecCreate(PETSC_COMM_WORLD, b, ierr))
PetscCallA(PetscObjectSetName(b, 'RHS', ierr))
PetscCallA(VecSetSizes(b, PETSC_DECIDE, n, ierr))
PetscCallA(VecSetFromOptions(b, ierr))
PetscCallA(VecDuplicate(b, x, ierr))
PetscCallA(VecDuplicate(b, u, ierr))

! Create matrix
PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr))
PetscCallA(MatSetFromOptions(A, ierr))
PetscCallA(MatSetUp(A, ierr))

! Assemble matrix, setting 4 on the diagonal and -1 everywhere else
DO i = 1, n
  Ii = i - 1
  DO j = 1, n
    Jj = j - 1
    IF (i == j) THEN
      v = 4.0
    ELSE
      v = -1.0
    END IF
    PetscCallA(MatSetValues(A, 1, Ii, 1, Jj, v, INSERT_VALUES, ierr))
  END DO
END DO
PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))

! View the matrix
PetscCallA(MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr))

! Set exact solution: u = (1,1,1)
PetscCallA(VecSet(u,one,ierr))
! Compute right-hand-side vector: b = (2,2,2)
PetscCallA(MatMult(A,u,b,ierr))

PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Right-hand side:\n',ierr))
PetscCallA(VecView(b, PETSC_VIEWER_STDOUT_WORLD, ierr))

! --------------------------------------------------------------------------------
! ### Test the CG solver
! Create KSP solver and set options
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, '----------------------------\n',ierr))
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Testing the CG solver:\n',ierr))
PetscCallA(KSPCreate(PETSC_COMM_WORLD, ksp, ierr))

! Select the CG solver by setting KSPCG
PetscCallA(KSPSetOperators(ksp, A, A, ierr))
PetscCallA(KSPSetType(ksp,KSPCG,ierr))
tol = 1.e-7
PetscCallA(KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr))

PetscCallA(KSPSetUp(ksp,ierr))

! Solve the linear system
PetscCallA(KSPSolve(ksp, b, x, ierr))

! View the solution vector
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Solution vector:\n',ierr))
PetscCallA(VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr))

! Check the error: x = neg_one * u + x
PetscCallA(VecAXPY(x, neg_one, u, ierr))
PetscCallA(VecNorm(x, NORM_2, norm, ierr))

! Check if the solution is correct
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Norm of error:\n',ierr))
WRITE(*,*) norm

IF(norm.GT.1e-12) THEN
  CALL Abort(__STAMP__,'ERROR in PETSc unit test: Norm of error is above 1e-12 for the CG solver!')
END IF

PetscCallA(KSPDestroy(ksp, ierr))
! --------------------------------------------------------------------------------
! ### Test the direct solver
! Create KSP solver and set options
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, '----------------------------\n',ierr))
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Testing the direct solver:\n',ierr))
PetscCallA(KSPCreate(PETSC_COMM_WORLD, ksp, ierr))

! Select the direct solve by setting KSPPREONLY
PetscCallA(KSPSetOperators(ksp, A, A, ierr))
PetscCallA(KSPSetType(ksp,KSPPREONLY,ierr))
PetscCallA(KSPSetUp(ksp,ierr))

! Solve the linear system
PetscCallA(KSPSolve(ksp, b, x, ierr))

! View the solution vector
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Solution vector:\n',ierr))
PetscCallA(VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr))

! Check the error: x = neg_one * u + x
PetscCallA(VecAXPY(x, neg_one, u, ierr))
PetscCallA(VecNorm(x, NORM_2, norm, ierr))

! Check if the solution is correct
PetscCallA(PetscPrintf(PETSC_COMM_WORLD, 'Norm of error:\n',ierr))
WRITE(*,*) norm

IF(norm.GT.1e-12) THEN
  CALL Abort(__STAMP__,'ERROR in PETSc unit test: Norm of error is above 1e-12 for the direct solver!')
END IF
! --------------------------------------------------------------------------------
! Final clean up
PetscCallA(KSPDestroy(ksp, ierr))
PetscCallA(MatDestroy(A, ierr))
PetscCallA(VecDestroy(x, ierr))
PetscCallA(VecDestroy(b, ierr))
PetscCallA(VecDestroy(u, ierr))

! Finalize PETSc
PetscCallA(PetscFinalize(ierr))

END PROGRAM PETScKSPUnitTest
