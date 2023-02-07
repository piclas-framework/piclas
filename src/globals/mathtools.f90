!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "piclas.h"

!==================================================================================================================================
!> Routines providing general math functions
!==================================================================================================================================
MODULE MOD_Mathtools
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE INVERSE
   MODULE PROCEDURE INVERSE
END INTERFACE

INTERFACE INVERSE_LU
  MODULE PROCEDURE INVERSE_LU
END INTERFACE
PUBLIC::INVERSE
PUBLIC::INVERSE_LU
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes matrix inverse using LAPACK
!> Input matrix should be a square matrix
!==================================================================================================================================
FUNCTION INVERSE(A) RESULT(AINV)
! MODULES
USE MOD_Globals, ONLY: abort
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: A(:,:)                      !< input matrix
REAL             :: AINV(SIZE(A,1),SIZE(A,2))   !< result: inverse of A
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI
! LOCAL VARIABLES
REAL    :: sdet
REAL    :: work(SIZE(A,1))  ! work array for LAPACK
INTEGER :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER :: n,info
!==================================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

IF (n.EQ.2) THEN
  sdet=1./(A(1,1) * A(2,2) - A(1,2)*A(2,1))
  AINV(1,1) = (  A(2,2) ) * sdet
  AINV(1,2) = (- A(1,2) ) * sdet
  AINV(2,1) = (- A(2,1) ) * sdet
  AINV(2,2) = (  A(1,1) ) * sdet
ELSE IF (n.EQ.3) THEN
  sdet = 1./ (  ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * A(3,3) &
              + ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * A(3,1) &
              + ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * A(3,2))
  AINV(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) * sdet
  AINV(1,2) = ( A(1,3) * A(3,2) - A(1,2) * A(3,3) ) * sdet
  AINV(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * sdet
  AINV(2,1) = ( A(2,3) * A(3,1) - A(2,1) * A(3,3) ) * sdet
  AINV(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1) ) * sdet
  AINV(2,3) = ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * sdet
  AINV(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) * sdet
  AINV(3,2) = ( A(1,2) * A(3,1) - A(1,1) * A(3,2) ) * sdet
  AINV(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * sdet
ELSE
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  CALL DGETRF(n, n, Ainv, n, ipiv, info)

  IF(info.NE.0)THEN
    CALL abort(__STAMP__,'INVERSE(A): Matrix is numerically singular! INFO = ',IntInfoOpt=INFO)
  END IF

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

  IF(info.NE.0)THEN
    CALL abort(__STAMP__,'INVERSE(A): Matrix inversion failed! INFO = ',IntInfoOpt=INFO)
  END IF
END IF
END FUNCTION INVERSE


!==================================================================================================================================
!> Computes matrix inverse using Doolittle LU factorization for Ax=b
!> Input matrix should be a square matrix
!==================================================================================================================================
FUNCTION INVERSE_LU(A) RESULT(AINV)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: A(:,:)                      !< Input matrix
REAL             :: AINV(SIZE(A,1),SIZE(A,2))   !< Result: inverse of A
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, DIMENSION(SIZE(A,1),SIZE(A,2)) :: U,L,A2  !< Upper and Lower part of A and copy
REAL, DIMENSION(SIZE(A,1))           :: b,d,x   !< RHS and aux vectors
INTEGER                              :: i,j,k   !< Loop indices
INTEGER                              :: n       !< Size of first dimension of A, assume that A is a square matirx
REAL                                 :: c       !< Auxiliary coefficient
!==================================================================================================================================
n = size(A,1)

IF (n.LE.3) THEN
  ! Call other routine that calculates the exact inverse
  AINV = INVERSE(A)
ELSE
  ! 0.) Store A in A2 to prevent it from being overwritten
  ! and init lower and upper matrices L and U and RHS b
  A2 = A
  L  = 0.
  U  = 0.
  b  = 0.

  ! 1.) Forward elimination
  DO k=1, n-1
    DO i=k+1,n
      c      = A2(i,k)/A2(k,k)
      L(i,k) = c
      DO j=k+1,n
        A2(i,j) = A2(i,j)-c*A2(k,j)
      END DO
    END DO
  END DO

  ! 2.) Set lower L and upper U matrices
  ! L matrix: A2 matrix of the elimination coefficient and diagonal elements are unity
  DO i=1,n
    L(i,i) = 1.
  END DO

  ! U matrix: Upper part of A2
  DO j=1,n
    DO i=1,j
      U(i,j) = A2(i,j)
    END DO
  END DO

  ! 3.) Columns of the inverse matrix AINV
  DO k=1,n
    b(k)=1.
    d(1) = b(1)

    ! 4.) Forward substitution: Solve L*d = b
    DO i=2,n
      d(i)=b(i)
      DO j=1,i-1
        d(i) = d(i) - L(i,j)*d(j)
      END DO
    END DO

    ! 5.) Backward substitution: Solve U*x = d
    x(n)=d(n)/U(n,n)
    DO i = n-1,1,-1
      x(i) = d(i)
      DO j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      END DO
      x(i) = x(i)/u(i,i)
    END DO

    ! 6.) Copy x into AINV
    DO i=1,n
      AINV(i,k) = x(i)
    END DO
    b(k)=0.
  END DO
END IF
END FUNCTION INVERSE_LU


END MODULE MOD_Mathtools
