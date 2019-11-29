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

MODULE MOD_Utils
!===================================================================================================================================
! Contains utils required by xy- modules
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE BubbleSortID
  MODULE PROCEDURE BubbleSortID
END INTERFACE

INTERFACE InsertionSort
  MODULE PROCEDURE InsertionSort
END INTERFACE

INTERFACE QSort1Doubleint1Pint
  MODULE PROCEDURE QSort1DoubleInt1Pint
END INTERFACE

INTERFACE RootsOfBesselFunctions
  MODULE PROCEDURE RootsOfBesselFunctions
END INTERFACE

INTERFACE QuadraticSolver
  MODULE PROCEDURE QuadraticSolver
END INTERFACE

PUBLIC:: BubbleSortID
PUBLIC:: InsertionSort
PUBLIC:: QSort1Doubleint1Pint
PUBLIC:: RootsOfBesselFunctions
PUBLIC:: QuadraticSolver
!===================================================================================================================================

CONTAINS

PURE SUBROUTINE InsertionSort(a,id,len)
!===================================================================================================================================
! Insertion sort
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: len
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                :: a(1:len)
INTEGER,INTENT(INOUT),OPTIONAL    :: id(1:len)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: tmpR
INTEGER                           :: i, j, tmpI
!===================================================================================================================================

IF(PRESENT(ID))THEN
  DO i=2,len
    j=i-1
    tmpR=a(i)
    tmpI=ID(i)
    DO WHILE (j.GE.1) !(j.GE.1 .AND. a(j).GT.tmpR)
      IF (a(j).LE.tmpR) EXIT
      a (j+1) = a(j)
      ID(j+1) = ID(j)
      j=j-1
    END DO
    a (j+1) =tmpR
    ID(j+1) =tmpI
  END DO ! i
ELSE
  DO i=2,len
    j=i-1
    tmpR=a(i)
    DO WHILE (j.GE.1) !(j.GE.1 .AND. a(j).GT.tmpR)
      IF (a(j).LE.tmpR) EXIT
      a (j+1) = a(j)
      j=j-1
    END DO
    a (j+1) =tmpR
  END DO ! i
END IF

END SUBROUTINE InsertionSort


SUBROUTINE BubbleSortID(a,id,len)
!===================================================================================================================================
! bubble sort, taken from rosetta-wiki and modified for own use
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: len
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                :: a(len)
INTEGER,INTENT(INOUT),OPTIONAL    :: id(len)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: temp
INTEGER                           :: iloop,jloop, temp2
LOGICAL                           :: swapped = .TRUE.
!===================================================================================================================================

IF(PRESENT(id))THEN
  DO jloop=len-1,1,-1
    swapped = .FALSE.
    DO iloop=1,jloop
      IF (a(iloop).GT.a(iloop+1))THEN
        ! switch entries
        temp=a(iloop)
        a(iloop) = a(iloop+1)
        a(iloop+1) = temp
        ! switch ids
        temp2=id(iloop)
        id(iloop) = id(iloop+1)
        id(iloop+1) = temp2
        swapped = .TRUE.
      END IF
    END DO ! iloop
    IF (.NOT. swapped) EXIT
  END DO ! jloop
ELSE
  DO jloop=len-1,1,-1
    swapped = .FALSE.
    DO iloop=1,jloop
      IF (a(iloop).GT.a(iloop+1))THEN
        ! switch entries
        temp=a(iloop)
        a(iloop) = a(iloop+1)
        a(iloop+1) = temp
        swapped = .TRUE.
      END IF
    END DO ! iloop
    IF (.NOT. swapped) EXIT
  END DO ! jloop
END IF
END SUBROUTINE BubbleSortID


RECURSIVE SUBROUTINE Qsort1DoubleInt1PInt(A,P)
!===================================================================================================================================
! QSort1int: (c) by Mark Haas
!  Uses the Quicksort algorithm to sort an INTEGER table with one relevant columns, along with a passive onedimensional array.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER(KIND=8),INTENT(INOUT)    :: A(:) ! active array, will be sorted
   INTEGER,INTENT(INOUT)            :: P(:) ! passive array, is sorted like A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: marker  ! ?
!===================================================================================================================================

   IF(SIZE(A) .GT. 1) THEN
      CALL Partition1DoubleInt1Pint(A,P,marker)
      CALL Qsort1DoubleInt1Pint(A(:marker-1),P(:marker-1))
      CALL Qsort1DoubleInt1Pint(A(marker:),P(marker:))
   END IF
 RETURN
END SUBROUTINE Qsort1DoubleInt1PInt


SUBROUTINE Partition1DoubleInt1PInt(A,P,marker)
!===================================================================================================================================
!  Sorting routine used by QSort1int above. This routine is PRIVATE
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER(KIND=8),INTENT(INOUT)    :: A(:) ! active array, will be sorted
   INTEGER,INTENT(INOUT)            :: P(:) ! passive array, is sorted like A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   INTEGER,INTENT(OUT)              :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: i,j  ! ?
   INTEGER(KIND=8)                  :: temp,x  ! ?
   INTEGER                          :: ptemp  ! ?
!===================================================================================================================================
   x = A(1)
   i = 0
   j  = SIZE(A)+1

   DO
      j = j-1
      DO
         IF (A(j) .LE. x)THEN
            EXIT
         END IF
         j = j-1
      END DO
      i = i+1
      DO
         IF(A(i) .GE. x)THEN
            EXIT
         END IF
         i = i+1
      END DO
      IF (i .LT. j) THEN
         ! exchange A(i) and A(j), P(i) and P(j)
         temp = A(i)
         A(i)  = A(j)
         A(j)  = temp
         ptemp = P(i)
         P(i)  = P(j)
         P(j)  = ptemp
      ELSE IF (i .EQ. j) THEN
         marker = i+1
         RETURN
      ELSE
         marker = i
         RETURN
      END IF
   END DO
 RETURN
END SUBROUTINE Partition1DoubleInt1PInt


SUBROUTINE RootsOfBesselFunctions( targetN, targetM, p, zo )
!===================================================================================================================================
! Computes the first targetM roots of the targetN-th Bessel function or derivative of Bessel function
! Required for TE/TE modes
! p-0 TE mode
! p-1 TM mode
!
! routine taken from
! http://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
! at
! 2017-10-13
! JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
!
!  Discussion:
!
!    This procedure computes the zeros of Bessel functions Jn(x) and
!    Jn'(x), and arrange them in the order of their magnitudes.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    13 October 2017
!    01 August  2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin, Philip Ortwein
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    INPUT, INTEGER(kind=4) targetN, targetM

!    INPUT, INTEGER(kind=4) p 0-TE mode, 1 TM mode

!    of Jn(x) correspond to TM modes and those of Jn'(x) correspond to TE modes.
!
!    Output, real ( kind = 8 ) ZO(*), the zeros of Jn(x) and Jn'(x).
!
!===================================================================================================================================
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=4),INTENT(IN) ::  p ! select TE-0 or TM-1 mode
INTEGER(KIND=4),INTENT(IN) ::  targetN
INTEGER(KIND=4),INTENT(IN) ::  targetM
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(KIND=8),INTENT(OUT)   ::  zo(1:targetM)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND= 4)  :: j
INTEGER(KIND= 4)  :: k
INTEGER(KIND= 4)  :: l
INTEGER(KIND= 4)  :: l0
INTEGER(KIND= 4)  :: l1
INTEGER(KIND= 4)  :: l2
REAL ( KIND = 8)  :: x
REAL ( KIND = 8)  :: x0
REAL ( KIND = 8)  :: x1
REAL ( KIND = 8)  :: x2
REAL ( KIND = 8)  :: zoc(1:targetM)
REAL ( KIND = 8)  :: dbessel
REAL ( KIND = 8)  :: d2bessel
REAL ( KIND = 8)  :: bessel
!===================================================================================================================================

  l0 = 0
  x1 = 0.407658D+00 + 0.4795504D+00 &
    * ( real ( targetN, kind = 8 ) ) ** 0.5D+00 + 0.983618D+00 * ( targetN )
  x2 = 1.99535D+00 + 0.8333883 * ( real ( targetN, kind = 8 ) ) ** 0.5D+00 &
    + 0.984584D+00 * ( targetN )
  l1 = 0

  do j = 1, targetM

    if(p.EQ.0)THEN ! TE modes
      if ( targetN+1 == 1 .and. j == 1 ) then

        l1 = l1 + 1
        zoc(l1) = x

        if ( targetN+1 <= 15 ) then
          x1 = x + 3.057D+00 + 0.0122D+00 * ( targetN ) &
            + ( 1.555D+00 + 0.41575D+00 * ( targetN ) ) / ( j + 1 ) ** 2
        else
          x1 = x + 2.918D+00 + 0.01924D+00 * ( targetN ) &
            + ( 6.26D+00 + 0.13205D+00 * ( targetN ) ) / ( j + 1 ) ** 2
        end if

      else

        x = x1

        do
          ! newton step

          ! call bjndd ( i, x, bj, dj, fj )
          ! WE use the inbuilt functions
          !print*,'i,x',i,x,bj(i),dj(i),fj(i)
          ! BESSEL_JN(i,x) -> BJ(i)
          ! dbessel -> dj
          ! d2bessel -> fj
          ! again, rememger  the shift for i
          !bessel=BESSEL_JN(targetN,x)
          dbessel=0.5*(BESSEL_JN(targetN-1,x)-BESSEL_JN(targetN+1,x))
          d2bessel=0.25*(BESSEL_JN(targetN-2,x)-2*Bessel_JN(targetN,x)+BESSEL_JN(targetN+2,x)) ! again shifted

          x0 = x
          x = x - dbessel / d2bessel

          if ( abs ( x - x0 ) <= 1.0D-10 ) then
            l1 = l1 + 1
            zoc(l1) = x

            if ( targetN+1 <= 15 ) then
              x1 = x + 3.057D+00 + 0.0122D+00 * ( targetN ) &
                + ( 1.555D+00 + 0.41575D+00 * ( targetN ) ) / ( j + 1 ) ** 2
            else
              x1 = x + 2.918D+00 + 0.01924D+00 * ( targetN ) &
                + ( 6.26D+00 + 0.13205D+00 * ( targetN ) ) / ( j + 1 ) ** 2
            end if
            exit
          end if

        end do
      end if

    ELSE
    !  TM modes

     x = x2

     do

       !call bjndd ( targetN+1, x, bj, dj, fj )
       ! maybe to be fixed
       bessel=BESSEL_JN(targetN,x)
       dbessel=0.5*(BESSEL_JN(targetN-1,x)-BESSEL_JN(targetN+1,x))
       !d2bessel=0.25*(BESSEL_JN(targetN-2,x)-2*Bessel_JN(targetN,x)+BESSEL_JN(targetN+2,x)) ! again shifted

       x0 = x
       !x = x - bj(targetN+1) / dj(targetN+1)
       x = x - bessel / dbessel

       if ( abs ( x - x0 ) <= 1.0D-10 ) then
         exit
       end if

     end do

     l1 = l1 + 1
     zoc(l1) = x
     if ( targetN+1 <= 15 ) then
       x2 = x + 3.11D+00 + 0.0138D+00 * ( targetN ) &
         + ( 0.04832D+00 + 0.2804D+00 * ( targetN ) ) / ( j + 1 ) ** 2
     else
       x2 = x + 3.001D+00 + 0.0105D+00 * ( targetN ) &
         + ( 11.52D+00 + 0.48525D+00 * ( targetN ) ) / ( j + 3 ) ** 2
     end if

    END IF

  end do

  l = l0 + l1
  l2 = l

  do

    if ( l0 == 0 ) then
      do k = 1, l
        zo(k) = zoc(k)
      end do
      l1 = 0
    else if ( l0 /= 0 ) then
      if ( zoc(l1) .le. zo(l0) ) then
        zo(l0+l1) = zo(l0)
        l0 = l0 - 1
      else
        zo(l0+l1) = zoc(l1)
        l1 = l1 - 1
      end if
    end if

    if ( l1 == 0 ) then
      exit
    end if

  end do

  l0 = l2
END SUBROUTINE RootsOfBesselFunctions


!================================================================================================================================
!> Stable algorithm to compute the number of (nRoot) non-complex solutions R1 and R2 for the quadratic equation A*x^2 + B*x + C 
!================================================================================================================================
SUBROUTINE QuadraticSolver(A,B,C,nRoot,r1,r2)
#ifdef CODE_ANALYZE
USE MOD_Globals,            ONLY:UNIT_stdOut,MyRank
#endif /*CODE_ANALYZE*/
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: A,B,C
!--------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(OUT)     :: nRoot
REAL,INTENT(OUT)        :: R1,R2
!--------------------------------------------------------------------------------------------------------------------------------
! local variables
REAL                    :: radicant
!================================================================================================================================

! Use P-Q-formula and calculate first solution R1
! Use Theorem of Vieta (R1*R2=C/A) to calculate R2
! cf: wikipedia 2017.06.13 https://de.wikipedia.org/wiki/Quadratische_Gleichung
IF (A.NE.0. .AND. B.EQ.0. .AND. C.EQ.0.) THEN
  nRoot=1
  R1=0.
  R2=0.
ELSE IF(A.NE.0.)THEN
  radicant = (0.5*B/A)**2 - (C/A)
  IF (radicant.LT.0.) THEN
    nRoot=0
    R1=0.
    R2=0.
  ELSE
    nRoot=2
    R1=-0.5*(B/A)-SIGN(1.,B/A)*SQRT(radicant)
    R2=(C/A)/R1
  END IF
ELSE
  IF(B.NE.0.)THEN
    nRoot=1
    R1=-C/B
    R2=0.
  ELSE
    nRoot=0
    R1=0.
    R2=0.
  END IF
END IF

#ifdef CODE_ANALYZE
IF(nRoot.GT.0)THEN
  IF((ABS(R1).LE.1.).AND.(ABS(A*R1**2+B*R1+C).GT.1e-10))THEN
    IPWRITE(UNIT_stdOut,'(I0,A,G0,A)')    ' WARNING!!!: RHS of R1 is ',A*R1**2+B*R1+C &
        ,' (.GT.1e-10) in quadratic solver.'
  END IF
END IF
IF(nRoot.GT.1)THEN
  IF((ABS(R2).LE.1.).AND.(ABS(A*R2**2+B*R2+C).GT.1e-10))THEN
    IPWRITE(UNIT_stdOut,'(I0,A,G0,A)')    ' WARNING!!!: RHS of R2 is ',A*R2**2+B*R2+C &
        ,' (.GT.1e-10) in quadratic solver.'
  END IF
END IF
#endif /*CODE_ANALYZE*/

END SUBROUTINE QuadraticSolver


END MODULE MOD_Utils
