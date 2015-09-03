#include "boltzplatz.h"

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

PUBLIC:: BubbleSortID
PUBLIC:: InsertionSort
PUBLIC:: QSort1Doubleint1Pint
!===================================================================================================================================

CONTAINS

SUBROUTINE InsertionSort(a,id,len)
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
INTEGER,INTENT(INOUT)             :: id(1:len)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: tmpR
INTEGER                           :: i, j, tmpI
!===================================================================================================================================

DO i=2,len
  j=i-1  
  tmpR=a(i)
  tmpI=ID(i)
  DO WHILE (j.GE.1 .AND. a(j).GT.tmpR)
    a (j+1) = a(j)
    ID(j+1) = ID(j)
    j=j-1
  END DO
  a (j+1) =tmpR
  ID(j+1) =tmpI
END DO ! i

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
INTEGER,INTENT(INOUT)             :: id(len)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: temp
INTEGER                           :: iloop,jloop, temp2
LOGICAL                           :: swapped = .TRUE.
!===================================================================================================================================

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


END MODULE MOD_Utils
