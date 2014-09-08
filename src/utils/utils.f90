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

PUBLIC:: BubbleSortID
!===================================================================================================================================

CONTAINS

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

END MODULE MOD_Utils
