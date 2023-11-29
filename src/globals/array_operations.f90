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

MODULE MOD_Array_Operations
!===================================================================================================================================
! Contains tools for array related operations.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ChangeSizeArray
  MODULE PROCEDURE ChangeSizeArrayLOG1,  ChangeSizeArrayLOG2, &
                   ChangeSizeArrayINT1,  ChangeSizeArrayINT2, &
                   ChangeSizeArrayREAL1, ChangeSizeArrayREAL2, ChangeSizeArrayREAL3
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ChangeSizeArray
!===================================================================================================================================

CONTAINS


SUBROUTINE ChangeSizeArrayLOG1(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Logical 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT) :: Vec(:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
LOGICAL,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE               :: TempVec(:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayLOG1


SUBROUTINE ChangeSizeArrayLOG2(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Logical 2D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT) :: Vec(:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
LOGICAL,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE               :: TempVec(:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(:,1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(:,OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayLOG2


SUBROUTINE ChangeSizeArrayINT1(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Integer 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(INOUT) :: Vec(:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
INTEGER,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE               :: TempVec(:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayINT1


SUBROUTINE ChangeSizeArrayINT2(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Integer 2D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(INOUT) :: Vec(:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
INTEGER,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE               :: TempVec(:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(:,1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(:,OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayINT2


SUBROUTINE ChangeSizeArrayREAL1(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Real 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,ALLOCATABLE,INTENT(INOUT)    :: Vec(:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
REAL,INTENT(IN),OPTIONAL          :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                  :: TempVec(:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayREAL1


SUBROUTINE ChangeSizeArrayREAL2(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Real 2D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,ALLOCATABLE,INTENT(INOUT)    :: Vec(:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
REAL,INTENT(IN),OPTIONAL          :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                  :: TempVec(:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(:,1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(:,OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayREAL2


SUBROUTINE ChangeSizeArrayREAL3(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Real 3D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,ALLOCATABLE,INTENT(INOUT)    :: Vec(:,:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
REAL,INTENT(IN),OPTIONAL          :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                  :: TempVec(:,:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),SIZE(Vec,2),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(:,:,1:OldSize)=Vec
  IF(PRESENT(Default)) TempVec(:,:,OldSize+1:NewSize) = Default
ELSE
  TempVec=Vec(:,:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayREAL3


END MODULE MOD_Array_Operations
