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

MODULE MOD_Vector
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

!===================================================================================================================================
CONTAINS


SUBROUTINE VNullify(nTotal,Vec)
!===================================================================================================================================
! Y=0
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: Vec(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Vec=0.
END SUBROUTINE VNullify


SUBROUTINE VSetConst(nTotal,Vec,Const)
!===================================================================================================================================
! Y=a
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal
REAL,INTENT(IN)       :: Const
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: Vec(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Vec=Const
END SUBROUTINE VSetConst


SUBROUTINE VCopy(nTotal,VecIn,VecOut)
!===================================================================================================================================
! Y=X
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal
REAL,INTENT(IN)       :: VecIn(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: VecOut(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
VecOut=VecIn
END SUBROUTINE VCopy

SUBROUTINE LVCopy(nTotal,VecIn,VecOut)
!===================================================================================================================================
! Y=X
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)    :: nTotal
REAL,INTENT(IN)               :: VecIn(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(KIND=8),INTENT(OUT)      :: VecOut(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
VecOut=VecIn
END SUBROUTINE LVCopy

SUBROUTINE VAX(nTotal,VecIn,ConstIn)
!===================================================================================================================================
! Y=AY+BX
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: nTotal
REAL,INTENT(IN)          :: ConstIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)      :: Vecin(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!===================================================================================================================================
DO i=1,nTotal
  VecIn(i)=VecIn(i)*ConstIn
END DO
END SUBROUTINE VAX

SUBROUTINE VAXPBY(nTotal,VecIn,VecOut,ConstIn,ConstOut)
!===================================================================================================================================
! Y=AY+BX
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal
REAL,INTENT(IN)       :: VecIn(nTotal)
REAL,INTENT(IN),OPTIONAL :: ConstIn
REAL,INTENT(IN),OPTIONAL :: ConstOut
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: VecOut(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!===================================================================================================================================
IF(    PRESENT(ConstIn).AND.(.NOT.PRESENT(ConstOut)))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)+VecIn(i)*ConstIn
  END DO
ELSEIF(PRESENT(ConstOut).AND.(.NOT.PRESENT(ConstIn)))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)*ConstOut+VecIn(i)
  END DO
ELSEIF(PRESENT(ConstIn).AND.PRESENT(ConstOut))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)*ConstOut+VecIn(i)*ConstIn
  END DO
ELSE
  DO i=1,nTotal
    VecOut=VecOut(i)+VecIn(i)
  END DO
END IF
END SUBROUTINE VAXPBY

SUBROUTINE LVAXPBY(nTotal,VecIn,VecOut,ConstIn,ConstOut)
!===================================================================================================================================
! Y=AY+BX
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)    :: nTotal
REAL,INTENT(IN)       :: VecIn(nTotal)
REAL,INTENT(IN),OPTIONAL :: ConstIn
REAL,INTENT(IN),OPTIONAL :: ConstOut
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: VecOut(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)           :: i
!===================================================================================================================================
IF(    PRESENT(ConstIn).AND.(.NOT.PRESENT(ConstOut)))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)+VecIn(i)*ConstIn
  END DO
ELSEIF(PRESENT(ConstOut).AND.(.NOT.PRESENT(ConstIn)))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)*ConstOut+VecIn(i)
  END DO
ELSEIF(PRESENT(ConstIn).AND.PRESENT(ConstOut))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)*ConstOut+VecIn(i)*ConstIn
  END DO
ELSE
  DO i=1,nTotal
    VecOut=VecOut(i)+VecIn(i)
  END DO
END IF
END SUBROUTINE LVAXPBY


END MODULE MOD_Vector
