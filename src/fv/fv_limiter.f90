!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!> Module for slope limiters of FV sub-cells.
!==================================================================================================================================
MODULE MOD_FV_Limiter
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: FV_LIMITERTYPE_NULL    = 0
INTEGER,PARAMETER :: FV_LIMITERTYPE_MINMOD  = 1
INTEGER,PARAMETER :: FV_LIMITERTYPE_SWEBY   = 2
INTEGER,PARAMETER :: FV_LIMITERTYPE_VANLEER = 3
INTEGER,PARAMETER :: FV_LIMITERTYPE_VKT = 4
INTEGER,PARAMETER :: FV_LIMITERTYPE_CENTRAL = 9

INTERFACE InitFV_Limiter
  MODULE PROCEDURE InitFV_Limiter
END INTERFACE

ABSTRACT INTERFACE
  PPURE SUBROUTINE LimiterInt(sL, sR, s)
    REAL,INTENT(IN)  :: sL(:),sR(:)
    REAL,INTENT(OUT) :: s(:)
  END SUBROUTINE
END INTERFACE

PROCEDURE(LimiterInt),POINTER :: FV_Limiter !< limiting function (see: fv_limiter.f90)

PUBLIC::InitFV_Limiter
PUBLIC::FV_Limiter
!==================================================================================================================================


CONTAINS

!==================================================================================================================================
!> Initialize pointer to chosen limiter type and readin of required parameter
!==================================================================================================================================
SUBROUTINE InitFV_Limiter()
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_FV_Vars     ,ONLY: LimiterType,FV_sweby_beta
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Read in LimiterType parameters
SELECT CASE(LimiterType)
CASE (FV_LIMITERTYPE_NULL) ! NullLimiter
  FV_Limiter => NullLimiter
  SWRITE(UNIT_stdOut,'(A)') '  Using "NullLimiter" limiter.'
CASE (FV_LIMITERTYPE_MINMOD) ! MinMod
  FV_Limiter => MinMod
  SWRITE(UNIT_stdOut,'(A)') '  Using "Minmod" limiter.'
CASE (FV_LIMITERTYPE_SWEBY) ! Sweby
  FV_Limiter => Sweby
  FV_sweby_beta = GETREAL('swebyb')
  SWRITE(UNIT_stdOut,'(A,F8.6)') '  Using "Sweby" limiter with beta = ', FV_sweby_beta
CASE (FV_LIMITERTYPE_VANLEER) ! Van Leer
  FV_Limiter => VanLeer
  SWRITE(UNIT_stdOut,'(A)') '  Using "van Leer" limiter.'
CASE (FV_LIMITERTYPE_VKT) ! Venkatakrishnan
  FV_Limiter => VanLeer
  SWRITE(UNIT_stdOut,'(A)') '  Using "Venkatakrishnan" limiter.'
CASE (FV_LIMITERTYPE_CENTRAL) ! Central
  FV_Limiter => CentralLimiter
  SWRITE(UNIT_stdOut,'(A)') '  Using "Central" limiter.'
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'FV Limiter-Type unknown.')
END SELECT
END SUBROUTINE InitFV_Limiter

!==================================================================================================================================
!> Limiter sets slope to zero.
!==================================================================================================================================
PPURE SUBROUTINE NullLimiter(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(:) !< left slope
REAL,INTENT(IN)  :: sR(:) !< right slope
REAL,INTENT(OUT) :: s(:)  !< limited slope
!==================================================================================================================================
! NullLimiter
s = 0.
END SUBROUTINE NullLimiter

!==================================================================================================================================
!> MinMod slope limiter.
!==================================================================================================================================
PPURE SUBROUTINE MinMod(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(:) !< left slope
REAL,INTENT(IN)  :: sR(:) !< right slope
REAL,INTENT(OUT) :: s(:)  !< limited slope
!==================================================================================================================================
! MinMod
s = MERGE(sL,sR, ABS(sL) .LT. ABS(sR))
s = MERGE(s,0., sL*sR .GT. 0.)
!if (sL(1)*sR(1) .LT.0) then
  !if (abs(sl(1)*sR(1)).gt.0.00001) then
     !WRITE (*,*) sl(1), sr(1), s(1)
  !end if
!end if
END SUBROUTINE MinMod

!==================================================================================================================================
!> Sweby slope limiter.
!==================================================================================================================================
PPURE SUBROUTINE Sweby(sL, sR, s)
! MODULES
USE MOD_PreProc ,ONLY: PP_nVar
USE MOD_FV_Vars ,ONLY: FV_sweby_beta
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(:) !< left slope
REAL,INTENT(IN)  :: sR(:) !< right slope
REAL,INTENT(OUT) :: s(:)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: sa(PP_nVar),sb(PP_nVar)
!==================================================================================================================================
CALL MinMod(sL*FV_sweby_beta,sR,sa)
CALL MinMod(sL,sR*FV_sweby_beta,sb)
s = SIGN(MAX(ABS(sa),ABS(sb)),sL)
END SUBROUTINE Sweby

!==================================================================================================================================
!> central limiter s = (sL + sR)/2  (ATTENTION: unstable and not TVD)
!==================================================================================================================================
PPURE SUBROUTINE CentralLimiter(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(:) !< left slope
REAL,INTENT(IN)  :: sR(:) !< right slope
REAL,INTENT(OUT) :: s(:)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
s = (sL+sR)/2.0
END SUBROUTINE CentralLimiter

!==================================================================================================================================
!> van Leer limiter
!==================================================================================================================================
PPURE SUBROUTINE VanLeer(sL, sR, s)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: sL(:) !< left slope
REAL,INTENT(IN)  :: sR(:) !< right slope
REAL,INTENT(OUT) :: s(:)  !< limited slope
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
s = (SIGN(1.,sL)+SIGN(1.,sR))*ABS(sL)*ABS(sR)/(ABS(sL)+ABS(sR)+EPSILON(sL))
END SUBROUTINE VanLeer


END MODULE MOD_FV_Limiter
