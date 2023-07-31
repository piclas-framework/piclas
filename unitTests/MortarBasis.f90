!==================================================================================================================================
! Copyright (c) 2023 Stephen M. Copplestone, Patrick Kopper
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

!==================================================================================================================================
!> Unit test 'MortarBasisUnitTest'
!> Test the module: MOD_Mortar, function InitMortar
!==================================================================================================================================
PROGRAM MortarBasisUnitTest
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mortar        ,ONLY: InitMortar,FinalizeMortar
USE MOD_MPI           ,ONLY: InitMPI
USE MOD_Interpolation ,ONLY: InitInterpolation,FinalizeInterpolation
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nArgs,i
INTEGER,PARAMETER :: MaxPolDeg=50
!==================================================================================================================================
CALL InitMPI()
! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) CALL abort(__STAMP__,'ERROR - Unknown command line argument.')

#if PP_N == N
DO i = 1, MaxPolDeg
  PP_N = i
#endif
  WRITE(*,'(A,I0)') " Testing MortarInit for N = ",PP_N
  CALL InitInterpolation(PP_N,2*(PP_N+1))
  CALL InitMortar()
  CALL FinalizeMortar()
  CALL FinalizeInterpolation()
#if PP_N == N
END DO ! i = 1, MaxPolDeg
#endif

#if USE_MPI
! we also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError.NE.0) CALL abort(__STAMP__,'MPI finalize error')
#endif

END PROGRAM MortarBasisUnitTest
