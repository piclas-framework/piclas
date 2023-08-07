!==================================================================================================================================
! Copyright (c) 2023 - 2023 Stephen Copplestone
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

MODULE MOD_Equation_Tools
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_MPI && defined(PARTICLES) && USE_HDG
PUBLIC :: SynchronizeCPP
#endif /*USE_MPI && defined(PARTICLES) && USE_HDG*/
!===================================================================================================================================
CONTAINS

#if USE_MPI && defined(PARTICLES) && USE_HDG
!===================================================================================================================================
!> Communicate the CPP values from MPIRoot to sub-communicator processes
!===================================================================================================================================
SUBROUTINE SynchronizeCPP()
! MODULES
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
USE MOD_HDG_Vars ,ONLY: CPPCOMM
USE MOD_HDG_Vars ,ONLY: CoupledPowerPotential,CPPDataLength
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(CPPCOMM%UNICATOR.NE.MPI_COMM_NULL)THEN
  ! Broadcast from root to other processors on the sub-communicator
  CALL MPI_BCAST(CoupledPowerPotential, CPPDataLength, MPI_DOUBLE_PRECISION, 0, CPPCOMM%UNICATOR, IERROR)
END IF
END SUBROUTINE SynchronizeCPP
#endif /*USE_MPI && defined(PARTICLES) && USE_HDG*/

END MODULE MOD_Equation_Tools
