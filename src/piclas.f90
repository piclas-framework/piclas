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

PROGRAM Piclas
!===================================================================================================================================
! Control program of the Piclas code.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Piclas
USE MOD_Piclas_init ,ONLY: FinalizePiclas
USE MOD_TimeDisc    ,ONLY: TimeDisc
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Time
LOGICAL                 :: userblockFound
!===================================================================================================================================
! Initialize
CALL InitializePiclas()

! Run Simulation
CALL TimeDisc()

! Finalize
CALL FinalizePiclas(IsLoadBalance=.FALSE.)

#if USE_MPI
! we also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif
END PROGRAM Piclas
