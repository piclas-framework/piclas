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

PROGRAM Piclas
!===================================================================================================================================
! PICLas main program
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Piclas
USE MOD_Piclas_init       ,ONLY: FinalizePiclas
USE MOD_TimeDisc          ,ONLY: TimeDisc
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI               ,ONLY: OutputMPIW8Time
USE MOD_MPI_Vars          ,ONLY: MPIW8Time,MPIW8TimeSim,MPIW8TimeField,MPIW8TimeBaS
#if defined(PARTICLES)
USE MOD_Particle_MPI_Vars ,ONLY: MPIW8TimePart
#endif /*defined(PARTICLES)*/
#endif /*defined(MEASURE_MPI_WAIT)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if defined(MEASURE_MPI_WAIT)
MPIW8TimeSim   = 0.
MPIW8TimeBaS   = 0.
#if defined(PARTICLES)
MPIW8TimePart  = 0.
#endif /*defined(PARTICLES)*/
MPIW8TimeField = 0.
MPIW8Time      = 0.
#endif /*defined(MEASURE_MPI_WAIT)*/

! Initialize
CALL InitializePiclas()

! Run Simulation
CALL TimeDisc()

! Finalize
CALL FinalizePiclas(IsLoadBalance=.FALSE.)

#if defined(MEASURE_MPI_WAIT)
! Collect the MPI_WAIT() over all procs and output
IF(nProcessors.GT.1) CALL OutputMPIW8Time()
#endif /*defined(MEASURE_MPI_WAIT)*/

! MPI
#if USE_MPI
! We also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif
END PROGRAM Piclas
