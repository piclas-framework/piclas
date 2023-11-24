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

!===================================================================================================================================
!> Message Passing Interface (MPI) variables
!===================================================================================================================================
MODULE MOD_MPI_Vars
#if USE_MPI
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: SendRequest_U(:),SendRequest_Flux(:),SendRequest_gradUx(:),SendRequest_gradUy(:),SendRequest_gradUz(:)
INTEGER,ALLOCATABLE :: SendRequest_U2(:),RecRequest_U2(:)
INTEGER,ALLOCATABLE :: RecRequest_U(:),RecRequest_Flux(:),RecRequest_gradUx(:),RecRequest_gradUy(:),RecRequest_gradUz(:)
INTEGER,ALLOCATABLE :: SendRequest_Geo(:),RecRequest_Geo(:)
INTEGER             :: iNbProc
INTEGER             :: nSendVal,nRecVal,DataSizeSide
INTEGER             :: SideID_start,SideID_end
#if USE_MPI
INTEGER               :: nNbProcs         ! number of neighbor procs
INTEGER,ALLOCATABLE   :: NbProc(:)        ! iProc list of neighbor procs; allocated from 1:nNbProcs
INTEGER,ALLOCATABLE   :: nMPISides_Proc(:)
INTEGER,ALLOCATABLE   :: nMPISides_MINE_Proc(:)
INTEGER,ALLOCATABLE   :: nMPISides_YOUR_Proc(:)
INTEGER,ALLOCATABLE   :: offsetMPISides_MINE(:)! gives position of send/recv block in *_MINE arrays,allocated from 0:nNbProcs
INTEGER,ALLOCATABLE   :: offsetMPISides_YOUR(:)! gives position of send/recv block in *_YOUR arrays,allocated from 0:nNbProcs
INTEGER,ALLOCATABLE   :: offsetElemMPI(:)      ! gives offset position of elements of all procs
INTEGER,ALLOCATABLE   :: nMPISides_send(:,:),nMPISides_rec(:,:)
INTEGER,ALLOCATABLE   :: OffsetMPISides_send(:,:),OffsetMPISides_rec(:,:)
#endif /*USE_MPI*/

#if defined(MEASURE_MPI_WAIT)
! Elapsed times
REAL(KIND=8)             :: MPIW8TimeSim                   !< measure global time in the simulation as reference
REAL(KIND=8)             :: MPIW8TimeBaS                   !< measure time on each proc it is in BARRIER_AND_SYNC
REAL(KIND=8)             :: MPIW8TimeMM                    !< measure time on each proc it is in REDUCE for RAM measurement
REAL(KIND=8)             :: MPIW8TimeField(MPIW8SIZEFIELD) !< measure time on each proc it is in MPI_WAIT() during the field solver
REAL(KIND=8)             :: MPIW8Time(MPIW8SIZE)           !< measure time on each proc it is in MPI_WAIT()
REAL(KIND=8)             :: MPIW8TimeGlobal(MPIW8SIZE)     !< measure time on each proc it is in MPI_WAIT() global over all ranks
REAL(KIND=8),ALLOCATABLE :: MPIW8TimeProc(:)               !< measure time on each proc it is in MPI_WAIT() proc local output
! Counter
INTEGER(KIND=8)             :: MPIW8CountBaS                   !< count the number of measurements on each proc it is in BARRIER_AND_SYNC
INTEGER(KIND=8)             :: MPIW8CountMM                    !< count the number of measurements on each proc it is in REDUCE for RAM measurement
INTEGER(KIND=8)             :: MPIW8CountField(MPIW8SIZEFIELD) !< count the number of measurements on each proc it is in MPI_WAIT() during the field solver
INTEGER(KIND=8)             :: MPIW8Count(MPIW8SIZE)           !< count the number of measurements on each proc it is in MPI_WAIT()
INTEGER(KIND=8)             :: MPIW8CountGlobal(MPIW8SIZE)     !< count the number of measurements on each proc it is in MPI_WAIT() global over all ranks
INTEGER(KIND=8),ALLOCATABLE :: MPIW8CountProc(:)               !< count the number of measurements on each proc it is in MPI_WAIT() proc local output
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#endif
END MODULE MOD_MPI_Vars
