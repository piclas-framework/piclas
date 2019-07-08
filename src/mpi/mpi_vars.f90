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
MODULE MOD_MPI_Vars
#ifdef MPI
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
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
LOGICAL             :: MPIInitIsDone=.FALSE.
#ifdef MPI
INTEGER               :: nNbProcs         ! number of neighbor procs
INTEGER,ALLOCATABLE   :: NbProc(:)        ! iProc list of neighbor procs; allocated from 1:nNbProcs
INTEGER,ALLOCATABLE   :: nMPISides_Proc(:)
INTEGER,ALLOCATABLE   :: nMPISides_MINE_Proc(:)
INTEGER,ALLOCATABLE   :: nMPISides_YOUR_Proc(:)
INTEGER,ALLOCATABLE   :: offsetMPISides_MINE(:)! gives position of send/recv block in *_MINE arrays,allocated from 0:nNbProcs
INTEGER,ALLOCATABLE   :: offsetMPISides_YOUR(:)! gives position of send/recv block in *_YOUR arrays,allocated from 0:nNbProcs
INTEGER,ALLOCATABLE   :: offsetElemMPI(:)      ! gives offsetposotion of elements of all procs
INTEGER,ALLOCATABLE   :: nMPISides_send(:,:),nMPISides_rec(:,:)
INTEGER,ALLOCATABLE   :: OffsetMPISides_send(:,:),OffsetMPISides_rec(:,:)
#endif /*MPI*/
!===================================================================================================================================
#endif
END MODULE MOD_MPI_Vars
