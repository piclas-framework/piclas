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

MODULE MOD_Particle_MPI_Vars
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
!USE mpi
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTEGER                                   :: nExchangeProcessors            ! number of MPI processes for particles exchange
INTEGER,ALLOCATABLE                       :: ExchangeProcToGlobalProc(:,:)  ! mapping from exchange proc ID to global proc ID
INTEGER,ALLOCATABLE                       :: GlobalProcToExchangeProc(:,:)  ! mapping from global proc ID to exchange proc ID
LOGICAL                                   :: ParticleMPIInitIsDone=.FALSE.

TYPE tPartMPIGROUP
  INTEGER                                 :: COMM                           ! MPI communicator for PIC GTS region
  INTEGER                                 :: Request                        ! MPI request for asynchronous communication
  INTEGER                                 :: nProcs                         ! number of MPI processes for particles
  INTEGER                                 :: MyRank                         ! MyRank of PartMPIVAR%COMM
  LOGICAL                                 :: MPIRoot                        ! Root, MPIRank=0
  INTEGER,ALLOCATABLE                     :: GroupToComm(:)                 ! list containing the rank in PartMPI%COMM
  INTEGER,ALLOCATABLE                     :: CommToGroup(:)                 ! list containing the rank in PartMPI%COMM
END TYPE

TYPE tPeriodicPtr
  INTEGER                  , ALLOCATABLE  :: BGMPeriodicBorder(:,:)         ! indices of periodic border nodes
END TYPE

#if USE_MPI
TYPE tPartMPIConnect
  TYPE(tPeriodicPtr)       , ALLOCATABLE  :: Periodic(:)                    ! data for different periodic borders for process
  LOGICAL                                 :: isBGMNeighbor                  ! Flag: which process is neighber wrt. bckgrnd mesh
  LOGICAL                                 :: isBGMPeriodicNeighbor          ! Flag: which process is neighber wrt. bckgrnd mesh
  INTEGER                  , ALLOCATABLE  :: BGMBorder(:,:)                 ! indices of border nodes (1=min 2=max,xyz)
  INTEGER                                 :: BGMPeriodicBorderCount         ! Number(#) of overlapping areas due to periodic bc
END TYPE
#endif /*USE_MPI*/

TYPE tPartMPIVAR
#if USE_MPI
  TYPE(tPartMPIConnect)    , ALLOCATABLE  :: DepoBGMConnect(:)              ! MPI connect for each process
#endif /*USE_MPI*/
  TYPE(tPartMPIGROUP),ALLOCATABLE         :: InitGroup(:)                   ! small communicator for initialization
  INTEGER                                 :: COMM                           ! MPI communicator for PIC GTS region
  INTEGER                                 :: nProcs                         ! number of MPI processes for particles
  INTEGER                                 :: MyRank                         ! MyRank of PartMPIVAR%COMM
  LOGICAL                                 :: MPIRoot                        ! Root, MPIRank=0
END TYPE

TYPE (tPartMPIVAR)                        :: PartMPI

REAL                                      :: SafetyFactor                   ! Factor to scale the halo region with MPI
REAL                                      :: halo_eps_velo                  ! halo_eps_velo
REAL                                      :: halo_eps                       ! length of halo-region
REAL                                      :: halo_eps2                      ! length of halo-region^2

#if USE_MPI
INTEGER                                   :: PartCommSize                   ! Number of REAL entries for particle communication
INTEGER                                   :: PartCommSize0                  ! Number of REAL entries for particle communication
                                                                            ! should think about own MPI-Data-Typ

TYPE tMPIMessage
  REAL,ALLOCATABLE                        :: content(:)                     ! message buffer real
  LOGICAL,ALLOCATABLE                     :: content_log(:)                 ! message buffer logical for BGM
  INTEGER,ALLOCATABLE                     :: content_int(:)                 ! message buffer for integer
END TYPE

TYPE(tMPIMessage),ALLOCATABLE             :: PartRecvBuf(:)                 ! PartRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE             :: PartSendBuf(:)                 ! PartSendBuf with all required types

TYPE(tMPIMessage),ALLOCATABLE             :: SurfRecvBuf(:)                 ! SurfRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE             :: SurfSendBuf(:)                 ! SurfSendBuf with all required types

TYPE(tMPIMessage),ALLOCATABLE             :: PorousBCRecvBuf(:)             ! SurfRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE             :: PorousBCSendBuf(:)             ! SurfSendBuf with all required types

TYPE(tMPIMessage),ALLOCATABLE             :: EmissionRecvBuf(:)             ! EmissionRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE             :: EmissionSendBuf(:)             ! EmissionSendBuf with all required types

TYPE tParticleMPIExchange
  INTEGER,ALLOCATABLE                     :: nPartsSend(:,:)                ! Only MPI neighbors
  INTEGER,ALLOCATABLE                     :: nPartsRecv(:,:)                ! Only MPI neighbors
  INTEGER                                 :: nMPIParticles                  ! Number of all received particles
  INTEGER,ALLOCATABLE                     :: SendRequest(:,:)               ! Send request message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE                     :: RecvRequest(:,:)               ! Receive request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE           :: send_message(:)                ! Message, required for particle emission
END TYPE

TYPE (tParticleMPIExchange)               :: PartMPIExchange
TYPE (tParticleMPIExchange)               :: PartMPIInsert
TYPE (tParticleMPIExchange)               :: PartMPILocate

INTEGER,ALLOCATABLE                       :: PartTargetProc(:)              ! local rank id for communication
REAL, ALLOCATABLE                         :: PartShiftVector(:,:)           ! store particle periodic map
#endif /*USE_MPI*/
!===================================================================================================================================

END MODULE MOD_Particle_MPI_Vars
