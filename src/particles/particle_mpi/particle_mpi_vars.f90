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
INTEGER,ALLOCATABLE :: PartHaloElemToProc(:,:)                               ! containing native elemid and native proc id
                                                                             ! 1 - Native_Elem_ID
                                                                             ! 2 - Rank of Proc
                                                                             ! 3 - local neighbor id

INTEGER,ALLOCATABLE :: PartHaloSideToProc(:,:)                               ! containing native sideid and native proc id
                                                                             ! 1 - Native_Side_ID
                                                                             ! 3 - Rank of Proc
                                                                             ! 4 - local neighbor id


INTEGER             :: myRealKind
LOGICAL                                  :: ParticleMPIInitIsDone=.FALSE.
LOGICAL                                  :: printMPINeighborWarnings         ! print warning messages or not
LOGICAL                                  :: printBezierControlPointsWarnings ! print warning messages or not
INTEGER                                  :: iMessage                         ! Number of MPI-Messages for Debug purpose

TYPE tPartMPIGROUP
!  TYPE(tPartMPIConnect)        , ALLOCATABLE :: MPIConnect(:)               ! MPI connect for each process
  !TYPE(MPI_Comm)                         :: COMM                            ! MPI communicator for PIC GTS region
  INTEGER                                :: COMM                             ! MPI communicator for PIC GTS region
  INTEGER                                :: nProcs                           ! number of MPI processes for particles
  INTEGER                                :: MyRank                           ! MyRank of PartMPIVAR%COMM
  LOGICAL                                :: MPIRoot                          ! Root, MPIRank=0
!  INTEGER                                :: nMPINeighbors                   ! number of MPI-Neighbors with HALO
!  LOGICAL,ALLOCATABLE                    :: isMPINeighbor(:)                ! list of possible neighbors
  INTEGER,ALLOCATABLE                    :: GroupToComm(:)                   ! list containing the rank in PartMPI%COMM
  INTEGER,ALLOCATABLE                    :: CommToGroup(:)                   ! list containing the rank in PartMPI%COMM
END TYPE

TYPE tPeriodicPtr
  INTEGER                  , ALLOCATABLE  :: BGMPeriodicBorder(:,:)          ! indices of periodic border nodes
END TYPE

#ifdef MPI
TYPE tPartMPIConnect
!  TYPE(tSidePtr)               , POINTER :: tagToSide(:)           =>NULL() ! gives side pointer for each MPI tag
  TYPE(tPeriodicPtr)       , ALLOCATABLE :: Periodic(:)                     ! data for different periodic borders for process
  LOGICAL                                :: isBGMNeighbor                    ! Flag: which process is neighber wrt. bckgrnd mesh
  LOGICAL                                :: isBGMPeriodicNeighbor            ! Flag: which process is neighber wrt. bckgrnd mesh
!  LOGICAL                      , POINTER :: myBGMPoint(:,:,:)      =>NULL()   ! Flag: does BGM point(i,j,k) belong to me?
!  LOGICAL                      , POINTER :: yourBGMPoint(:,:,:)    =>NULL()   ! Flag: does BGM point(i,j,k) belong to process?
  INTEGER                  , ALLOCATABLE :: BGMBorder(:,:)            ! indices of border nodes (1=min 2=max,xyz)
!  INTEGER                                :: nmyBGMPoints                      ! Number of BGM points in my part of the border
!  INTEGER                                :: nyourBGMPoints                    ! Number of BGM points in your part of border
  INTEGER                                :: BGMPeriodicBorderCount            ! Number(#) of overlapping areas due to periodic bc
END TYPE
#endif /*MPI*/


TYPE tPartMPIVAR
#ifdef MPI
  TYPE(tPartMPIConnect)        , ALLOCATABLE :: DepoBGMConnect(:)            ! MPI connect for each process
#endif /*MPI*/
  TYPE(tPartMPIGROUP),ALLOCATABLE        :: InitGroup(:)                     ! small communicator for initialization
  !TYPE(MPI_Comm)                         :: COMM                            ! MPI communicator for PIC GTS region
  INTEGER                                :: COMM                             ! MPI communicator for PIC GTS region
  INTEGER                                :: nProcs                           ! number of MPI processes for particles
  INTEGER                                :: MyRank                           ! MyRank of PartMPIVAR%COMM
  LOGICAL                                :: MPIRoot                          ! Root, MPIRank=0
  INTEGER                                :: nMPINeighbors                    ! number of MPI-Neighbors with HALO
  LOGICAL,ALLOCATABLE                    :: isMPINeighbor(:)                 ! list of possible neighbors
  INTEGER,ALLOCATABLE                    :: MPINeighbor(:)                   ! list containing the rank of MPI-neighbors
  INTEGER,ALLOCATABLE                    :: GlobalToLocal(:)                 ! map from global proc id to local
END TYPE

TYPE (tPartMPIVAR)                       :: PartMPI

REAL                                     :: SafetyFactor                     ! Factor to scale the halo region with MPI
REAL                                     :: halo_eps_velo                    ! halo_eps_velo
REAL                                     :: halo_eps                         ! length of halo-region
REAL                                     :: halo_eps2                        ! length of halo-region^2

#ifdef MPI
INTEGER                                  :: PartCommSize                     ! Number of REAL entries for particle communication
INTEGER                                  :: PartCommSize0                    ! Number of REAL entries for particle communication
                                                                             ! should think about own MPI-Data-Typ

TYPE tMPIMessage
  REAL,ALLOCATABLE                      :: content(:)                        ! message buffer real
  LOGICAL,ALLOCATABLE                   :: content_log(:)                    ! message buffer logical for BGM
  INTEGER,ALLOCATABLE                   :: content_int(:)                    ! message buffer for integer for adsorption
END TYPE

TYPE(tMPIMessage),ALLOCATABLE  :: PartRecvBuf(:)                             ! PartRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: PartSendBuf(:)                             ! PartSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: SurfRecvBuf(:)                             ! SurfRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: SurfSendBuf(:)                             ! SurfSendBuf with all requried types

TYPE tParticleMPIExchange
  INTEGER,ALLOCATABLE            :: nPartsSend(:,:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nPartsRecv(:,:)     ! only mpi neighbors
  INTEGER                        :: nMPIParticles     ! number of all received particles
  INTEGER,ALLOCATABLE            :: SendRequest(:,:)  ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:,:)  ! recv request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)   ! message, required for particle emission
!  INTEGER                       ,POINTER :: MPINbrOfParticles(:)
!  INTEGER                       ,POINTER :: MPIProcNbr(:)
!  INTEGER                       ,POINTER :: MPITags(:)
!  INTEGER                       ,POINTER :: nbrOfSendParticles(:,:)  ! (1:nProcs,1:2) 1: pure MPI part, 2: shape part
!  INTEGER                       ,POINTER :: NbrArray(:)  ! (1:nProcs*2)
!  INTEGER                       ,POINTER :: nbrOfSendParticlesEmission(:)  ! (1:nProcs)
END TYPE
 TYPE (tParticleMPIExchange)     :: PartMPIExchange

TYPE tParticleMPIExchange2
  INTEGER,ALLOCATABLE            :: nPartsSend(:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nPartsRecv(:)     ! only mpi neighbors
  INTEGER                        :: nMPIParticles     ! number of all received particles
  INTEGER,ALLOCATABLE            :: SendRequest(:,:)  ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:,:)  ! recv request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)   ! message, required for particle emission
END TYPE
TYPE (tParticleMPIExchange2)     :: PartMPIInsert 

TYPE tSurfMPIExchange
  INTEGER,ALLOCATABLE            :: nSidesSend(:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nSidesRecv(:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: SendRequest(:)   ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:)   ! recv request message handle,  1 - Number, 2-Message
END TYPE
TYPE (tSurfMPIExchange)          :: SurfExchange


INTEGER,ALLOCATABLE                      :: PartTargetProc(:)                ! local proc id for communication
LOGICAL,ALLOCATABLE                      :: PartMPIDepoSend(:)               ! index of part number, if particle has to be send
                                                                             ! for deposition, e.g. shape-function
LOGICAL                                  :: DoExternalParts                  ! external particles, required for 
                                                                             ! shape-function or b-spline or valume weighting
INTEGER                                  :: NbrOfExtParticles                ! number of external particles
LOGICAL                                  :: ExtPartsAllocated                ! logical,if exp parts are allocated 
REAL, ALLOCATABLE                        :: ExtPartState(:,:)                ! external particle state
INTEGER, ALLOCATABLE                     :: ExtPartSpecies(:)                ! species of external particles
INTEGER, ALLOCATABLE                     :: ExtPartToFIBGM(:,:)              ! mapping form particle to bounding box in FIBGM
REAL, ALLOCATABLE                        :: ExtPartMPF(:)                    ! macro-particle factor of external particles
INTEGER                                  :: ExtPartCommSize                  ! number of entries for ExtParticles
!INTEGER                                  :: nPartShape

REAL, ALLOCATABLE                        :: PartShiftVector(:,:)             ! store particle periodic map
#endif /*MPI*/
!===================================================================================================================================

END MODULE MOD_Particle_MPI_Vars
