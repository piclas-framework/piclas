#include "boltzplatz.h"

MODULE MOD_SurfaceModel_MPI_Vars
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
#ifdef MPI
TYPE tMPIMessage
  REAL,ALLOCATABLE             :: content(:)                                 ! message buffer real
  LOGICAL,ALLOCATABLE          :: content_log(:)                             ! message buffer logical for BGM
  INTEGER,ALLOCATABLE          :: content_int(:)                             ! message buffer for integer for adsorption
END TYPE

TYPE(tMPIMessage),ALLOCATABLE  :: SurfRecvBuf(:)                             ! SurfRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: SurfSendBuf(:)                             ! SurfSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: SurfDistRecvBuf(:)                         ! SurfDistRecvBuf with all requried types
TYPE(tMPIMessage),ALLOCATABLE  :: SurfDistSendBuf(:)                         ! SurfDistSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: AdsorbRecvBuf(:)
TYPE(tMPIMessage),ALLOCATABLE  :: AdsorbSendBuf(:)

TYPE(tMPIMessage),ALLOCATABLE  :: SurfCoverageRecvBuf(:)
TYPE(tMPIMessage),ALLOCATABLE  :: SurfCoverageSendBuf(:)

TYPE(tMPIMessage),ALLOCATABLE  :: CondensRecvBuf(:)
TYPE(tMPIMessage),ALLOCATABLE  :: CondensSendBuf(:)

TYPE tDistNbrComm
  INTEGER,ALLOCATABLE :: nPosSend(:) ! number of distribution site to send for surface (nSidesSend,nCoordination)
  INTEGER,ALLOCATABLE :: nPosRecv(:) ! number of distribution sites to receive for surface (nSidesRecv,nCoordination)
END TYPE

TYPE tSurfModelMPIExchange
  INTEGER,ALLOCATABLE            :: nSidesSend(:)            ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nSidesRecv(:)            ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: SendRequest(:)           ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:)           ! recv request message handle,  1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: nSurfDistSidesSend(:)    ! number of mpi sides to send (nProcs)
  INTEGER,ALLOCATABLE            :: nSurfDistSidesRecv(:)    ! number of sides received from mpi (nProcs)
  INTEGER,ALLOCATABLE            :: nCoverageSidesSend(:)    ! number of mpi sides to send (nProcs)
  INTEGER,ALLOCATABLE            :: nCoverageSidesRecv(:)    ! number of sides received from mpi (nProcs)
  INTEGER,ALLOCATABLE            :: SurfDistSendRequest(:,:) ! send request message handle,  1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: SurfDistRecvRequest(:,:) ! recv request message handle,  1 - Number, 2-Message
  TYPE(tDistNbrComm),ALLOCATABLE :: NbrOfPos(:)              ! array for number of distribution sites sending per proc 
END TYPE
TYPE (tSurfModelMPIExchange)          :: SurfModelExchange

#endif /*MPI*/
!===================================================================================================================================

END MODULE MOD_SurfaceModel_MPI_Vars
