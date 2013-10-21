MODULE MOD_part_MPI_Vars
!===================================================================================================================================
! Contains the Particles' variables (general for all modules: PIC, DSMC, FP)
!===================================================================================================================================
! MODULES
USE MOD_mesh_vars, ONLY : tSidePtr   ! THIS NEEDS TO BE REMOVED!!!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE tMPIMessage
  REAL                         , POINTER :: content(:)             =>NULL()   ! message buffer real
  LOGICAL                      , POINTER :: content_log(:)         =>NULL()   ! message buffer logical
END TYPE

TYPE tPeriodicPtr
!  LOGICAL                      , POINTER :: myBGMPeriodicPoint(:,:,:)   =>NULL() ! same as below for periodic borders +#
!  LOGICAL                      , POINTER :: yourBGMPeriodicPoint(:,:,:) =>NULL() ! same as below for periodic borders +#
  INTEGER                      , POINTER :: BGMPeriodicBorder(:,:)          ! indices of periodic border nodes
!  INTEGER                                :: nmyBGMPeriodicPoints           ! same as below for periodic border for each border
!  INTEGER                                :: nyourBGMPeriodicPoints         ! same as below for periodic border for each border
END TYPE

TYPE tPartMPIConnect
!  TYPE(tSidePtr)               , POINTER :: tagToSide(:)           =>NULL()   ! gives side pointer for each MPI tag
  TYPE(tPeriodicPtr)       , ALLOCATABLE :: Periodic(:)                       ! data for different periodic borders for process
  LOGICAL                                :: isBGMNeighbor                     ! Flag: which process is neighber wrt. bckgrnd mesh
  LOGICAL                                :: isBGMPeriodicNeighbor             ! Flag: which process is neighber wrt. bckgrnd mesh
!  LOGICAL                      , POINTER :: myBGMPoint(:,:,:)      =>NULL()   ! Flag: does BGM point(i,j,k) belong to me?
!  LOGICAL                      , POINTER :: yourBGMPoint(:,:,:)    =>NULL()   ! Flag: does BGM point(i,j,k) belong to process?
  INTEGER                  , ALLOCATABLE :: BGMBorder(:,:)            ! indices of border nodes (1=min 2=max,xyz)
!  INTEGER                                :: nmyBGMPoints                      ! Number of BGM points in my part of the border
!  INTEGER                                :: nyourBGMPoints                    ! Number of BGM points in your part of border
  INTEGER                                :: BGMPeriodicBorderCount            ! Number(#) of overlapping areas due to periodic bc
END TYPE

TYPE tPartMPIVAR
  TYPE(tPartMPIConnect)        , ALLOCATABLE :: MPIConnect(:)             ! MPI connect for each process
  INTEGER                                :: COMM                              ! MPI communicator for PIC GTS region
  INTEGER                                :: GROUP                             ! MPI group for PIC GTS region
  INTEGER                                :: GROUPWORLD                        ! MPI group of MPI_COMM_WORLD
!  REAL                                   :: PTime,MTime,PMPITime,MMPITime,tmpTime
!  LOGICAL                                :: capo
!  INTEGER                      , POINTER :: MPIproclist(:)         =>NULL()   ! List of MPI ranks connected in PIC
  INTEGER                                :: iProc                             ! MPIRank for GTSRegions communicator
  INTEGER                                :: nProcs                            ! MPIRank for GTSRegions communicator
  INTEGER                       ,POINTER :: cango_send_request(:)
  INTEGER                       ,POINTER :: cango_recv_request(:)
  LOGICAL                       ,POINTER :: MPINeighbor(:)
END TYPE

TYPE (tPartMPIVAR)                       :: PMPIVAR

TYPE tPartMPIExchange
  INTEGER                       ,POINTER :: MPINbrOfParticles(:)
  INTEGER                       ,POINTER :: MPIProcNbr(:)
  INTEGER                       ,POINTER :: MPITags(:)
  TYPE(tMPIMessage)             ,POINTER :: send_message(:)
  INTEGER                       ,POINTER :: send_request(:,:)
  INTEGER                       ,POINTER :: recv_request(:,:)
  INTEGER                       ,POINTER :: nbrOfSendParticles(:,:)  ! (1:nProcs,1:2) 1: pure MPI part, 2: shape part
  INTEGER                       ,POINTER :: NbrArray(:)  ! (1:nProcs*2)
  INTEGER                       ,POINTER :: nbrOfSendParticlesEmission(:)  ! (1:nProcs)
  INTEGER                       ,POINTER :: nbrOfRecvParticlesEmission(:)  ! (1:nProcs)
  INTEGER                                :: nMPIParticles
END TYPE

TYPE (tPartMPIExchange)                  :: PMPIInsert
TYPE (tPartMPIExchange)                  :: PMPIExchange

TYPE tSurfMPIExchange
  INTEGER                       ,POINTER :: MPIProcNbr(:)
  INTEGER                       ,POINTER :: MPITags(:)
  TYPE(tMPIMessage)             ,POINTER :: send_message(:)
  INTEGER                       ,POINTER :: NbrArray(:)     
  INTEGER                       ,POINTER :: send_request(:,:)
  INTEGER                       ,POINTER :: recv_request(:,:)
END TYPE

TYPE (tSurfMPIExchange)                  :: SurfMPIExchang

TYPE tMPIGeometry
  INTEGER, ALLOCATABLE                   :: haloMPINbSide(:)                  ! maps from ourSide to halo-cell side at bound
                                                                              ! (1:nMPISides)
  INTEGER, ALLOCATABLE                   :: ElemToNodeID(:,:)                 ! ElemToNodeID(1:nElemNodes,1:nElems)
  INTEGER, ALLOCATABLE                   :: PeriodicElemSide(:,:)             ! (1:6,1:nElems), see GEO
  LOGICAL, ALLOCATABLE                   :: ConcaveElemSide(:,:)              ! Side Concave? See GEO
  INTEGER, ALLOCATABLE                   :: ElemSideNodeID(:,:,:)             ! ElemSideNodeID(1:nSideNodes,1:nLocSides,1:nElems)
                                                                              ! From element sides to node IDs
  REAL,    ALLOCATABLE                   :: NodeCoords(:,:)                   ! Node Coordinates (1:nDim,1:nNodes)
  INTEGER, ALLOCATABLE                   :: ElemToSide(:,:,:)                 ! (1:2,1:nLocSides,1:nElems) 
  INTEGER, ALLOCATABLE                   :: SideToElem(:,:)                   ! (1:5,1:nSides)
  INTEGER, ALLOCATABLE                   :: BC(:,:)                           ! (1:4,1:nSides)
                                                                              ! 1: BC / -1 for MPI-Side / 424242 for halo-boundary
                                                                              ! 2: 0 / if MPI-Side: NbProc
                                                                              ! 3: 0 / if MPI-Side: 1:MINE, 2:YOUR
                                                                              ! 4: 0 / if MPI-Side: MPI-Tag (init) or 
                                                                              !        (-1)*connection-Side if halo-MPI-side 
                                                                              ! 5: 0 / if MPI-Side: from which process am I?
  INTEGER, ALLOCATABLE                   :: NativeElemID(:)                   ! (1:nElems)
                                                                              ! nElems(0:PMPIVAR%nProcs-1)
  INTEGER, ALLOCATABLE                   :: ElemMPIID(:)                      ! (1:nElems)
END TYPE

TYPE (tMPIGeometry), POINTER             :: MPIGEO                            ! Container for HALO data

INTEGER, ALLOCATABLE                     :: ElemToSide(:,:,:)
INTEGER, ALLOCATABLE                     :: SideToElem(:,:)
INTEGER, ALLOCATABLE                     :: casematrix(:,:)                   ! matrix to compute periodic cases
INTEGER                                  :: NbrOfCases                        ! Number of periodic cases
REAL, ALLOCATABLE                        :: partShiftVector(:,:)              ! tells periodic shift vector for each particle

INTEGER                                  :: NbrOfAllocatedExtParts
INTEGER                                  :: NbrOfExtParticles
LOGICAL                                  :: ExtPartsAllocated

REAL, ALLOCATABLE                        :: ExtPartState(:,:)
INTEGER, ALLOCATABLE                        :: ExtPartSpecies(:)

INTEGER                                  :: FIBGMCellPadding(1:3)
REAL, PARAMETER                          :: SafetyFactor = 1.0   
REAL                                     :: halo_eps               ! halo region radius
!===================================================================================================================================
END MODULE MOD_part_MPI_Vars
