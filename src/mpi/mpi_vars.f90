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
