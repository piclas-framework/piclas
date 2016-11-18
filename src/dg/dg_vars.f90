MODULE MOD_DG_Vars
!===================================================================================================================================
! Contains global variables used by the DG modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! DG basis
REAL,ALLOCATABLE                      :: D(:,:)
REAL,ALLOCATABLE                      :: D_T(:,:)    ! D^T
REAL,ALLOCATABLE                      :: D_Hat(:,:)
REAL,ALLOCATABLE                      :: D_Hat_T(:,:) ! D_Hat^T
REAL,ALLOCATABLE                      :: L_HatMinus(:)
REAL,ALLOCATABLE                      :: L_HatPlus(:)
! DG solution (reference / physical)
REAL,ALLOCATABLE,TARGET               :: U(:,:,:,:,:) ! computed from JU
! DG time derivative
REAL,ALLOCATABLE                      :: Ut(:,:,:,:,:)
! number of array items in U, Ut, gradUx, gradUy, gradUz after allocated
INTEGER                               :: nTotalU
INTEGER                               :: nTotal_vol    !loop i,j,k
INTEGER                               :: nTotal_face   !loop i,j
! interior face values for all elements
REAL,ALLOCATABLE                      :: U_master(:,:,:,:),U_slave(:,:,:,:)
REAL,ALLOCATABLE                      :: Flux(:,:,:,:)
LOGICAL                               :: DGInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_DG_Vars
