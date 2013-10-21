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
REAL,ALLOCATABLE                      :: D_Hat(:,:)
REAL,ALLOCATABLE                      :: L_HatMinus(:)
REAL,ALLOCATABLE                      :: L_HatPlus(:)
! DG solution
REAL,ALLOCATABLE                      :: U(:,:,:,:,:)
! DG time derivative
REAL,ALLOCATABLE                      :: Ut(:,:,:,:,:)
! number of array items in U, Ut, gradUx, gradUy, gradUz after allocated
INTEGER                               :: nTotalU
! interior face values for all elements
REAL,ALLOCATABLE                      :: U_Minus(:,:,:,:),U_Plus(:,:,:,:)
REAL,ALLOCATABLE                      :: Flux(:,:,:,:)
LOGICAL                               :: DGInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_DG_Vars
