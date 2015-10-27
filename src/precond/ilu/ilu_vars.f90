MODULE MOD_ILU_Vars
!===================================================================================================================================
! Contains global variables used by the Jac_Ex module.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: DiagBILU0
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)   :: XiBILU0,EtaBILU0,ZetaBILU0
INTEGER                                 :: nBlockEntries,nBDOF
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: BlockAA
INTEGER,ALLOCATABLE,DIMENSION(:)        :: BlockIA,BlockJA,BlockDiag
!===================================================================================================================================
END MODULE MOD_ILU_Vars
