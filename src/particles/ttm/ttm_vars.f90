MODULE MOD_TTM_Vars
!===================================================================================================================================
! Contains the TTM' variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                 :: TTMInitIsDone=.FALSE.                             !> initialization of TTM read-in
LOGICAL                 :: DoImportTTMFile                                   !> read IMD Two-Temperature Model (TTM) data (FD grid)
CHARACTER(255)          :: TTMFile                                           !> TTW Data file
! TTM-DG solution (reference / physical)
REAL,ALLOCATABLE,TARGET :: TTM(:,:,:,:,:)                                    !> TTM DG solution
! TTM-FD solution
REAL,ALLOCATABLE,TARGET :: TTM_FD(:,:,:,:)                                   !> TTM FD bary centre solution
INTEGER                 :: TTMGridFDdim(3)                                   !> number of FD grid cells in each direction
INTEGER                 :: TTMNumber                                         !> file time index
INTEGER                 :: FD_Nelems                                         !> number of TTM FD grid cells

!===================================================================================================================================
END MODULE MOD_TTM_Vars
