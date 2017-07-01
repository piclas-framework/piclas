#ifdef PARTICLES
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
LOGICAL                    :: TTMInitIsDone=.FALSE.                             !> initialization of TTM read-in
LOGICAL                    :: DoImportTTMFile                                   !> read IMD Two-Temperature Model (TTM) data (FD grid)
CHARACTER(255)             :: TTMFile                                           !> TTW Data file
CHARACTER(255)             :: TTMLogFile                                        !> TTW Data file
! TTM-DG solution (reference / physical)
REAL,ALLOCATABLE,TARGET    :: TTM_DG(:,:,:,:,:)                                 !> TTM DG solution
! TTM-FD solution
REAL,ALLOCATABLE,TARGET    :: TTM_FD(:,:,:,:)                                   !> TTM FD bary center solution
REAL,ALLOCATABLE,TARGET    :: ElemBaryFD(:,:)                                   !> TTM FD bary center position
INTEGER,ALLOCATABLE,TARGET :: ElemIndexFD(:,:)                                  !> TTM FD index position of FD cell in DG grid
REAL                       :: TTMElemBaryTolerance                              !> TTM FD bary center tolerance to DG bary center
LOGICAL,ALLOCATABLE,TARGET :: ElemIsDone(:)                                     !> TTM FD bary center has been found
INTEGER                    :: TTMGridFDdim(3)                                   !> number of FD grid cells in each direction
INTEGER                    :: TTMNumber                                         !> file time index
INTEGER                    :: FD_nElems                                         !> number of TTM FD grid cells

!===================================================================================================================================
END MODULE MOD_TTM_Vars
#endif /*PARTICLES*/
