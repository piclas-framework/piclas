MODULE MOD_Filter_Vars
!===================================================================================================================================
! Contains global variables used for/by the filter module
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE   :: Vdm_Leg(:,:), sVdm_Leg(:,:) ! 1D Vandermondematrix to Legendre polynomials and its inverse
INTEGER            :: FilterType
REAL               :: HestFilterParam(3)
REAL,ALLOCATABLE   :: FilterMat(:,:)        ! 1D nodal filter matrix
LOGICAl            :: FilterInitIsDone = .FALSE.
!===================================================================================================================================
END MODULE MOD_Filter_Vars
