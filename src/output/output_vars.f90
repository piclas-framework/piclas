MODULE MOD_Output_Vars
!===================================================================================================================================
! Contains global variables provided by the output routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                      :: NVisu                        ! number of visualisation points is NVisu+1
REAL,ALLOCATABLE             :: Vdm_GaussN_NVisu(:,:)        ! for direct interpolation from computation grid to visu grid
REAL,PARAMETER               :: FileVersion=0.1
CHARACTER(LEN=255),PARAMETER :: ProgramName='Boltzplatz'
CHARACTER(LEN=255)           :: ProjectName
INTEGER                      :: outputFormat=0           ! =0: visualization off, >0 visualize
LOGICAL                      :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Output_Vars
