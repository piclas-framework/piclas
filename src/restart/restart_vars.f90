MODULE MOD_Restart_Vars
!===================================================================================================================================
! Contains global variables used by the restart module
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE   :: Vdm_GaussNRestart_GaussN(:,:)! for interpolation from restart grid to computation grid
INTEGER            :: nVar_Restart
INTEGER            :: N_Restart = 0
INTEGER            :: nElems_Restart
LOGICAl            :: RestartInitIsDone   = .FALSE.
LOGICAl            :: DoRestart           = .FALSE.
LOGICAl            :: DoInitialAutoRestart= .FALSE.
INTEGER            :: InitialAutoRestartSample
LOGICAl            :: BuildNewMesh        = .TRUE.
LOGICAl            :: WriteNewMesh        = .TRUE.
LOGICAL            :: InterpolateSolution =.FALSE.
CHARACTER(LEN=300) :: RestartFile = ""
CHARACTER(LEN=255) :: NodeType_Restart
REAL               :: RestartTime
REAL               :: RestartWallTime ! wall time at the beginning of a simulation OR when a restart is performed via Load Balance
!===================================================================================================================================
END MODULE MOD_Restart_Vars
