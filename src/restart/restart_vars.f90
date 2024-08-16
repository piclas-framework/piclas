!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
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
LOGICAl            :: FlushInitialState   = .FALSE. !< During restart delete the state restart file when FlushInitialState=T
LOGICAl            :: DoInitialAutoRestart= .FALSE.
LOGICAl            :: InitialAutoRestartPartWeight= .FALSE.
INTEGER            :: InitialAutoRestartSample
LOGICAL            :: InterpolateSolution =.FALSE.
CHARACTER(LEN=300) :: RestartFile = ""
CHARACTER(LEN=255) :: NodeType_Restart
REAL               :: RestartTime
REAL               :: RestartWallTime ! wall time at the beginning of a simulation OR when a restart is performed via Load Balance
LOGICAL            :: RestartNullifySolution ! set the DG solution to zero (ignore the DG solution in the state file)
! Restart from sampled macroscopic variables
LOGICAL            :: DoMacroscopicRestart
CHARACTER(LEN=300) :: MacroRestartFileName
REAL, ALLOCATABLE  :: MacroRestartValues(:,:,:)
!===================================================================================================================================
END MODULE MOD_Restart_Vars
