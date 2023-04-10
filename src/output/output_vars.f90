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
INTEGER                    :: NVisu                    !< number of visualisation points is NVisu+1
CHARACTER(LEN=6),PARAMETER :: ProgramName='PICLas'
LOGICAL                    :: OutputInitIsDone=.FALSE.
LOGICAL                    :: doPrintStatusLine        !< flag indicating if status line should be printed
INTEGER                    :: userblock_len            !< length of userblock file in bytes
INTEGER                    :: userblock_total_len      !< length of userblock file + length of ini-file (with header) in bytes
CHARACTER(LEN=255)         :: UserBlockTmpFile='userblock.tmp' !< name of user block temp file
LOGICAL                    :: DoWriteStateToHDF5           !< only write HDF5 output if this is true
!===================================================================================================================================
END MODULE MOD_Output_Vars
