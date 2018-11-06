!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
!===================================================================================================================================
!> Contains global variables provided by the output routines
!===================================================================================================================================
MODULE MOD_OutputRPVisu_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: nSamples_out                !< number of visualisation points is NVisu+1
INTEGER                       :: nCoords                     !< number of visualisation points is NVisu+1
CHARACTER(LEN=255),ALLOCATABLE:: CoordNames(:)               !< including space and time coordinates
REAL,ALLOCATABLE              :: RPData_out(:,:,:)           !< output state 
REAL,ALLOCATABLE              :: RPDataTimeAvg_out(:,:)      !< time average state 
LOGICAL                       :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_OutputRPVisu_Vars
