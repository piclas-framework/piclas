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
MODULE MOD_FV_Vars
!===================================================================================================================================
! Contains global variables used by the FV modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! FV solution
REAL,ALLOCATABLE,TARGET               :: U_FV(:,:,:,:,:)
! FV time derivative
REAL,ALLOCATABLE                      :: Ut_FV(:,:,:,:,:)

! interior face values for all elements
REAL,ALLOCATABLE                      :: U_master_FV(:,:,:,:),U_slave_FV(:,:,:,:)
REAL,ALLOCATABLE                      :: Flux_Master_FV(:,:,:,:)
REAL,ALLOCATABLE                      :: Flux_Slave_FV(:,:,:,:)

LOGICAL                               :: FVInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_FV_Vars
