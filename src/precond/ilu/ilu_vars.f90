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
