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
MODULE MOD_Gradient_Vars
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
! Dimension of the Variable of Interest
INTEGER                   :: Grad_DIM

! interior face values for all elements
REAL,ALLOCATABLE                      :: Var_master(:,:)
REAL,ALLOCATABLE                      :: Var_slave(:,:)
REAL,ALLOCATABLE                      :: Diff_side(:,:)
REAL,ALLOCATABLE                      :: Gradient_elem(:,:,:)

! Distances for Gradient Calculation
REAL,ALLOCATABLE                      :: Grad_dx_slave(:,:)
REAL,ALLOCATABLE                      :: Grad_dx_master(:,:)
REAL,ALLOCATABLE                      :: Grad_SysSol_slave(:,:)
REAL,ALLOCATABLE                      :: Grad_SysSol_master(:,:)
REAL,ALLOCATABLE                      :: Grad_SysSol_BC(:,:)

! Limiter
INTEGER                               :: GradLimiterType
REAL                                  :: GradLimVktK
!===================================================================================================================================
END MODULE MOD_Gradient_Vars