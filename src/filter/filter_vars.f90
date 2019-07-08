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
