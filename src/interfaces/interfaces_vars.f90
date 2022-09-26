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
MODULE MOD_Interfaces_Vars
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: InterfacesInitIsDone = .FALSE.
INTEGER,ALLOCATABLE :: InterfaceRiemann(:)            ! face identifier for switching between different Riemann solvers
!-----------------------------------------------------------------------------------------------------------------------------------
! GEOMETRY - PRE-DEFINED COORDINATES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                        :: GeometryIsSet=.FALSE.  ! use pre-defined coordinates, e.g., a gyrotron tube
REAL,ALLOCATABLE               :: Geometry(:,:)          ! Coordinates of the pre-defined geometry
INTEGER                        :: GeometryNPoints        ! Number of Coordinates (first array dimension of Geometry(:,:))
REAL,ALLOCATABLE               :: GeometryMin(:)         ! Minimum value of Geometry coordinates of each column
REAL,ALLOCATABLE               :: GeometryMax(:)         ! Maximum value of Geometry coordinates of each column
END MODULE MOD_Interfaces_Vars
