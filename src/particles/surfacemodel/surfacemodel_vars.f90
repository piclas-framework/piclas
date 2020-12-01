!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
MODULE MOD_SurfaceModel_Vars
!===================================================================================================================================
!> Contains the SurfaceModel variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! definition of surfacmodel mean info container
TYPE tMeanInfo
  INTEGER                                :: WallCollCount           ! counter of wallcollisions
  INTEGER                                :: NumOfAds                ! Number of Adsorptions on surfaces
  INTEGER                                :: NumOfDes                ! Number of Desorptions on surfaces
END TYPE

TYPE tSurfaceModel
  TYPE(tMeanInfo), ALLOCATABLE           :: Info(:)                 ! surfacemodel info for species n (nSpecies)
END TYPE
TYPE(tSurfaceModel)                      :: SurfModel               ! SurfModel container

TYPE tAdsorption
  INTEGER , ALLOCATABLE                  :: ResultSpec(:,:)         ! Resulting species after surfacemodel treatment (nPartBound,nSpecies)
END TYPE
TYPE(tAdsorption)                        :: Adsorption              ! Adsorption-container



!===================================================================================================================================
END MODULE MOD_SurfaceModel_Vars
