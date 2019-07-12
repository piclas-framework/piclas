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
MODULE MOD_SurfaceModel_Analyze_Vars
!===================================================================================================================================
! Contains global variables used by the Analyze modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: SurfModelAnalyzeInitIsDone = .FALSE.
INTEGER                       :: SurfaceAnalyzeStep                  ! Analyze of surface is performed each Nth time step
LOGICAL                       :: CalcSurfNumSpec                     ! Calculate the number of simulated particles per species
                                                                     ! on surfaces
LOGICAL                       :: CalcEvaporation                     ! Calculate rate of evaporation [kg/s]
LOGICAL                       :: CalcSurfCoverage                    ! Calculate the surface coverages for each species
LOGICAL                       :: CalcAccomodation                    ! Calculate the surface accommodation coefficient
LOGICAL                       :: CalcAdsorbRates                     ! Calculate the adsorption probabilities of species
LOGICAL                       :: CalcSurfRates                       ! Calculate the surface reaction rate per reaction (k_r)
!===================================================================================================================================
END MODULE MOD_SurfaceModel_Analyze_Vars
