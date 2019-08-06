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
LOGICAL                       :: CalcCollCounter                     ! Calculate the number of surface collision and number of
                                                                     ! adsorbed particles per species
LOGICAL                       :: CalcDesCounter                      ! Calculate the number of desorption particle per species
LOGICAL                       :: CalcSurfNumSpec                     ! Calculate the number of simulated particles per species
                                                                     ! on surfaces
#if (PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)
LOGICAL                       :: CalcEvaporation                     ! Calculate rate of evaporation [kg/s]
LOGICAL                       :: CalcSurfCoverage                    ! Calculate the surface coverages for each species
LOGICAL                       :: CalcAccomodation                    ! Calculate the surface accommodation coefficient
LOGICAL                       :: CalcAdsorbRates                     ! Calculate the adsorption probabilities of species (k_r)
LOGICAL                       :: CalcAdsorbProb                      ! Calculate the surface reaction rate per reaction (P_r)
LOGICAL                       :: CalcAdsorbnu                        ! Calculate the surface reaction rate per reaction (nu_r)
LOGICAL                       :: CalcAdsorbE                         ! Calculate the surface reaction rate per reaction (E_r)
LOGICAL                       :: CalcSurfRates                       ! Calculate the surface reaction rate per reaction (k_r)
LOGICAL                       :: CalcSurfProb                        ! Calculate the surface reaction rate per reaction (P_r)
LOGICAL                       :: CalcSurfnu                          ! Calculate the surface reaction rate per reaction (nu_r)
LOGICAL                       :: CalcSurfE                           ! Calculate the surface reaction rate per reaction (E_r)
LOGICAL                       :: CalcHeatflux                        ! Calculate the surface reaction rate per reaction (q)
#endif
!===================================================================================================================================
END MODULE MOD_SurfaceModel_Analyze_Vars
