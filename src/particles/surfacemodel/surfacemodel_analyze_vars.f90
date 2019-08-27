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
LOGICAL                       :: CalcAdsProb                         ! Calculate the probability of adsorption per species
LOGICAL                       :: CalcDesProb                         ! Calculate the probability of desorption per species
LOGICAL                       :: CalcSurfNumSpec                     ! Calculate the number of simulated particles per species
                                                                     ! on surfaces
#if (PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)
LOGICAL                       :: CalcSurfCoverage                    ! Calculate the surface coverages for each species
LOGICAL                       :: CalcAccomodation                    ! Calculate the surface accommodation coefficient
#if (PP_TimeDiscMethod==42)
LOGICAL                       :: CalcAdsorbRates                     ! Flag activating Analyze every refined rate data of gas-surface reactions
LOGICAL                       :: CalcAdsorbProb                      ! Flag activating Analyze eaction probabilities per reaction and species
LOGICAL                       :: CalcAdsorbnu                        ! Flag activating Analyze reaction frequencies (nu_r) per reaction and species'
LOGICAL                       :: CalcAdsorbE                         ! Flag activating Analyze activation barriers per reaction and species
LOGICAL                       :: CalcSurfRates                       ! Flag activating Analyze every refined rate data on the surfaces
LOGICAL                       :: CalcSurfProb                        ! Flag activating Analyze eaction probabilities per reaction and species
LOGICAL                       :: CalcSurfnu                          ! Flag activating Analyze reaction frequencies (nu_r) per reaction and species
LOGICAL                       :: CalcSurfE                           ! Flag activating Analyze activation barriers per reaction and species
LOGICAL                       :: CalcHeatflux                        ! Flag activating Analyze the the heat fluxes onto surface and corresponding
                                                                     ! reaction counters per reaction
#endif
#endif
!===================================================================================================================================
END MODULE MOD_SurfaceModel_Analyze_Vars
