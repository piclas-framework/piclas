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
