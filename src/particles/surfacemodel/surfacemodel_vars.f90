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
LOGICAL                                  :: ModelERSpecular         ! Flag defining case for ER reflection (diffuse, specular)
LOGICAL                                  :: BlockingNeigh(3,3)      ! defines which of neighbour sites can block current site
                                                                    ! (nCurrentCoords,nNeighCoords) relevant for SMCR
#if (PP_TimeDiscMethod==42)
! definition of surfacmodel mean info container
TYPE tMeanInfo
  REAL                                   :: MeanProbAds             ! mean adsorption probability
  REAL                                   :: MeanProbDes             ! mean desorption probability
  INTEGER                                :: WallSpecNumCount        ! counter of Particles on Surface
  INTEGER                                :: NumOfAds                ! Number of Adsorptions on surfaces
  INTEGER                                :: NumOfDes                ! Number of Desorptions on surfaces
  REAL                                   :: Accomodation            ! Accomodation coeffcient calculated from Hard-Cube-Model
  INTEGER                                :: WallCollCount           ! counter of wallcollisions
END TYPE

! definition of proper info container allocated for each reflective surface
TYPE tProperInfo
  REAL     , ALLOCATABLE                 :: NumAdsReact(:)          ! mean probability for reaction at surface collision
  REAL     , ALLOCATABLE                 :: NumSurfReact(:)         ! mean probability for reaction on surface
  REAL     , ALLOCATABLE                 :: MeanAdsActE(:)          ! mean activation energy of adsorption reaction
  REAL     , ALLOCATABLE                 :: MeanAdsnu(:)            ! mean pre-exponential factor of desorption reaction
  REAL     , ALLOCATABLE                 :: MeanSurfActE(:)         ! mean activation energy of desorption reaction
  REAL     , ALLOCATABLE                 :: MeanSurfnu(:)           ! mean pre-exponential factor of desorption reaction
  REAL     , ALLOCATABLE                 :: ProperAdsActE(:)        ! activation energy of accepted adsorption reaction
  REAL     , ALLOCATABLE                 :: ProperAdsnu(:)          ! pre-exponential factor of accepted adsorption reaction
  REAL     , ALLOCATABLE                 :: ProperSurfActE(:)       ! activation energy of accepted desorption reaction
  REAL     , ALLOCATABLE                 :: ProperSurfnu(:)         ! pre-exponential factor accepted desorption reaction
  INTEGER  , ALLOCATABLE                 :: AdsReactCount(:)        ! Number of reactive adsorption probability calculations
  INTEGER  , ALLOCATABLE                 :: SurfReactCount(:)       ! Number of reactive desorption probability caluclations
  INTEGER  , ALLOCATABLE                 :: ProperAdsReactCount(:)  ! Number of reactive adsorptions
  INTEGER  , ALLOCATABLE                 :: ProperSurfReactCount(:) ! Number of reactive desorptions
  REAL     , ALLOCATABLE                 :: HeatFlux(:)             ! heatflux on surface due to species reacting on surface
                                                                    ! 1: adsorption process ; 2: desorption process
  REAL     , ALLOCATABLE                 :: HeatFluxDesCount(:)     ! heatflux on surface due to species reacting on surface
  REAL     , ALLOCATABLE                 :: HeatFluxAdsCount(:)     ! heatflux on surface due to species reacting on adsorption
END TYPE
#endif

TYPE tSurfaceModel
#if (PP_TimeDiscMethod==42)
  TYPE(tMeanInfo), ALLOCATABLE           :: Info(:)                 ! surfacemodel info for species n (nSpecies)
  TYPE(tProperInfo), ALLOCATABLE         :: ProperInfo(:)           ! proper surfacemodel info for species n (nSpecies)
#endif
  ! are all reset in updateSurfaceVars except SumEvapPart, which is reset in particle emission after particle inserting
  INTEGER , ALLOCATABLE                  :: SumDesorbPart(:,:,:,:)  ! Number of Particles of Species iSpec desorbed from Surface
                                                                    ! (nSurfSample,nSurfSample,nlocalSurfSide,nSpecies)
  INTEGER , ALLOCATABLE                  :: SumAdsorbPart(:,:,:,:)  ! Number of Particles of Species iSpec absorbed by the Surface
                                                                    ! (nSurfSample,nSurfSample,ntotalSurfSide,nSpecies)
  INTEGER , ALLOCATABLE                  :: SumReactPart(:,:,:,:)   ! Number of Particles that reacted on the surface
                                                                    ! (nSurfSample,nSurfSample,ntotalSurfSide,nSpecies)
  INTEGER , ALLOCATABLE                  :: SumERDesorbed(:,:,:,:)  ! Number of Particles desorbed through ER-type-reaction
                                                                    ! (nSurfSample,nSurfSample,ntotalSurfSide,nSpecies)
  INTEGER , ALLOCATABLE                  :: SumEvapPart(:,:,:,:)    ! number of evaporated particles
                                                                    ! (nSurfSample,nSurfSample,nlocalSurfSide,nSpecies)
END TYPE
TYPE(tSurfaceModel)                      :: SurfModel               ! SurfModel container

TYPE tAdsorption
#if (PP_TimeDiscMethod==42)
  LOGICAL                                :: LateralInactive         ! Flag for deactivation of lateral interactions in Q_a
  LOGICAL                                :: CoverageReduction       ! Flag for activating coverage reduction per iteration
  INTEGER, ALLOCATABLE                   :: CovReductionStep(:)     ! Step size for coverage reduction
  LOGICAL                                :: TPD                     ! Flag for TPD spectrum calculation
  REAL                                   :: TPD_beta                ! temperature slope for TPD [K/s]
  REAL                                   :: TPD_Temp                ! Walltemperature for TPD [K]
#endif
  INTEGER                                :: NumCovSamples           ! number of times coverage was sampled since last macrooutput
  LOGICAL , ALLOCATABLE                  :: SurfaceSpec(:,:)        ! set species as species of the surface
  REAL    , ALLOCATABLE                  :: Coverage(:,:,:,:)       ! coverage of surface i with species n
                                                                    ! (nSurfSample,nSurfSample,nSurfSide,nSpecies)
  REAL    , ALLOCATABLE                  :: ProbAds(:,:,:,:)        ! Adsorption probability of surface n
                                                                    ! (nSurfSample,nSurfSample,nSurfSide,nSpecies)
  REAL    , ALLOCATABLE                  :: ProbDes(:,:,:,:)        ! Desorption probability of surface n
                                                                    ! (nSurfSample,nSurfSample,nSurfSide,nSpecies)
  REAL    , ALLOCATABLE                  :: DensSurfAtoms(:)        ! density of surfaceatoms
  REAL    , ALLOCATABLE                  :: AreaIncrease(:)         ! Factor for increasing surface density
  INTEGER , ALLOCATABLE                  :: CrystalIndx(:)          ! Number of binding atoms in hollow site

  REAL    , ALLOCATABLE                  :: ReactCoeff(:,:)         ! Recaction coeffiecient (nPartBound,nSpecies)
  INTEGER , ALLOCATABLE                  :: ResultSpec(:,:)         ! Resulting species after surfacemodel treatment (nPartBound,nSpecies)
  REAL    , ALLOCATABLE                  :: ReactEnergy(:,:)        ! Reaction energy on surface (nPartBound,nSpecies)
  REAL    , ALLOCATABLE                  :: ReactAccomodation(:,:)  ! Energy Accomodation coefficient (nPartBound,nSpecies)
  ! parameters for Kisliuk and Polanyi Wigner model (surfacemodel=1)
  REAL    , ALLOCATABLE                  :: MaxCoverage(:,:)        ! maximum coverage of surface i with species n
  REAL    , ALLOCATABLE                  :: InitStick(:,:)          ! initial sticking coefficient (S_0) for surface n
  REAL    , ALLOCATABLE                  :: PrefactorStick(:,:)     ! prefactor of sticking coefficient for surface n
  INTEGER , ALLOCATABLE                  :: Adsorbexp(:,:)          ! Adsorption exponent for surface n
  REAL    , ALLOCATABLE                  :: Nu_a(:,:)               ! Nu exponent a for surface n
  REAL    , ALLOCATABLE                  :: Nu_b(:,:)               ! Nu exponent b for surface n
  REAL    , ALLOCATABLE                  :: DesorbEnergy(:,:)       ! Desorption energy (K) for surface n
  REAL    , ALLOCATABLE                  :: Intensification(:,:)    ! Intensification energy (K) for surface n
  ! parameters for UBI-QEP model (surfacemodel=3)
  LOGICAL                                :: EnableAdsAttraction
  LOGICAL                                :: NoDiffusion
  REAL    , ALLOCATABLE                  :: HeatOfAdsZero(:,:)      ! heat of adsorption (K) on clear surfaces
  INTEGER                                :: DissNum                 ! max number of dissociative surface reactions per species
  INTEGER                                :: RecombNum               ! max number of associative surface reactions per species
  INTEGER                                :: ReactNum                ! max number of diss/assoc surface reactions per species
  INTEGER , ALLOCATABLE                  :: DissocReact(:,:,:)      ! Resulting species for given dissoc (2,MaxDissNum,nSpecies)
  REAL    , ALLOCATABLE                  :: Diss_Prefactor(:,:)
  REAL    , ALLOCATABLE                  :: Diss_Powerfactor(:,:)
  REAL    , ALLOCATABLE                  :: ER_Prefactor(:,:)
  REAL    , ALLOCATABLE                  :: ER_Powerfactor(:,:)
  REAL    , ALLOCATABLE                  :: EDissBond(:,:)          ! Bond dissociation energy (K) for diss into resulting species
                                                                    ! (ReactNum,nspecies)
  REAL    , ALLOCATABLE                  :: EDissBondAdsorbPoly(:,:)! Bond dissociation energy (K) for diss into resulting species
                                                                    ! (ReactNum,nspecies)
  INTEGER , ALLOCATABLE                  :: RecombReact(:,:,:)      ! Partner/Result species for associative reaction (2,ReactNum,nSpecies)
  INTEGER , ALLOCATABLE                  :: ChemReactant(:,:)
  INTEGER , ALLOCATABLE                  :: ChemProduct(:,:)
  REAL    , ALLOCATABLE                  :: Reactant_DissBond_K(:,:)
  REAL    , ALLOCATABLE                  :: Product_DissBond_K(:,:)
  INTEGER                                :: nDissocReactions
  INTEGER                                :: nExchReactions
  INTEGER , ALLOCATABLE                  :: Coordination(:,:)       ! site bound coordination (1=hollow 2=bridge 3=on-top)(nSpecies)
  INTEGER , ALLOCATABLE                  :: DiCoord(:,:)            ! (1:nSpecies) bound via bridge bonding (=1) or chelating (=2)
  REAL    , ALLOCATABLE                  :: Ads_Powerfactor(:)
  REAL    , ALLOCATABLE                  :: Ads_Prefactor(:)
  INTEGER                                :: NumOfDissocReact        ! sum of all possible dissociation reactions on surfaces
  INTEGER                                :: NumOfRecombReact        ! sum of all possible recombination reactions on surfaces
  INTEGER                                :: NumOfExchReact          ! sum of all possible exchange reactions on surfaces
  ! TST Factor calculation variables
  LOGICAL , ALLOCATABLE                  :: TST_Calc(:,:)
  REAL    , ALLOCATABLE                  :: IncidentNormalTempAtSurf(:,:)
  REAL    , ALLOCATABLE                  :: SurfaceNormalVelo(:,:)
  REAL    , ALLOCATABLE                  :: SurfaceNormalVelo2(:,:)
  INTEGER , ALLOCATABLE                  :: CollSpecPartNum(:,:)
END TYPE
TYPE(tAdsorption)                        :: Adsorption              ! Adsorption-container

TYPE tSpeciesSurface
  REAL                                   :: ParamAntoine(1:3)                ! Parameter for Anointe Eq (vapor pressure)
  INTEGER                                :: condensCase
  REAL                                   :: liquidTkrit
  REAL                                   :: liquidTmelt
  REAL                                   :: liquidAlpha
  REAL                                   :: liquidBeta
  REAL                                   :: liquidBetaCoeff(1:6)
END TYPE
TYPE(tSpeciesSurface), ALLOCATABLE       :: SpecSurf(:)

TYPE tAdsorbateMapping
  INTEGER , ALLOCATABLE                  :: UsedSiteMap(:)          ! Mapping of adsorbateindex to surfposindex
                                                                    ! (1:SitesRemain) --> free site positions
                                                                    ! (SitesRemain+1:nSites) --> vacant site positions
  INTEGER , ALLOCATABLE                  :: Species(:)              ! species of adsorbate on sitepos (1:nSites)
  REAL    , ALLOCATABLE                  :: EVib(:)                 ! vibrational energy of adsorbate on sitepos (1:nSites)
  INTEGER , ALLOCATABLE                  :: BondAtomIndx(:,:)       ! adjacent surfatoms index x (1:nSites,1:nInterAtom)
  INTEGER , ALLOCATABLE                  :: BondAtomIndy(:,:)       ! adjacent surfatoms index y (1:nSites,1:nInterAtom)
  INTEGER                                :: nInterAtom              ! number of adjacent surface atoms
  INTEGER , ALLOCATABLE                  :: NeighPos(:,:)           ! pos of adjacent Neigbhour (1:nSites,1:nNeigbhours)
  INTEGER , ALLOCATABLE                  :: NeighSite(:,:)          ! site of adjacent Neigbhour (1:nSites,1:nNeigbhours)
  INTEGER                                :: nNeighbours             ! number of adjacent Neigbours sites
                                                                    ! (all possible Coordinations incl.)
  LOGICAL , ALLOCATABLE                  :: IsNearestNeigh(:,:)     ! Flag for defining nearest neighbour of binding site
#if USE_MPI
  LOGICAL , ALLOCATABLE                  :: Changed(:)              ! Flag if position changed during iteartion
#endif /*USE_MPI*/
END TYPE

TYPE tSurfaceDistributionInfo
  ! variables for surface distribution calculation
#if USE_MPI
  INTEGER , ALLOCATABLE                  :: Nbr_changed(:)
#endif /*USE_MPI*/
  INTEGER , ALLOCATABLE                  :: nSites(:)               ! max number of sites for site coordination (1:nCoordination=3)
  INTEGER , ALLOCATABLE                  :: SitesRemain(:)          ! number of empty sites for site coordination(1:nCoordination=3)
  INTEGER , ALLOCATABLE                  :: SurfAtomBondOrder(:,:,:)! bond order of surface atoms ((1:nSpecies,1:nXPos,1:nYPos)
                                                                    ! nXPos = nYPos = sqrt(nSites(3)) -> number of topsites
                                                                    ! applies for fcc(100) or similar surfaces
  REAL    , ALLOCATABLE                  :: desorbnum_tmp(:)
  REAL    , ALLOCATABLE                  :: adsorbnum_tmp(:)
  REAL    , ALLOCATABLE                  :: reactnum_tmp(:)
  TYPE(tAdsorbateMapping), ALLOCATABLE   :: AdsMap(:)               ! Mapping for adsorbates, adjacent surfatoms and neighbours
                                                                    ! (1:nCoordination)
END TYPE
TYPE(tSurfaceDistributionInfo),ALLOCATABLE   :: SurfDistInfo(:,:,:) ! Surface distribution tracking container

!===================================================================================================================================
END MODULE MOD_SurfaceModel_Vars
