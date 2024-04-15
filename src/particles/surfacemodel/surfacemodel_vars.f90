#include "piclas.h"
!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
MODULE MOD_SurfaceModel_Vars
!===================================================================================================================================
!> Contains the SurfaceModel variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING

USE MOD_DSMC_Vars                   ,ONLY: tCollCaseInfo
USE MOD_Particle_SurfaceFlux_Vars   ,ONLY: tSurfaceFlux

IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER , ALLOCATABLE            :: SurfModResultSpec(:,:)          ! Resulting species after surface model treatment
                                                                    ! (nPartBound,nSpecies)
CHARACTER(LEN=50) , ALLOCATABLE  :: SurfModEnergyDistribution(:)    ! Energy distribution of the reflected/created particle(s)
REAL , ALLOCATABLE               :: SurfModEmissionEnergy(:)        ! Energy of emitted particle for surface emission model (only available for SurfaceModel=7)
REAL , ALLOCATABLE               :: SurfModEmissionYield(:)         ! Emission yield factor for surface emission model (only changable for SurfaceModel=7)
REAL                             :: BackupVeloABS                   ! Backup of velocity during double-ARMfor 2nd SEE
! === Porous BC ====================================================================================================================
INTEGER                          :: nPorousBC                       ! Number of porous BCs
TYPE tPorousBC
  INTEGER                        :: BC                              ! Number of the reflective BC to be used as a porous BC
  REAL                           :: Pressure                        ! Pressure at the BC [Pa], user-given
  CHARACTER(LEN=50)              :: Type
  REAL                           :: PumpingSpeed                    ! Given/calculated pumping speed [m3/s]
  REAL                           :: DeltaPumpingSpeedKp             ! Proportional factor for the pumping speed controller
  REAL                           :: DeltaPumpingSpeedKi             ! Integral factor for the pumping speed controller
  CHARACTER(LEN=50)              :: Region                          ! Form of the porous BC: 'circular'
  LOGICAL                        :: UsingRegion                     ! Use only a smaller region on the BC as a porous BC (e.g. pump)
  INTEGER                        :: dir(3)                          ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                           :: origin(2)                       ! origin in orth. coordinates of polar system
  REAL                           :: rmax                            ! max radius of to-be inserted particles
  REAL                           :: rmin                            ! min radius of to-be inserted particles
END TYPE
TYPE(tPorousBC), ALLOCATABLE     :: PorousBC(:)                     ! Container for the porous BC, allocated with nPorousBC
!=== Heterogenous Surface BC ========================================================================================================

TYPE tPureSurf
  LOGICAL, ALLOCATABLE                   :: PureSurfReac(:)        ! List of boundaries on which LH/D reactions occur
END TYPE

TYPE, EXTENDS(tCollCaseInfo) :: tEventProbInfo
  REAL, ALLOCATABLE                      :: ProdTransACC(:)        ! Reaction-specific thermal accommodation
END TYPE

LOGICAL                                  :: DoChemSurface          ! Call the surface catalysis routines

TYPE tSurfChemistry ! General surface chemistry parameter
  INTEGER                                :: NumOfReact             ! Number of catalytic reactions
  LOGICAL                                :: OverwriteCatParameters ! Flag to read the cat parameters manually
  INTEGER                                :: SurfSpecies            ! Bulk species of the surface, involved in the reactions
  LOGICAL                                :: Diffusion              ! Activates instantaneous diffussion over the whole boundary
  LOGICAL                                :: TotDiffusion           ! Activates instantaneous diffussion over all boundaries
  INTEGER                                :: CatBoundNum            ! Number of catalytic boundaries
  LOGICAL, ALLOCATABLE                   :: BoundisChemSurf(:)     ! Boundary with catalytic activity
  TYPE(tPureSurf), ALLOCATABLE           :: PSMap(:)               ! Map for reactions occurring only on the surface
  TYPE(tCollCaseInfo), ALLOCATABLE       :: CollCaseInfo(:)        ! Information of collision cases (nCase)
  TYPE(tSurfaceFlux), POINTER            :: SurfaceFlux(:)         ! Surface flux data (using the regular surface flux type)
  TYPE(tSFAux), ALLOCATABLE              :: SFAux(:)               ! Additional surface flux data, where variables differ from the regular surface flux type
  ! Event probability
  TYPE(tEventProbInfo), ALLOCATABLE      :: EventProbInfo(:)       ! Number of reaction paths and their probability per species
END TYPE
TYPE(tSurfChemistry)                     :: SurfChem

TYPE tSurfReactions
  CHARACTER(LEN=64)                      :: CatName
  CHARACTER(LEN=255)                     :: ReactType              ! Type of Reaction (reaction num)
                                                                   !    P (event probability)
                                                                   !    A (adsorption)
                                                                   !    D (desorption)
                                                                   !    LH (Langmuir-Hinshelwood)
                                                                   !    LHD (LH with instant desorption)
                                                                   !    ER (Eley-Rideal)
  INTEGER                                :: Reactants(2)           ! Reactants: indices of the species starting the reaction [NumOfReact,2]
  INTEGER                                :: Products(3)            ! Products: indices of the species resulting from the reaction [NumOfReact,3]
  INTEGER                                :: Inhibition             ! Reaction number of inhibiting reactions
  INTEGER                                :: Promotion              ! Reaction number of promoting reactions
  INTEGER                                :: NumOfBounds            ! Number of catalytic boundaries
  INTEGER, ALLOCATABLE                   :: Boundaries(:)          ! Map of the reactions to the boundaries
  ! Surface energy accommodation
  REAL                                   :: EReact                 ! Reaction energy [K]
  REAL                                   :: EScale                 ! Scaling factor for E_reac [K]
  REAL                                   :: HeatAccommodation      ! Beta coefficient for the energy accommodation
  ! Parameters for the adsorption
  REAL                                   :: S_initial              ! Initial sticking coefficient at zero coverage
  REAL                                   :: MaxCoverage            ! Maximal surface coverage
  REAL                                   :: DissOrder              ! Molecular (1) or dissociative (2) adsorption
  REAL                                   :: EqConstant             ! Equilibrium constant for adsorption/desorption
  REAL                                   :: StickCoeff             ! Sticking coefficient
  ! Parameter for the dissociative-adsorption
  LOGICAL                                :: DissociativeAds        ! Dissociative adsorption where the other molecule half is desorbed
  INTEGER                                :: AdsorbedProduct        ! Species ID of the particle that stays adsorbed on the surface
  INTEGER                                :: GasProduct             ! Species ID of the particle that remains in the gas phase
  ! Parameters for the desorption
  REAL                                   :: E_initial              ! Desorption energy at zero coverage [K]
  REAL                                   :: W_interact             ! Scaling factor for Edes [K]
  REAL                                   :: C_a                    ! Pre-exponential factor
  REAL                                   :: C_b                    ! Pre-exponential factor
  ! General Parameters
  REAL                                   :: Rate                   ! Catalytic reaction rate [Cov/s*m^2]
  REAL                                   :: Prob                   ! Catalytic reaction probability
  REAL                                   :: Prefactor              ! Pre-exponential factor [1/s]
  REAL                                   :: ArrheniusEnergy        ! Catalytic reaction energy [K]
END TYPE
TYPE(tSurfReactions), ALLOCATABLE        :: SurfChemReac(:)

TYPE tSFAux
  REAL, ALLOCATABLE                      :: a_nIn(:,:,:,:)       ! Speed ratio projected to inwards normal (additionally to regular surface flux variable due to missing species dependency)
END TYPE

REAL,ALLOCATABLE                         :: ChemDesorpWall(:,:,:,:,:)     ! Desorption numbers
REAL,ALLOCPOINT                          :: ChemWallProp(:,:,:,:,:)       ! Adsorption count / heat flux
REAL,ALLOCPOINT                          :: ChemSampWall(:,:,:,:,:)       ! Sampling direct impact mechanism

#if USE_MPI
REAL,ALLOCPOINT                          :: ChemWallProp_Shared(:,:,:,:,:)  ! Adsorption count / heat flux
INTEGER                                  :: ChemWallProp_Shared_Win
REAL,ALLOCPOINT                          :: ChemSampWall_Shared(:,:,:,:,:)  ! Sampling direct impact mechanism
INTEGER                                  :: ChemSampWall_Shared_Win
#endif

! === SEE BC ====================================================================================================================
REAL                             :: BulkElectronTempSEE             ! Bulk electron temperature for SEE model by Morozov2004
                                                                    ! read-in in Kelvin (when using the SEE mode), but is directly
                                                                    ! converted to eV for usage in the code
LOGICAL                          :: SurfModSEEelectronTempAutoamtic ! BulkElectronTempSEE = BulkElectronTemp, which is calculated
                                                                    ! automatically for the first species ID for electrons
REAL, ALLOCATABLE                :: SurfModSEEPowerFit(:,:)         ! Power-fit coefficients (1=a, 2=b) of the form: a*T(ev)^b

! === Sticking coefficient from simple models/interpolation
REAL, ALLOCATABLE                :: StickingCoefficientData(:,:)    ! Data for the model using non-bounce and condensation probability
                                                                    ! [:,1]: Maximum impact angle for model parameters
                                                                    ! [:,2]: Cut-off angle for non-bounce probability
                                                                    ! [:,3:4]: Temperature limits for linear interpolation of condensation probability
!===================================================================================================================================
END MODULE MOD_SurfaceModel_Vars