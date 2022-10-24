#include "piclas.h"
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

USE MOD_DSMC_Vars,                ONLY:tCollCaseInfo

IMPLICIT NONE 
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER , ALLOCATABLE            :: SurfModResultSpec(:,:)          ! Resulting species after surface model treatment
                                                                    ! (nPartBound,nSpecies)
CHARACTER(LEN=50) , ALLOCATABLE  :: SurfModEnergyDistribution(:)    ! Energy distribution of the reflected particles
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

TYPE tBoundMap
  INTEGER, ALLOCATABLE                   :: Boundaries(:)          ! Map of the reactions to the boundaries
END TYPE

TYPE tPureSurf
  LOGICAL, ALLOCATABLE                   :: PureSurfReac(:)        ! List of boundaries on which LH/D reactions occur
END TYPE

LOGICAL                                  :: DoChemSurface          ! Call the surface catalysis routines

  CHARACTER(LEN=256)                     :: SpeciesDatabase        ! Name of the species database

TYPE tSurfReactions
  INTEGER                                :: NumOfReact             ! Number of catalytic reactions
  CHARACTER(LEN=64),ALLOCATABLE          :: CatName(:)
  LOGICAL                                :: OverwriteCatParameters ! Flag to read the cat parameters manually
  INTEGER                                :: SurfSpecies            ! Bulk species of the surface, involved in the reactions
  CHARACTER(LEN=255),ALLOCATABLE         :: ReactType(:)           ! Type of Reaction (reaction num)
                                                                   !    A (adsorption)
                                                                   !    D (desorption)
                                                                   !    LH (Langmuir-Hinshlewood)
                                                                   !    LHD (LH with instant desorption)
                                                                   !    ER (Eley-Rideal)
  INTEGER, ALLOCATABLE                   :: Reactants(:,:)         ! Reactants: indices of the species starting the reaction [NumOfReact,3]
  INTEGER, ALLOCATABLE                   :: Products(:,:)          ! Products: indices of the species resulting from the reaction [NumOfReact,4]
  INTEGER, ALLOCATABLE                   :: Inhibition(:)          ! Reaction number of inhibiting reactions
  INTEGER, ALLOCATABLE                   :: Promotion(:)           ! Reaction number of promoting reactions
  INTEGER, ALLOCATABLE                   :: NumOfBounds(:)         ! Number of catalytic boundaries
  ! Surface energy accomodation
  REAL, ALLOCATABLE                      :: EReact(:)              ! Reaction energy [K]
  REAL, ALLOCATABLE                      :: EScale(:)              ! Scaling factor for E_reac [K]
  REAL, ALLOCATABLE                      :: HeatAccomodation(:)    ! Beta coefficient for the energy accomodation
  ! Parameters for the adsorption
  REAL, ALLOCATABLE                      :: S_initial(:)           ! Initial sticking coefficient at zero coverage
  REAL, ALLOCATABLE                      :: MaxCoverage(:)         ! Maximal surface coverage
  REAL, ALLOCATABLE                      :: DissOrder(:)           ! Molecular (1) or dissociative (2) adsorption
  REAL, ALLOCATABLE                      :: EqConstant(:)          ! Equilibrium constant for adsorption/desorption
  REAL, ALLOCATABLE                      :: StickCoeff(:)          ! Sticking coefficient 
  ! Parameters for the desorption
  REAL, ALLOCATABLE                      :: E_initial(:)           ! Desorption energy at zero coverage [K]
  REAL, ALLOCATABLE                      :: W_interact(:)          ! Scaling factor for Edes [K]
  REAL, ALLOCATABLE                      :: C_a(:)                 ! Pre-exponential factor
  REAL, ALLOCATABLE                      :: C_b(:)                 ! Pre-exponential factor
  ! General Parameters
  REAL, ALLOCATABLE                      :: Rate(:)                ! Catalytic reaction rate [Cov/s*m^2]
  REAL, ALLOCATABLE                      :: Prob(:)                ! Catalytic reaction probability
  REAL, ALLOCATABLE                      :: Prefactor(:)           ! Pre-exponential factor [1/s]
  REAL, ALLOCATABLE                      :: ArrheniusEnergy(:)     ! Catalytic reaction energy [K]
  LOGICAL, ALLOCATABLE                   :: BoundisChemSurf(:)     ! Boundary with catalytic activity
  LOGICAL                                :: Diffusion              ! Activates instantaneous diffussion over the whole boundary
  LOGICAL                                :: TotDiffusion           ! Activates instantaneous diffussion over all boundaries
  INTEGER                                :: CatBoundNum            ! Number of catalytic boundaries
  TYPE(tBoundMap), ALLOCATABLE           :: BoundMap(:)            ! Map of the reactions to the boundaries
  TYPE(tPureSurf), ALLOCATABLE           :: PSMap(:)               ! Map for reactions occurring only on the surface   
  TYPE(tCollCaseInfo), ALLOCATABLE       :: CollCaseInfo(:)        ! Information of collision cases (nCase) 
  TYPE(tSurfaceflux), ALLOCATABLE        :: SurfaceFlux(:)         ! Surface flux data
END TYPE
TYPE(tSurfReactions)                     :: SurfChemReac

TYPE tSurfaceflux                                                           ! Surface flux properties on reactive boundaries
  INTEGER                                :: BC                              ! Catalytic boundary         
  CHARACTER(30)                          :: velocityDistribution            ! keyword for the velocity distribution        
  REAL                                   :: VeloIC                          ! velocity for inital Data
  REAL                                   :: VeloVecIC(3)                    ! normalized velocity vector
  REAL                                   :: MWTemperatureIC                 ! Temperature for Maxwell Distribution     
  LOGICAL                                :: VeloIsNormal                    ! VeloIC is in Surf-Normal instead of VeloVecIC    
  LOGICAL                                :: AcceptReject                    ! perform ARM for skewness of RefMap-positioning    
  INTEGER                                :: ARM_DmaxSampleN                 ! number of sample intervals in xi/eta for Dmax-calc.  
  REAL                                   :: VFR_total                       ! Total Volumetric flow rate through surface    
  REAL                     , ALLOCATABLE :: VFR_total_allProcs(:)           ! -''-, all values for root in ReduceNoise-case   
  REAL                                   :: VFR_total_allProcsTotal         !     -''-, total    
  REAL                                   :: totalAreaSF                     ! Total area of the respective surface flux 
  INTEGER(KIND=8)                        :: InsertedParticle                ! Number of all already inserted Particles     
  INTEGER(KIND=8)                        :: tmpInsertedParticle             ! tmp Number of all already inserted Particles     
  INTEGER(KIND=8)                        :: tmpInsertedParticleSurplus      ! tmp Number of all already inserted Particles    
  TYPE(tSurfFluxSubSideData), ALLOCATABLE :: SurfFluxSubSideData(:,:,:)     ! SF-specific Data of Sides (1:N,1:N,1:SideNumber)   
  INTEGER                                :: dir(3)                          ! axial (1) and orth. coordinates (2,3) of polar system    
  REAL                                   :: origin(2)                       ! origin in orth. coordinates of polar system     
  REAL                                   :: rmax                            ! max radius of to-be inserted particles   
  REAL                                   :: rmin                            ! min radius of to-be inserted particles    
  LOGICAL                                :: Adaptive                        ! Is the surface flux an adaptive boundary?    
  INTEGER                                :: AdaptiveType                    ! Chose the adaptive type, description in DefineParams
  REAL, ALLOCATABLE                      :: nVFRSub(:,:)                    ! normal volume flow rate through subsubside  
END TYPE

TYPE tSurfFluxSubSideData                                                   ! Reactive surface flux sub sides
  REAL                                   :: projFak                         ! VeloVecIC projected to inwards normal                           
  REAL                                   :: Velo_t1                         ! Velo comp. of first orth. vector
  REAL                                   :: Velo_t2                         ! Velo comp. of second orth. vector
  REAL                                   :: Dmax                            ! maximum Jacobian determinant of subside for opt. ARM              
  REAL,ALLOCATABLE                       :: nVFR(:)                         ! normal volume flow rate through subside     
  REAL,ALLOCATABLE                       :: a_nIn(:)                        ! speed ratio projected to inwards normal
END TYPE tSurfFluxSubSideData

REAL,ALLOCATABLE                         :: ChemSampWall(:,:,:,:,:)         ! Sampling direct impact mechanism
REAL,ALLOCATABLE                         :: ChemDesorpWall(:,:,:,:,:)       ! Desorption numbers
REAL,ALLOCPOINT                          :: ChemWallProp(:,:,:,:,:)         ! Adsorption count / heat flux
! INTEGER,ALLOCATABLE                      :: ChemCountReacWall(:,:,:,:,:)    ! Count the number of catalytic reactions on the subside

#if USE_MPI
INTEGER                                  :: ChemWallProp_Shared_Win         ! Adsorption count / heat flux
REAL,ALLOCPOINT                          :: ChemWallProp_Shared(:,:,:,:,:)    
REAL,POINTER                             :: ChemSampWall_Shared(:,:,:,:,:)  ! Sampling direct impact mechanism    
INTEGER                                  :: ChemSampWall_Shared_Win
#endif

!===================================================================================================================================
END MODULE MOD_SurfaceModel_Vars