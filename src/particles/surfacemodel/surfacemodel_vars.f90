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
  INTEGER, ALLOCATABLE            :: Boundaries(:)
END TYPE

TYPE tSurfFluxMap
  TYPE(typeSurfaceflux),ALLOCATABLE  :: Surfaceflux(:)                   
END TYPE

  LOGICAL                         :: DoChemSurface 

TYPE tSurfReactions
  INTEGER                         :: NumOfReact             ! Number of possible reactions
  CHARACTER(LEN=5),ALLOCATABLE    :: ReactType(:)           ! Type of Reaction (reaction num)
                                                            !    A (adsorption)
                                                            !    D (desorption)
                                                            !    LH (Langmuir-Hinshlewood)
                                                            !    ER (Eley-Rideal)
  INTEGER, ALLOCATABLE            :: Reactants(:,:)         ! Reactants: indices of the species starting the reaction [NumOfReact,3]
  INTEGER, ALLOCATABLE            :: Products(:,:)          ! Products: indices of the species resulting from the reaction [NumOfReact,4]
  INTEGER, ALLOCATABLE            :: Inhibition(:)      ! Inhibition reaction
  INTEGER, ALLOCATABLE            :: NumOfBounds(:)         
  !REAL, ALLOCATABLE               :: ReactProb(:)
  REAL, ALLOCATABLE               :: EForm(:)
  ! Parameters for the adsorption
  REAL, ALLOCATABLE               :: S_initial(:)           ! Initial sticking coefficient
  REAL, ALLOCATABLE               :: MaxCoverage(:)         ! Maximal surface coverage
  REAL, ALLOCATABLE               :: DissOrder(:)           ! molecular or dissociative adsorption
  REAL, ALLOCATABLE               :: EqConstant(:)          ! adsorption/dissociation
  REAL, ALLOCATABLE               :: StickCoeff(:)         
  ! Parameters for the desorption
  REAL, ALLOCATABLE               :: E_initial(:)
  REAL, ALLOCATABLE               :: W_interact(:)
  REAL, ALLOCATABLE               :: C_a(:)
  REAL, ALLOCATABLE               :: C_b(:)
  ! General Parameters
  REAL, ALLOCATABLE               :: Rate(:)
  REAL, ALLOCATABLE               :: Prob(:)
  REAL, ALLOCATABLE               :: Prefactor(:)
  REAL, ALLOCATABLE               :: ArrheniusEnergy(:)
  LOGICAL, ALLOCATABLE            :: BoundisChemSurf(:)  
  TYPE(tBoundMap), ALLOCATABLE    :: BoundMap(:)        
  TYPE(tCollCaseInfo), ALLOCATABLE:: CollCaseInfo(:)        ! Information of collision cases (nCase) 
  TYPE(tSurfFluxMap), ALLOCATABLE :: SFMap(:)
END TYPE
TYPE(tSurfReactions)              :: SurfChemReac

TYPE typeSurfaceflux
  INTEGER                                :: BC                              
  CHARACTER(30)                          :: velocityDistribution           
  REAL                                   :: VeloIC             
  REAL                                   :: VeloVecIC(3)            
  REAL                                   :: MWTemperatureIC                  
  LOGICAL                                :: VeloIsNormal                    
  LOGICAL                                :: AcceptReject                    
  INTEGER                                :: ARM_DmaxSampleN               
  REAL                                   :: VFR_total                       
  REAL                     , ALLOCATABLE :: VFR_total_allProcs(:)          
  REAL                                   :: VFR_total_allProcsTotal         
  REAL                                   :: totalAreaSF                  
  INTEGER(KIND=8)                        :: InsertedParticle                 
  INTEGER(KIND=8)                        :: tmpInsertedParticle              
  INTEGER(KIND=8)                        :: tmpInsertedParticleSurplus      
  TYPE(tSurfFluxSubSideData), ALLOCATABLE :: SurfFluxSubSideData(:,:,:)     
  INTEGER                                :: dir(3)                          
  REAL                                   :: origin(2)                        
  REAL                                   :: rmax                           
  REAL                                   :: rmin                            
  LOGICAL                                :: Adaptive                        
  INTEGER                                :: AdaptiveType                
  REAL, ALLOCATABLE                      :: nVFRSub(:,:)                  
END TYPE

TYPE tSurfFluxSubSideData
  REAL                                   :: projFak                         
  REAL                                   :: a_nIn                           
  REAL                                   :: Velo_t1                       
  REAL                                   :: Velo_t2                         
  REAL                                   :: nVFR                            
  REAL                                   :: Dmax                            
END TYPE tSurfFluxSubSideData

! TYPE tRecombinationModel
! INTEGER                         :: NumOfReact             ! Number of possible reactions
! !CHARACTER(LEN=5),ALLOCATABLE    :: ReactType(:)           ! Type of Reaction (reaction num)
! !                                                            !    A (adsorption)
! !                                                            !    D (desorption)
! !                                                            !    LH (Langmuir-Hinshlewood)
! !                                                            !    ER (Eley-Rideal)
! INTEGER, ALLOCATABLE            :: Reactants(:,:)         ! Reactants: indices of the species starting the reaction [NumOfReact,3]
! INTEGER, ALLOCATABLE            :: Products(:,:)          ! Products: indices of the species resulting from the reaction [NumOfReact,4]
! INTEGER, ALLOCATABLE            :: NumOfBounds(:)         
! REAL, ALLOCATABLE               :: RecombCoeff(:)
! LOGICAL, ALLOCATABLE            :: BoundisCatSurf(:)  
! TYPE(tBoundMap), ALLOCATABLE    :: BoundMap(:)        
! TYPE(tCollCaseInfo), ALLOCATABLE:: CollCaseInfo(:)        ! Information of collision cases (nCase) 
! !TYPE(tSurfFluxMap), ALLOCATABLE :: SFMap(:)             
! END TYPE tRecombinationModel
! TYPE(tRecombinationModel),ALLOCATABLE   :: RecombModel(:) 

!===================================================================================================================================
END MODULE MOD_SurfaceModel_Vars
