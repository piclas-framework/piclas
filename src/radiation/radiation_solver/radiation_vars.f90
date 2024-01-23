#include "piclas.h"
!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
MODULE MOD_Radiation_Vars
!===================================================================================================================================
! Contains the radiation variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tRadiationInput               ! DSMC output
  REAL                  :: Ttrans(4)                        ! Temperature (Tx, Ty, Tz, Tt)
  REAL                  :: NumDens                          ! Particle density
  REAL                  :: Tvib                             ! Vibrational Temp
  REAL                  :: Trot                             ! Rotational Temp
  REAL                  :: Telec                            ! Electronic Temp
  REAL                  :: IonizationEn                     ! ionization energy [1/cm]
  REAL                  :: Radius
  REAL                  :: Starkex                          ! Exponent for Stark broadening
  INTEGER               :: NuclCharge                       ! nuclear charge (0: atom, 1: single-ionized atom, ...)
  CHARACTER(LEN=256)    :: RadiationSpectraFileName
  LOGICAL               :: DoRadiation                      ! Flag to consider/cycle species
END TYPE

TYPE(tRadiationInput), ALLOCATABLE :: RadiationInput(:)

REAL                    :: NumDensElectrons                 ! Electron Density
REAL                    :: TElectrons                       ! Electron Temperature

TYPE tRadiationParameter           ! Radiation Wavelength Parameter
  REAL                  :: MinWaveLen                       ! minimum spectral wavelength
  REAL                  :: MaxWaveLen                       ! maximum spectral wavelength
  INTEGER               :: WaveLenDiscr                     ! number of points in calculated spectrum
  INTEGER               :: WaveLenDiscrCoarse               ! number of points in calculated spectrum
  INTEGER               :: WaveLenDiscrOutput               ! number of points in calculated spectrum
  INTEGER               :: WaveLenReductionFactorOutput     ! number of points in calculated spectrum
  INTEGER               :: WaveLenReductionFactor
  REAL                  :: WaveLenIncr                      ! wavelength increments
  REAL                  :: WaveLenIncrCoarse
  REAL, ALLOCATABLE     :: WaveLen(:)                       ! wavelength array
  REAL, ALLOCATABLE     :: WaveLenCoarse(:)
END TYPE tRadiationParameter

TYPE(tRadiationParameter) :: RadiationParameter

TYPE tRadiationSwitches            ! Radiation types and mechanisms
  INTEGER               :: RadType                          ! 1: particle radiation, 2: black body radiation
  LOGICAL               :: ff                               ! Switch for free-free radiation
  LOGICAL               :: bf                               ! Switch for bound-free radiation
  LOGICAL               :: bb_at                            ! Switch for atomic line radiation
  LOGICAL               :: bb_mol                           ! Switch for molecular band radiation
  LOGICAL               :: MacroRadInput                    ! Switch for input of DSMC files
  LOGICAL               :: SortCellsY                       ! Sorts Cells in y-direction for manually created input files (e.g. Laux's test case)
  LOGICAL               :: UseElectronicExcitation          ! Switch for using electronic excitation energies (t) OR T_electron for atoms and sqrt(T_vib*T_electron) for molecules (f)
END TYPE tRadiationSwitches

TYPE(tRadiationSwitches) :: RadiationSwitches

TYPE tMolecBands          ! Radiation Wavelength Parameter
  REAL,ALLOCATABLE      :: MolTransLines(:,:)               ! (no transition lines, (vu, vl, sumre/fcf, Franck-Condon-Factor, SRe2))
END TYPE tMolecBands

TYPE tSpeciesRadiation
  TYPE(tMolecBands),ALLOCATABLE  :: Bands(:)
  REAL, ALLOCATABLE     :: Level(:,:)                       ! (length, (degeneracy ge, level energy, Energy level index , quant num of released elec, Gaunt fac))
  REAL, ALLOCATABLE     :: LinesReal(:,:)                   ! (length, (center wl, A_ul, stark HW))
  INTEGER, ALLOCATABLE  :: LinesInt(:,:)                    ! (length, (lower level, upper level, degeneracy lower level, degeneracy  upper level))
  REAL, ALLOCATABLE     :: NumDensExc(:)                    ! (nLevels)
  REAL                  :: PartFunc                         ! PartitionFunction
  INTEGER               :: nLevels
  INTEGER               :: nLines
  INTEGER               :: nBands                           !number of Bands
  CHARACTER(LEN=256)    :: RadiationSpectraFileName
  CHARACTER(LEN=256),ALLOCATABLE :: EnergyLevelName(:)
  CHARACTER(LEN=256),ALLOCATABLE :: BandName(:)
  INTEGER, ALLOCATABLE  :: BandProperties(:,:)              !(numBands)(up level, low level, Type Index, smf)
  REAL, ALLOCATABLE     :: EnergyLevelProperties(:,:)       !(numEnergyLevel)(degen, te(elev_mol/eterm), D0, we, wexe, weye, weze, be, alpha, mu, nuspin, altnat, de, betae, re, A, Lambda)
  INTEGER, ALLOCATABLE  :: NumMolecularTransitions(:)       ! number of transitions(numBands)
END TYPE tSpeciesRadiation

TYPE(tSpeciesRadiation), ALLOCATABLE     :: SpeciesRadiation(:)         ! (nSpec)

REAL, ALLOCATABLE       :: Radiation_Absorption_SpeciesWave(:,:)
!REAL, ALLOCATABLE       :: Radiation_NumDens
REAL, ALLOCPOINT       :: Radiation_ElemEnergy_Species(:,:,:)! (number of species, number of mesh elements, 2(Emission,Absorption))



REAL,ALLOCPOINT                 :: Radiation_Emission_spec(:,:)     ! (WaveLen(:), number of mesh elements)
REAL,ALLOCPOINT                 :: Radiation_Absorption_spec(:,:)     ! (WaveLen(:), number of mesh elements)
INTEGER(KIND = 2), ALLOCPOINT   :: Radiation_Absorption_SpecPercent(:,:,:)  ! 1:RadiationParameter%WaveLenDiscrCoarse ,1:nSpecies, 1:nGlobalElems, KIND=2? TODO

REAL,ALLOCPOINT                 :: MacroRadInputParameters(:,:,:)   ! DSMC Output file (iElem, iSpec, 5 (density, Tvib, Trot, Telec, Ttrans_mean))

#if USE_MPI
INTEGER                         :: MacroRadInputParameters_Shared_Win
REAL,ALLOCPOINT                 :: MacroRadInputParameters_Shared(:,:,:)
INTEGER                         :: Radiation_Emission_Spec_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Shared(:,:)
INTEGER                         :: Radiation_Absorption_Spec_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_Absorption_Spec_Shared(:)
INTEGER                         :: Radiation_Absorption_SpecPercent_Shared_Win
INTEGER(KIND = 2),ALLOCPOINT    :: Radiation_Absorption_SpecPercent_Shared(:)
INTEGER                         :: Radiation_ElemEnergy_Species_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_ElemEnergy_Species_Shared(:,:,:)
#endif
!===================================================================================================================================
END MODULE MOD_Radiation_Vars
