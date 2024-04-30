!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
#include "piclas.h"

MODULE MOD_Particle_Emission_Vars
!===================================================================================================================================
!> Variables and types for the particle emission, used directly in MOD_Particle_Vars as types are part of the Species type
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE tInit                                                                   ! Particle Data for each init emission for each species
  !Specific Emission/Init values
  CHARACTER(40)                      :: SpaceIC                          ! specifying Keyword for Particle Space condition
  CHARACTER(30)                      :: velocityDistribution             ! specifying keyword for velocity distribution
  REAL                               :: Area                             ! Area for IC Rectangle
  REAL                               :: RadiusIC                         ! Radius for IC circle
  REAL                               :: Radius2IC                        ! Radius2 for IC cylinder (ring)
  REAL                               :: RadiusICGyro                     ! Radius for Gyrotron gyro radius
  REAL                               :: InflowRiseTime                   ! time to ramp the number of inflow particles
                                                                         ! linearly from zero to unity
  REAL                               :: NormalIC(3)                      ! Normal / Orientation of circle
  REAL                               :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                               :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                               :: NormalVector1IC(3)               ! 1st base vector normalized
  REAL                               :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                               :: NormalVector2IC(3)               ! 2nd base vector normalized
  REAL                               :: CuboidHeightIC                   ! third measure of cuboid
                                                                         ! (set 0 for flat rectangle),
                                                                         ! negative value = opposite direction
  REAL                               :: CylinderHeightIC                 ! third measure of cylinder
                                                                         ! (set 0 for flat rectangle),
                                                                         ! negative value = opposite direction
  REAL                               :: MinLocation(3)                   ! Minimal location for cell_local
  REAL                               :: MaxLocation(3)                   ! Maximal location for cell_local
  REAL                               :: VeloIC                           ! Velocity magnitude [m/s]
  REAL                               :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                               :: Amplitude                        ! Amplitude for sin-deviation initiation.
  REAL                               :: WaveNumber                       ! WaveNumber for sin-deviation initiation.
  INTEGER                            :: maxParticleNumberX               ! Maximum Number of all Particles in x direction
  INTEGER                            :: maxParticleNumberY               ! Maximum Number of all Particles in y direction
  INTEGER                            :: maxParticleNumberZ               ! Maximum Number of all Particles in z direction
  REAL                               :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                               :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                               :: PartDensity                      ! PartDensity (real particles per m^3)
  INTEGER                            :: ParticleEmissionType             ! Emission Type 0 = only initial,
                                                                         !               1 = emission rate in 1/s,
                                                                         !               2 = emission rate 1/iteration
  REAL                               :: ParticleNumber                   ! Initial, Emission in [1/s] or [1/Iteration]
  INTEGER(KIND=8)                    :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                    :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=4)                    :: InsertedParticleMisMatch=0       ! error in number of inserted particles of last step
#if USE_MPI
  INTEGER                            :: InitComm                         ! number of init-communicator
#endif /*USE_MPI*/
  INTEGER                            :: PartBCIndex                      ! Associated particle boundary ID
  REAL                               :: MacroParticleFactor              ! Emission-specific MPF
!=== photo ionization
  LOGICAL                            :: FirstQuadrantOnly  ! Only insert particles in the first quadrant that is spanned by the
                                                           ! vectors x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R]
  REAL                               :: PulseDuration      ! Pulse duration tau for a Gaussian-type pulse with
                                                           ! I~exp(-(t/tau)^2) [s]
  REAL                               :: WaistRadius        ! Beam waist radius (in focal spot) w_b for Gaussian-type pulse with
                                                           ! I~exp(-(r/w_b)^2) [m]
  REAL                               :: IntensityAmplitude ! Beam intensity maximum I0 Gaussian-type pulse with
                                                           ! I=I0*exp(-(t/tau)^2)exp(-(r/w_b)^2) [W/m^2]
  REAL                               :: WaveLength         ! Beam wavelength [m]
  REAL                               :: YieldSEE           ! Secondary photoelectron yield [-]
  REAL                               :: RepetitionRate     ! Pulse repetition rate [Hz]
  REAL                               :: Power              ! Average pulse power (energy of a single pulse times repetition rate) [J]
  REAL                               :: Energy             ! Single pulse energy (used when RepetitionRate and Power are not supplied [J]
  REAL                               :: Period             ! Time between the maximum intensity of two pulses [s]
  REAL                               :: tActive            ! Pulse will end at tActive (pulse time) [s]
  REAL                               :: tShift             ! Time shift for pulse corresponding to half of the Pulse width (pulse time) [s]
  INTEGER                            :: NbrOfPulses        ! Number of pulses [-]
  REAL                               :: NINT_Correction    ! nearest integer correction factor due to cut-off when converting
                                                           ! the number of particles calculated as real to integer for the
                                                           ! actual emission
  REAL                               :: WorkFunctionSEE    ! Photoelectron work function [eV]
  !REAL                               :: AngularBetaSEE
  REAL                               :: EffectiveIntensityFactor ! Scaling factor that increases I0 [-]
  INTEGER                            :: sumOfMatchedParticles    ! Sum of matched particles on all procs
  INTEGER                            :: sumOfRequestedParticles  ! Sum of requested particles on all procs
  INTEGER                            :: mySumOfMatchedParticles  ! Sum of matched particles on current proc
!=== Background gas regions
  INTEGER                            :: BGGRegion         ! Region number to be used for the species init
!=== Emission distribution
  CHARACTER(30)                      :: EmissionDistributionName  ! Species name, e.g., "electron" or "ArIon1" for particle emission
  REAL,ALLOCATABLE                   :: EmissionDistribution(:,:) !< pos (r,z or x,y,z) and particle properties (n, T, vx, vy, vz)
END TYPE tInit

! 2D Landmark
REAL, ALLOCATABLE :: PartPosLandmark(:,:)        ! Store particle positions during emission for placing
!                                                ! Electrons and ions at the exact same position
INTEGER           :: NbrOfParticleLandmarkMax    ! Array maximum size for storing positions
INTEGER           :: FractNbrOld,chunkSizeOld    ! Auxiliary integers for storing positions
LOGICAL              :: UseNeutralization           ! Flag for counting the charged particles impinging on a surface
CHARACTER(255)       :: NeutralizationSource        ! Name of the boundary for calculating the particle balance
INTEGER              :: nNeutralizationElems        ! Number of elements used for neutralization source (if required)
LOGICAL, ALLOCATABLE :: isNeutralizationElem(:)     ! Flag each element if it is a neutralization element
INTEGER, ALLOCATABLE :: NeutralizationBalanceElem(:)! Number of particles to be emitted within each neutralization element
INTEGER              :: NeutralizationBalance       ! Counter for charged particles (processor local): Add +1 for electrons and -1 for ions
INTEGER              :: NeutralizationBalanceGlobal ! Counter for charged particles (global): Add +1 for electrons and -1 for ions

! Bulk electron temperature
REAL              :: BulkElectronTemp            ! Bulk electron temperature for SEE model by Morozov2004
                                                 ! read-in in Kelvin (when using the SEE mode), but is directly converted
                                                 ! to eV for  usage in the code OR for neutralization BC (e.g. landmark)
LOGICAL           :: CalcBulkElectronTemp        ! Automatic bulk electron calculation
INTEGER           :: BulkElectronTempSpecID      ! Species ID (electron) for Automatic bulk electron calculation

! Emission distribution
LOGICAL              :: UseEmissionDistribution       !< Flag for activation particle emission by interpolation n, T and v (equidistant)
CHARACTER(255)       :: EmissionDistributionFileName  !< File name form which the data is read
INTEGER              :: EmissionDistributionN         !< Polynomial degree for particle emission in each element
INTEGER              :: EmissionDistributionDim       !< Spatial dimension of variable external field data: 1D, 2D or 3D
LOGICAL              :: EmissionDistributionAxisSym   !< True if the data is axis symmetric, e.g., B(r,z)
INTEGER              :: EmissionDistributionRadInd    !< Index of radial r-coordinate when using 2D data and axis symmetric
INTEGER              :: EmissionDistributionAxisDir   !< Direction that is used for the axial symmetric direction (1,2 or 3)
INTEGER              :: EmissionDistributionNum(1:3)  !< Number of points in x, y and z-direction
REAL                 :: EmissionDistributionMin(1:3)  !< Minimum values in x,y,z
REAL                 :: EmissionDistributionMax(1:3)  !< Maximum values in x,y,z
REAL                 :: EmissionDistributionDelta(1:3)!< equidistant z-spacing for the VariableExternalField (fast computation)
!===================================================================================================================================
END MODULE MOD_Particle_Emission_Vars
