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
#include "piclas.h"

MODULE MOD_Particle_Emission_Init
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DefineParametersParticleEmission, InitializeVariablesSpeciesInits, InitialParticleInserting
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle initialization
!==================================================================================================================================
SUBROUTINE DefineParametersParticleEmission()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================

CALL prms%SetSection("Particle Initialization")

CALL prms%CreateIntOption(      'Part-Species[$]-nInits'  &
                                , 'Number of different initial particle placements for Species [$]', '0', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Reset'  &
                                , 'Flag for resetting species distribution with init during restart' &
                                , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-ChargeIC' &
                                , 'Particle charge of species [$], multiple of an elementary charge [C]' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MassIC'  &
                                , 'Atomic mass of species [$] [kg]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MacroParticleFactor' &
                                , 'Particle weighting factor: number of simulation particles per real particle for species [$]' &
                                , numberedmulti=.TRUE.)
#if defined(IMPA)
CALL prms%CreateLogicalOption(  'Part-Species[$]-IsImplicit'  &
                                , 'Flag if specific particle species is implicit', '.FALSE.', numberedmulti=.TRUE.)
#endif
CALL prms%CreateLogicalOption(  'Part-Species[$]-IsIMDSpecies' &
                                , 'TODO-DEFINE-PARAMETER', '.FALSE.', numberedmulti=.TRUE.)

CALL prms%SetSection("Particle Initialization (Inits)")

CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-SpaceIC' &
                                , 'Specifying Keyword for particle space condition of species [$] in case of multiple inits' &
                                , 'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-velocityDistribution'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Specifying keyword for velocity distribution', 'constant'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-InflowRiseTime' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Time to ramp the number of inflow particles linearly from zero to unity'&
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-RadiusIC'  , 'Radius for IC circle'                 , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-Radius2IC' , 'Radius2 for IC cylinder (ring)' , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-RadiusICGyro','Radius for Gyrotron gyro radius','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-NormalIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Normal / Orientation of circle', '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-BasePointIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Base point for IC cuboid and IC sphere ', '0. , 0. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-BaseVector1IC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'First base vector for IC cuboid', '1. , 0. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-BaseVector2IC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Second base vector for IC cuboid', '0. , 1. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-CuboidHeightIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Height of cuboid if SpaceIC = cuboid. (set 0 for flat rectangle)'//&
                                  ',negative value = opposite direction', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-CylinderHeightIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Third measure of cylinder  (set 0 for flat rectangle),'//&
                                  ' negative value = opposite direction', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-VeloIC'  &
                                , 'Velocity magnitude [m/s]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-VeloVecIC'  &
                                , 'Normalized velocity vector', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-Amplitude'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Amplitude for sin-deviation initiation.', '0.01', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-WaveNumber'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'WaveNumber for sin-deviation initiation', '2.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-maxParticleNumber-x'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Maximum Number of all Particles in x direction', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-maxParticleNumber-y'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Maximum Number of all Particles in y direction', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-maxParticleNumber-z'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Maximum Number of all Particles in z direction', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-Alpha' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'WaveNumber for sin-deviation initiation.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-MWTemperatureIC' &
                                , 'Temperature for Maxwell distribution [K]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-PartDensity' &
                                , 'Number density (real particles per m^3)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-ParticleEmissionType'  &
                                , 'Emission Type \n'//&
                                  '0 = only initial,\n'//&
                                  '1 = emission rate in 1/s,\n'//&
                                  '2 = emission rate 1/iteration\n', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-ParticleNumber' &
                                , 'Particle number, initial, in [1/s] or [1/Iteration]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-NumberOfExcludeRegions'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Number of different regions to be excluded', '0', numberedmulti=.TRUE.)

CALL prms%SetSection("Particle Species Init RegionExculdes")
! some inits or exluded in some regions
CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-ExcludeRegion[$]-SpaceIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Specified keyword for excluded particle space condition of'//&
                                  ' species[$] in case of multiple inits  ', 'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-ExcludeRegion[$]-RadiusIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Radius for excluded IC circle', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-ExcludeRegion[$]-Radius2IC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Radius2 for excluded IC cylinder (ring)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-ExcludeRegion[$]-NormalIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Normal orientation of excluded circle', '0. , 0. , 1.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-ExcludeRegion[$]-BasePointIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Base point for excluded IC cuboid and IC sphere', '0. , 0. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-ExcludeRegion[$]-BaseVector1IC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'First base vector for excluded IC cuboid', '1. , 0. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-ExcludeRegion[$]-BaseVector2IC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Second base vector for excluded IC cuboid', '0. , 1. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-ExcludeRegion[$]-CuboidHeightIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Height of excluded cuboid, if'//&
                                  ' Part-Species[$]-Init[$]-ExcludeRegion[$]-SpaceIC=cuboid (set 0 for flat rectangle),'//&
                                  ' negative value = opposite direction', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-ExcludeRegion[$]-CylinderHeightIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Height of excluded cylinder, if'//&
                                  ' Part-Species[$]-Init[$]-ExcludeRegion[$]-SpaceIC=cylinder (set 0 for flat circle),'//&
                                  'negative value = opposite direction ', '1.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-NeutralizationSource'  &
                                , 'Name of the boundary used for calculating the charged particle balance used for thruster'//&
                                  ' neutralization (no default).' ,numberedmulti=.TRUE.)
! ====================================== photoionization =================================================================
CALL prms%CreateLogicalOption('Part-Species[$]-Init[$]-FirstQuadrantOnly','Only insert particles in the first quadrant that is'//&
                              ' spanned by the vectors x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R]',  '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-PulseDuration',&
                           'Pulse duration tau for a Gaussian-tpye pulse with I~exp(-(t/tau)^2) [s]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-WaistRadius',&
                           'Beam waist radius (in focal spot) w_b for Gaussian-tpye pulse with I~exp(-(r/w_b)^2) [m]',&
                            numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-IntensityAmplitude',&
                           'Beam intensity maximum I0 Gaussian-tpye pulse with I=I0*exp(-(t/tau)^2)exp(-(r/w_b)^2) [W/m^2]','-1.0',&
                            numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-WaveLength','Beam wavelength [m]',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-YieldSEE','Secondary photoelectron yield [-]',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-RepetitionRate','Pulse repetition rate (pulses per second) [Hz]',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-Power','Average pulse power (energy of a single pulse times repetition rate) [W]',&
                           '-1.0',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-Energy','Single pulse energy [J]','-1.0',numberedmulti=.TRUE.)
CALL prms%CreateIntOption( 'Part-Species[$]-Init[$]-NbrOfPulses','Number of pulses [-]','1',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-WorkFunctionSEE','Photoelectron work function [eV]', numberedmulti=.TRUE.)
!CALL prms%CreateRealOption('Part-Species[$]-Init[$]-AngularBetaSEE',&
                           !'Orbital configuration of the solid from which the photoelectrons emerge','0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-EffectiveIntensityFactor', 'Scaling factor that increases I0 [-]',&
                            '1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Species[$]-Init[$]-TraceSpecies','Flag background species as trace element.'//&
                              ' Different weighting factor can be used',  '.FALSE.', numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersParticleEmission


SUBROUTINE InitializeVariablesSpeciesInits()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas
USE MOD_DSMC_BGGas             ,ONLY: BGGas_Initialize
#if USE_MPI
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iInit
CHARACTER(32)         :: hilf, hilf2
!===================================================================================================================================
BGGas%NumberOfSpecies = 0
ALLOCATE(BGGas%BackgroundSpecies(nSpecies))
BGGas%BackgroundSpecies = .FALSE.
ALLOCATE(BGGas%TraceSpecies(nSpecies))
BGGas%TraceSpecies = .FALSE.
ALLOCATE(BGGas%NumberDensity(nSpecies))
BGGas%NumberDensity = 0.
ALLOCATE(SpecReset(1:nSpecies))
SpecReset=.FALSE.
UseNeutralization = .FALSE.

DO iSpec = 1, nSpecies
  SWRITE (UNIT_stdOut,'(66(". "))')
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(hilf)//'-nInits')
#if USE_MPI
  IF(.NOT.PerformLoadBalance) THEN
#endif /*USE_MPI*/
    SpecReset(iSpec)                     = GETLOGICAL('Part-Species'//TRIM(hilf)//'-Reset')
#if USE_MPI
  END IF
#endif /*USE_MPI*/
  Species(iSpec)%ChargeIC              = GETREAL('Part-Species'//TRIM(hilf)//'-ChargeIC')
  Species(iSpec)%MassIC                = GETREAL('Part-Species'//TRIM(hilf)//'-MassIC')
  Species(iSpec)%MacroParticleFactor   = GETREAL('Part-Species'//TRIM(hilf)//'-MacroParticleFactor')
#if defined(IMPA)
  Species(iSpec)%IsImplicit            = GETLOGICAL('Part-Species'//TRIM(hilf)//'-IsImplicit')
#endif
  ALLOCATE(Species(iSpec)%Init(1:Species(iSpec)%NumberOfInits))
  DO iInit = 1, Species(iSpec)%NumberOfInits
    WRITE(UNIT=hilf2,FMT='(I0)') iInit
    hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
    ! Read-in of type and particle number for emission per iteration
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(hilf2)//'-ParticleEmissionType')
    Species(iSpec)%Init(iInit)%ParticleNumber        = GETREAL('Part-Species'//TRIM(hilf2)//'-ParticleNumber')
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity')
    Species(iSpec)%Init(iInit)%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC')
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-SpaceIC','cell_local'))
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution'&
      ,'maxwell_lpn'))
    Species(iSpec)%Init(iInit)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC')
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3)
    !--- Normalize VeloVecIC
    IF(.NOT.ALL(Species(iSpec)%Init(iInit)%VeloVecIC(:).EQ.0.)) THEN
      Species(iSpec)%Init(iInit)%VeloVecIC = Species(iSpec)%Init(iInit)%VeloVecIC / VECNORM(Species(iSpec)%Init(iInit)%VeloVecIC)
    END IF
    ! Additional read-in for circular/cuboid cases
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('disc','circle','circle_equidistant','gyrotron_circle','cylinder','sphere','photon_cylinder','photon_SEE_disc')
      Species(iSpec)%Init(iInit)%RadiusIC               = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusIC')
      Species(iSpec)%Init(iInit)%Radius2IC              = GETREAL('Part-Species'//TRIM(hilf2)//'-Radius2IC')
      Species(iSpec)%Init(iInit)%CylinderHeightIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-CylinderHeightIC')
    CASE('cuboid')
      Species(iSpec)%Init(iInit)%CuboidHeightIC         = GETREAL('Part-Species'//TRIM(hilf2)//'-CuboidHeightIC')
    END SELECT
    ! Space-ICs requiring basic geometry information and other options (all excluding cell_local and background)
    IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cell_local').AND. &
       (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'background')) THEN
      Species(iSpec)%Init(iInit)%NormalIC               = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-NormalIC',3)
      Species(iSpec)%Init(iInit)%BasePointIC            = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BasePointIC',3)
      Species(iSpec)%Init(iInit)%BaseVector1IC          = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector1IC',3)
      Species(iSpec)%Init(iInit)%BaseVector2IC          = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector2IC',3)
      !--- Normalize NormalIC (and BaseVector 1 & 3 IC for cylinder/sphere) for Inits
      Species(iSpec)%Init(iInit)%NormalIC = Species(iSpec)%Init(iInit)%NormalIC / VECNORM(Species(iSpec)%Init(iInit)%NormalIC)
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere')) THEN
        Species(iSpec)%Init(iInit)%BaseVector1IC = Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector1IC &
                                                    / VECNORM(Species(iSpec)%Init(iInit)%BaseVector1IC)
        Species(iSpec)%Init(iInit)%BaseVector2IC = Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector2IC &
                                                    / VECNORM(Species(iSpec)%Init(iInit)%BaseVector2IC)
      END IF
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere')) THEN
        Species(iSpec)%Init(iInit)%NormalIC = Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%NormalIC &
                                              / VECNORM(Species(iSpec)%Init(iInit)%NormalIC)
      END IF
      Species(iSpec)%Init(iInit)%InflowRiseTime         = GETREAL('Part-Species'//TRIM(hilf2)//'-InflowRiseTime')
      Species(iSpec)%Init(iInit)%NumberOfExcludeRegions = GETINT('Part-Species'//TRIM(hilf2)//'-NumberOfExcludeRegions')
    ELSE
      Species(iSpec)%Init(iInit)%InflowRiseTime         = 0.
      Species(iSpec)%Init(iInit)%NumberOfExcludeRegions = 0
    END IF
    ! Additional read-in for specific cases
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sin_deviation') THEN
      Species(iSpec)%Init(iInit)%Amplitude              = GETREAL('Part-Species'//TRIM(hilf2)//'-Amplitude')
      Species(iSpec)%Init(iInit)%WaveNumber             = GETREAL('Part-Species'//TRIM(hilf2)//'-WaveNumber')
      Species(iSpec)%Init(iInit)%maxParticleNumberX     = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-x')
      Species(iSpec)%Init(iInit)%maxParticleNumberY     = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-y')
      Species(iSpec)%Init(iInit)%maxParticleNumberZ     = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-z')
    END IF
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'gyrotron_circle') THEN
      Species(iSpec)%Init(iInit)%RadiusICGyro           = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusICGyro')
      Species(iSpec)%Init(iInit)%Alpha                  = GETREAL('Part-Species'//TRIM(hilf2)//'-Alpha')
    END IF
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'photon_SEE_energy')THEN
      Species(iSpec)%Init(iInit)%WorkFunctionSEE        = GETREAL('Part-Species'//TRIM(hilf2)//'-WorkFunctionSEE')
    END IF
    ! Photoionization in cylindrical volume and SEE based on photon impact on a surface
    IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_SEE_disc').OR. &
       (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_cylinder')) THEN
      CALL InitializeVariablesPhotoIonization(iSpec,iInit,hilf2)
    END IF
    !--- Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
    ! for low-temperature partially magnetized plasmas (2019)
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('2D_landmark','2D_landmark_copy')
      Species(iSpec)%Init(iInit)%ParticleEmissionType = 8
      Species(iSpec)%Init(iInit)%NINT_Correction      = 0.0
    CASE('2D_landmark_neutralization')
      Species(iSpec)%Init(iInit)%ParticleEmissionType = 9
      NeutralizationSource = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-NeutralizationSource'))
      NeutralizationBalance = 0
      UseNeutralization = .TRUE.
    END SELECT ! SpaceIC = '2D_landmark' or '2D_landmark_copy'
    !---Inserted Particles
    Species(iSpec)%Init(iInit)%InsertedParticle         = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus  = 0
    !----------- various checks/calculations after read-in of Species(iSpec)%Init(iInit)%-data ----------------------------------!
    !--- Check if initial ParticleInserting is really used
    IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0) THEN
      IF ( (Species(iSpec)%Init(iInit)%ParticleNumber.EQ.0) .AND. (Species(iSpec)%Init(iInit)%PartDensity.EQ.0.)) THEN
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled as neither ParticleNumber"
        SWRITE(*,*) "nor PartDensity detected for Species, Init ", iSpec, iInit
        Species(iSpec)%Init(iInit)%ParticleEmissionType = -1
      END IF
    END IF
    ! Check if number density and particle number have been defined
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.)) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleNumber.GT.0) THEN
        CALL abort(__STAMP__, &
          'Either ParticleNumber or PartDensity can be defined for selected parameters, not both!')
      END IF
    END IF
    ! 2D simulation/variable time step only with cell_local and/or surface flux
    IF((Symmetry%Order.EQ.2).OR.VarTimeStep%UseVariableTimeStep) THEN
      IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cell_local') THEN
        CALL abort(__STAMP__&
          ,'ERROR: Particle insertion/emission for 2D/axisymmetric or variable time step only possible with'//&
            'cell_local-SpaceIC and/or surface flux!')
      END IF
    END IF
    ! 1D Simulation with cuboid-SpaceIC
    IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').AND.Symmetry%Order.EQ.1) THEN
      IF(Species(iSpec)%Init(iInit)%BasePointIC(2).NE.-0.5.AND.Species(iSpec)%Init(iInit)%BasePointIC(3).NE.-0.5 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector1IC(2).NE.1 .AND.Species(iSpec)%Init(iInit)%BaseVector1IC(3).NE.0 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector2IC(1).NE.0 .AND.Species(iSpec)%Init(iInit)%BaseVector2IC(2).NE.0 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector1IC(1).NE.0 .AND.Species(iSpec)%Init(iInit)%BaseVector2IC(3).NE.1 ) THEN
        SWRITE(*,*) 'For 1D simulations with SpaceIC=cuboid, the vectors have to be defined in the following form:'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-BasePointIC=(/x,-0.5,-0.5/), with x as the basepoint in x direction'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-BaseVector1IC=(/0.,1.,0/)'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-BaseVector2IC=(/0.,0.,1/)'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-CuboidHeightIC is the extension of the insertion region and has to be positive'
        CALL abort(__STAMP__,'See above')
      END IF
    END IF
    !--- integer check for ParticleEmissionType 2
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         (ABS(Species(iSpec)%Init(iInit)%ParticleNumber-INT(Species(iSpec)%Init(iInit)%ParticleNumber,8)).GT.0.0)) THEN
      CALL abort(__STAMP__, &
        ' If ParticleEmissionType = 2 (parts per iteration), ParticleNumber has to be an integer number')
    END IF
    !--- ExcludeRegions
    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      CALL InitializeVariablesExcludeRegions(iSpec,iInit,hilf2)
    END IF
    !--- Background gas
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'background') THEN
      IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN
        BGGas%NumberOfSpecies = BGGas%NumberOfSpecies + 1
        BGGas%BackgroundSpecies(iSpec) = .TRUE.
        BGGas%NumberDensity(iSpec)     = Species(iSpec)%Init(iInit)%PartDensity
        BGGas%TraceSpecies(iSpec)      = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-TraceSpecies')
        Species(iSpec)%Init(iInit)%ParticleEmissionType = -1
      ELSE
        CALL abort(__STAMP__, &
          'ERROR: Only one background definition per species is allowed!')
      END IF
    END IF
    !--- InflowRise
    IF(Species(iSpec)%Init(iInit)%InflowRiseTime.GT.0.)THEN
      IF(.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow)  CALL CollectiveStop(__STAMP__, &
        ' Linearly ramping of inflow-number-of-particles is only possible with PoissonRounding or DoTimeDepInflow!')
    END IF
  END DO ! iInit
END DO ! iSpec
IF(nSpecies.GT.0)THEN
  SWRITE (UNIT_stdOut,'(66(". "))')
END IF ! nSpecies.GT.0

!-- reading BG Gas stuff
!   (moved here from dsmc_init for switching off the initial emission)
IF (useDSMC) THEN
  IF (BGGas%NumberOfSpecies.GT.0) THEN
    CALL BGGas_Initialize()
  ELSE
    DEALLOCATE(BGGas%NumberDensity)
  END IF ! BGGas%NumberOfSpecies.GT.0
END IF !useDSMC

END SUBROUTINE InitializeVariablesSpeciesInits


SUBROUTINE InitialParticleInserting()
!===================================================================================================================================
!> Insert initial particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Dielectric_Vars         ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
USE MOD_DSMC_Vars               ,ONLY: useDSMC, DSMC
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SetParticleTimeStep
USE MOD_Part_Pos_and_Velo       ,ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_SetInitElectronVelo
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Vars           ,ONLY: Species,nSpecies,PDM,PEM, usevMPF, SpecReset, VarTimeStep
USE MOD_Restart_Vars            ,ONLY: DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, NbrOfParticle,iInit,iPart,PositionNbr
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' Initial particle inserting... '

CALL UpdateNextFreePosition()

CALL DetermineInitialParticleNumber()

DO iSpec = 1,nSpecies
  IF (DSMC%DoAmbipolarDiff) THEN
    IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
  END IF
  IF (DoRestart .AND. .NOT.SpecReset(iSpec)) CYCLE
  DO iInit = 1, Species(iSpec)%NumberOfInits
    IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0) THEN
      IF(Species(iSpec)%Init(iInit)%ParticleNumber.GT.HUGE(1)) CALL abort(__STAMP__,&
        ' Integer of initial particle number larger than max integer size: ',HUGE(1))
      NbrOfParticle = INT(Species(iSpec)%Init(iInit)%ParticleNumber,4)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle position for species ',iSpec,' ... '
      CALL SetParticlePosition(iSpec,iInit,NbrOfParticle)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',iSpec,' ... '
      CALL SetParticleVelocity(iSpec,iInit,NbrOfParticle)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle charge and mass for species ',iSpec,' ... '
      CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)
      IF (usevMPF) CALL SetParticleMPF(iSpec,NbrOfParticle)
      IF (VarTimeStep%UseVariableTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
      IF (useDSMC) THEN
        IF (DSMC%DoAmbipolarDiff) CALL AD_SetInitElectronVelo(iSpec,iInit,NbrOfParticle)
        DO iPart = 1, NbrOfParticle
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .NE. 0) THEN
            PDM%PartInit(PositionNbr) = iInit
          ELSE
            CALL abort(__STAMP__,&
              'ERROR in InitialParticleInserting: No free particle index - maximum nbr of particles reached?')
          END IF
        END DO
      END IF
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      CALL UpdateNextFreePosition()
    END IF  ! Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0
  END DO    ! Species(iSpec)%NumberOfInits
END DO      ! nSpecies

!--- set last element to current element (needed when ParticlePush is not executed, e.g. "delay")
DO iPart = 1,PDM%ParticleVecLength
  PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
END DO

!--- Remove particles from dielectric regions if DielectricNoParticles=.TRUE.
IF(DoDielectric)THEN
  IF(DielectricNoParticles)THEN
    DO iPart = 1,PDM%ParticleVecLength
      ! Remove particles in dielectric elements
      IF(isDielectricElem(PEM%LocalElemID(iPart)))THEN
        PDM%ParticleInside(iPart) = .FALSE.
      END IF
    END DO
  END IF
END IF

SWRITE(UNIT_stdOut,'(A)') ' ...DONE '

END SUBROUTINE InitialParticleInserting


SUBROUTINE InitializeVariablesExcludeRegions(iSpec,iInit,hilf2)
!===================================================================================================================================
!> Exclude regions allow to skip particle initialization in a defined region
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars   ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                 :: iSpec, iInit
CHARACTER(32)           :: hilf2
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iEx
CHARACTER(32)           :: hilf3
!===================================================================================================================================
ASSOCIATE(SpecInit => Species(iSpec)%Init(iInit))
ALLOCATE(SpecInit%ExcludeRegion(1:SpecInit%NumberOfExcludeRegions))
IF (((TRIM(SpecInit%SpaceIC).EQ.'cuboid').OR.(TRIM(SpecInit%SpaceIC).EQ.'cylinder').OR.(TRIM(SpecInit%SpaceIC).EQ.'sphere'))) THEN
  DO iEx=1,SpecInit%NumberOfExcludeRegions
    WRITE(UNIT=hilf3,FMT='(I0)') iEx
    hilf3=TRIM(hilf2)//'-ExcludeRegion'//TRIM(hilf3)
    SpecInit%ExcludeRegion(iEx)%SpaceIC = TRIM(GETSTR('Part-Species'//TRIM(hilf3)//'-SpaceIC','cuboid'))
    SpecInit%ExcludeRegion(iEx)%RadiusIC = GETREAL('Part-Species'//TRIM(hilf3)//'-RadiusIC','1.')
    SpecInit%ExcludeRegion(iEx)%Radius2IC = GETREAL('Part-Species'//TRIM(hilf3)//'-Radius2IC','0.')
    SpecInit%ExcludeRegion(iEx)%NormalIC = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-NormalIC',3,'0. , 0. , 1.')
    SpecInit%ExcludeRegion(iEx)%BasePointIC = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BasePointIC',3,'0. , 0. , 0.')
    SpecInit%ExcludeRegion(iEx)%BaseVector1IC = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector1IC',3,'1. , 0. , 0.')
    SpecInit%ExcludeRegion(iEx)%BaseVector2IC = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector2IC',3,'0. , 1. , 0.')
    SpecInit%ExcludeRegion(iEx)%CuboidHeightIC = GETREAL('Part-Species'//TRIM(hilf3)//'-CuboidHeightIC','1.')
    SpecInit%ExcludeRegion(iEx)%CylinderHeightIC = GETREAL('Part-Species'//TRIM(hilf3)//'-CylinderHeightIC','1.')
    !--normalize and stuff
    IF((TRIM(SpecInit%ExcludeRegion(iEx)%SpaceIC).EQ.'cuboid') .OR. &
          ((((.NOT.ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%BaseVector1IC(1),1.) &
        .OR. .NOT.ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%BaseVector1IC(2),0.)) &
        .OR. .NOT.ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%BaseVector1IC(3),0.)) &
      .OR. ((.NOT.ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%BaseVector2IC(1),0.) &
        .OR. .NOT.ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%BaseVector2IC(2),1.)) &
        .OR. .NOT.ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%BaseVector2IC(3),0.))) &
      .AND. (((ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%NormalIC(1),0.)) &
        .AND. (ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%NormalIC(2),0.))) &
        .AND. (ALMOSTEQUAL(SpecInit%ExcludeRegion(iEx)%NormalIC(3),1.))))) THEN
      !-- cuboid; or BV are non-default and NormalIC is default: calc. NormalIC for ExcludeRegions from BV1/2
      !   (for def. BV and non-def. NormalIC; or all def. or non-def.: Use User-defined NormalIC when ExclRegion is cylinder)
      SpecInit%ExcludeRegion(iEx)%NormalIC(1) &
        = SpecInit%ExcludeRegion(iEx)%BaseVector1IC(2) * SpecInit%ExcludeRegion(iEx)%BaseVector2IC(3) &
        - SpecInit%ExcludeRegion(iEx)%BaseVector1IC(3) * SpecInit%ExcludeRegion(iEx)%BaseVector2IC(2)
      SpecInit%ExcludeRegion(iEx)%NormalIC(2) &
        = SpecInit%ExcludeRegion(iEx)%BaseVector1IC(3) * SpecInit%ExcludeRegion(iEx)%BaseVector2IC(1) &
        - SpecInit%ExcludeRegion(iEx)%BaseVector1IC(1) * SpecInit%ExcludeRegion(iEx)%BaseVector2IC(3)
      SpecInit%ExcludeRegion(iEx)%NormalIC(3) &
        = SpecInit%ExcludeRegion(iEx)%BaseVector1IC(1) * SpecInit%ExcludeRegion(iEx)%BaseVector2IC(2) &
        - SpecInit%ExcludeRegion(iEx)%BaseVector1IC(2) * SpecInit%ExcludeRegion(iEx)%BaseVector2IC(1)
    ELSE IF ( (TRIM(SpecInit%ExcludeRegion(iEx)%SpaceIC).NE.'cuboid').AND. &
              (TRIM(SpecInit%ExcludeRegion(iEx)%SpaceIC).NE.'cylinder') )THEN
      CALL abort(__STAMP__&
        ,'Error in ParticleInit, ExcludeRegions must be cuboid or cylinder!')
    END IF
    IF (DOTPRODUCT(SpecInit%ExcludeRegion(iEx)%NormalIC(1:3)).GT.0.) THEN
      SpecInit%ExcludeRegion(iEx)%NormalIC = SpecInit%ExcludeRegion(iEx)%NormalIC / VECNORM(SpecInit%ExcludeRegion(iEx)%NormalIC)
      SpecInit%ExcludeRegion(iEx)%ExcludeBV_lenghts(1) = VECNORM(SpecInit%ExcludeRegion(iEx)%BaseVector1IC)
      SpecInit%ExcludeRegion(iEx)%ExcludeBV_lenghts(2) = VECNORM(SpecInit%ExcludeRegion(iEx)%BaseVector2IC)
    ELSE
      CALL abort(__STAMP__&
        ,'Error in ParticleInit, NormalIC Vector must not be zero!')
    END IF
  END DO !iEx
ELSE
  CALL abort(__STAMP__&
    ,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid, sphere or cylinder!')
END IF
END ASSOCIATE

END SUBROUTINE InitializeVariablesExcludeRegions


SUBROUTINE InitializeVariablesPhotoIonization(iSpec,iInit,hilf2)
!===================================================================================================================================
!> Read-in and initialize variables for the case of photo-ionization in a volume (cylinder) and SEE at a surface (disc)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars    ,ONLY: PI
USE MOD_ReadInTools
USE MOD_Particle_Vars   ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                 :: iSpec, iInit
CHARACTER(32)           :: hilf2
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: factor
!===================================================================================================================================
Species(iSpec)%Init(iInit)%ParticleEmissionType = 7
! Check coordinate system of normal vector and two tangential vectors (they must form an orthogonal basis)
ASSOCIATE( v1 => UNITVECTOR(Species(iSpec)%Init(iInit)%NormalIC)      ,&
            v2 => UNITVECTOR(Species(iSpec)%Init(iInit)%BaseVector1IC) ,&
            v3 => UNITVECTOR(Species(iSpec)%Init(iInit)%BaseVector2IC))
  IF(DOT_PRODUCT(v1,v2).GT.1e-4) CALL abort(&
      __STAMP__&
      ,'NormalIC and BaseVector1IC are not perpendicular! Their dot product yields ',RealInfoOpt=DOT_PRODUCT(v1,v2))
  IF(DOT_PRODUCT(v1,v3).GT.1e-4) CALL abort(&
      __STAMP__&
      ,'NormalIC and BaseVector2IC are not perpendicular! Their dot product yields ',RealInfoOpt=DOT_PRODUCT(v1,v3))
  IF(DOT_PRODUCT(v2,v3).GT.1e-4) CALL abort(&
      __STAMP__&
      ,'BaseVector1IC and BaseVector2IC are not perpendicular! Their dot product yields ',RealInfoOpt=DOT_PRODUCT(v2,v3))
END ASSOCIATE
Species(iSpec)%Init(iInit)%FirstQuadrantOnly = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-FirstQuadrantOnly')
Species(iSpec)%Init(iInit)%PulseDuration      = GETREAL('Part-Species'//TRIM(hilf2)//'-PulseDuration')
Species(iSpec)%Init(iInit)%tShift = SQRT(8.0) * Species(iSpec)%Init(iInit)%PulseDuration
Species(iSpec)%Init(iInit)%WaistRadius        = GETREAL('Part-Species'//TRIM(hilf2)//'-WaistRadius')
Species(iSpec)%Init(iInit)%WaveLength         = GETREAL('Part-Species'//TRIM(hilf2)//'-WaveLength')
Species(iSpec)%Init(iInit)%NbrOfPulses        = GETINT('Part-Species'//TRIM(hilf2)//'-NbrOfPulses')
Species(iSpec)%Init(iInit)%NINT_Correction    = 0.0
Species(iSpec)%Init(iInit)%Power              = GETREAL('Part-Species'//TRIM(hilf2)//'-Power')
Species(iSpec)%Init(iInit)%Energy             = GETREAL('Part-Species'//TRIM(hilf2)//'-Energy')
Species(iSpec)%Init(iInit)%IntensityAmplitude = GETREAL('Part-Species'//TRIM(hilf2)//'-IntensityAmplitude')
! Set dummy value as it might not be read
Species(iSpec)%Init(iInit)%RepetitionRate = -1.0
IF(Species(iSpec)%Init(iInit)%Power.GT.0.0)THEN
  Species(iSpec)%Init(iInit)%RepetitionRate = GETREAL('Part-Species'//TRIM(hilf2)//'-RepetitionRate')
  Species(iSpec)%Init(iInit)%Period = 1./Species(iSpec)%Init(iInit)%RepetitionRate
  SWRITE(*,*) 'Photoionization in cylindrical volume: Selecting mode [RepetitionRate and Power]'

  Species(iSpec)%Init(iInit)%Energy = Species(iSpec)%Init(iInit)%Power / Species(iSpec)%Init(iInit)%RepetitionRate
  CALL PrintOption('Single pulse energy: Part-Species'//TRIM(hilf2)//'-Energy [J]','CALCUL.',&
                    RealOpt=Species(iSpec)%Init(iInit)%Energy)

  Species(iSpec)%Init(iInit)%IntensityAmplitude = Species(iSpec)%Init(iInit)%Energy / &
    (Species(iSpec)%Init(iInit)%WaistRadius**2 * Species(iSpec)%Init(iInit)%PulseDuration * PI**(3.0/2.0))

  CALL PrintOption('Intensity amplitude: I0 [W/m^2]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%IntensityAmplitude)
ELSEIF(Species(iSpec)%Init(iInit)%Energy.GT.0.0)THEN
  ! Check if more than one pulse is required
  IF(Species(iSpec)%Init(iInit)%NbrOfPulses.GT.1)THEN
    Species(iSpec)%Init(iInit)%RepetitionRate     = GETREAL('Part-Species'//TRIM(hilf2)//'-RepetitionRate')
    Species(iSpec)%Init(iInit)%Period = 1./Species(iSpec)%Init(iInit)%RepetitionRate
  ELSE
    Species(iSpec)%Init(iInit)%Period = 2.0 * Species(iSpec)%Init(iInit)%tShift
  END IF ! Species(iSpec)%Init(iInit)%NbrOfPulses
  SWRITE(*,*) 'Photoionization in cylindrical volume: Selecting mode [Energy]'

  Species(iSpec)%Init(iInit)%IntensityAmplitude = Species(iSpec)%Init(iInit)%Energy / &
    (Species(iSpec)%Init(iInit)%WaistRadius**2 * Species(iSpec)%Init(iInit)%PulseDuration * PI**(3.0/2.0))

  CALL PrintOption('Intensity amplitude: I0 [W/m^2]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%IntensityAmplitude)
ELSEIF(Species(iSpec)%Init(iInit)%IntensityAmplitude.GT.0.0)THEN
  ! Check if more than one pulse is required
  IF(Species(iSpec)%Init(iInit)%NbrOfPulses.GT.1)THEN
    Species(iSpec)%Init(iInit)%RepetitionRate     = GETREAL('Part-Species'//TRIM(hilf2)//'-RepetitionRate')
    Species(iSpec)%Init(iInit)%Period = 1./Species(iSpec)%Init(iInit)%RepetitionRate
  ELSE
    Species(iSpec)%Init(iInit)%Period = 2.0 * Species(iSpec)%Init(iInit)%tShift
  END IF ! Species(iSpec)%Init(iInit)%NbrOfPulses
  SWRITE(*,*) 'Photoionization in cylindrical volume: Selecting mode [IntensityAmplitude]'

  ! Calculate energy: E = I0*w_b**2*tau*PI**(3.0/2.0)
  Species(iSpec)%Init(iInit)%Energy = Species(iSpec)%Init(iInit)%IntensityAmplitude*Species(iSpec)%Init(iInit)%WaistRadius**2&
                                      *Species(iSpec)%Init(iInit)%PulseDuration*PI**(3.0/2.0)
  CALL PrintOption('Single pulse energy: Part-Species'//TRIM(hilf2)//'-Energy [J]','CALCUL.',&
                    RealOpt=Species(iSpec)%Init(iInit)%Energy)
ELSE
  CALL abort(&
    __STAMP__&
    ,'Photoionization in cylindrical volume: Supply either power P and repetition rate f, or energy E or intensity maximum I0!')
END IF ! use RepetitionRate and Power

! Sanity check: overlapping of pulses is not implemented (use multiple emissions for this)
IF(2.0*Species(iSpec)%Init(iInit)%tShift.GT.Species(iSpec)%Init(iInit)%Period) CALL abort(&
  __STAMP__&
  ,'Pulse length (2*tShift) is greater than the pulse period. This is not implemented!')

! Calculate the corrected intensity amplitude (due to temporal "-tShift to tShift" and spatial cut-off "0 to R")
factor = PI**(3.0/2.0) * Species(iSpec)%Init(iInit)%WaistRadius**2 * Species(iSpec)%Init(iInit)%PulseDuration * &
    (1.0-EXP(-Species(iSpec)%Init(iInit)%RadiusIC**2/(Species(iSpec)%Init(iInit)%WaistRadius**2)))*&
    ERF(Species(iSpec)%Init(iInit)%tShift/Species(iSpec)%Init(iInit)%PulseDuration)
Species(iSpec)%Init(iInit)%IntensityAmplitude = Species(iSpec)%Init(iInit)%Energy / factor
CALL PrintOption('Corrected Intensity amplitude: I0_corr [W/m^2]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%IntensityAmplitude)

CALL PrintOption('Pulse period (Time between maximum of two pulses) [s]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%Period)

CALL PrintOption('Temporal pulse width (pulse time 2x tShift) [s]','CALCUL.',RealOpt=2.0*Species(iSpec)%Init(iInit)%tShift)
Species(iSpec)%Init(iInit)%tActive = REAL(Species(iSpec)%Init(iInit)%NbrOfPulses - 1)*Species(iSpec)%Init(iInit)%Period &
                                        + 2.0*Species(iSpec)%Init(iInit)%tShift
CALL PrintOption('Pulse will end at tActive (pulse final time) [s]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%tActive)

IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_cylinder') THEN
  Species(iSpec)%Init(iInit)%EffectiveIntensityFactor = GETREAL('Part-Species'//TRIM(hilf2)//'-EffectiveIntensityFactor')
ELSE
  Species(iSpec)%Init(iInit)%YieldSEE           = GETREAL('Part-Species'//TRIM(hilf2)//'-YieldSEE')
END IF

END SUBROUTINE InitializeVariablesPhotoIonization


SUBROUTINE DetermineInitialParticleNumber()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Vars   ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: PI
USE MOD_DSMC_Vars           ,ONLY: RadialWeighting, DSMC
USE MOD_Particle_Mesh_Vars  ,ONLY: LocalVolume
USE MOD_Particle_Vars       ,ONLY: PDM,Species,nSpecies,SpecReset,Symmetry
USE MOD_ReadInTools
USE MOD_Restart_Vars        ,ONLY: DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iSpec, iInit
INTEGER(KIND=8)             :: insertParticles
REAL                        :: A_ins
!===================================================================================================================================

! Do sanity check of max. particle number compared to the number that is to be inserted for certain insertion types
insertParticles = 0
DO iSpec=1,nSpecies
  IF (DSMC%DoAmbipolarDiff) THEN
    IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
  END IF
  IF (DoRestart .AND. .NOT.SpecReset(iSpec)) CYCLE
  DO iInit = 1, Species(iSpec)%NumberOfInits
    ! Skip inits utilized for emission per iteration
    IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.0) CYCLE
    ! Calculate the ParticleNumber from PartDensity
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.)) THEN
      SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
      CASE('cuboid','sphere','cylinder')
        SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution))
        CASE('constant','maxwell','maxwell_lpn')
          !--- cross-product of the base-vectors to determine the base area
          A_ins = DOTPRODUCT(CROSS(Species(iSpec)%Init(iInit)%BaseVector1IC,Species(iSpec)%Init(iInit)%BaseVector2IC))
          IF (A_ins.GT.0.) THEN
            !---calculation of initial (macro)particle number
            IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
              Species(iSpec)%Init(iInit)%ParticleNumber = INT(Species(iSpec)%Init(iInit)%PartDensity &
                / Species(iSpec)%MacroParticleFactor * Species(iSpec)%Init(iInit)%CuboidHeightIC * SQRT(A_ins))
            ELSE IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
              Species(iSpec)%Init(iInit)%ParticleNumber = INT(Species(iSpec)%Init(iInit)%PartDensity &
                / Species(iSpec)%MacroParticleFactor * Species(iSpec)%Init(iInit)%CylinderHeightIC &
                * PI * (Species(iSpec)%Init(iInit)%RadiusIC**2-Species(iSpec)%Init(iInit)%Radius2IC**2))
            ELSE IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere') THEN
              Species(iSpec)%Init(iInit)%ParticleNumber = INT(Species(iSpec)%Init(iInit)%PartDensity &
                / Species(iSpec)%MacroParticleFactor * (Species(iSpec)%Init(iInit)%RadiusIC**3 * 4./3. * PI))
            END IF
          ELSE
            CALL abort(__STAMP__, &
              'BaseVectors are parallel or zero!')
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__, &
            'Given velocity distribution is not supported with the SpaceIC cuboid/sphere/cylinder!')
        END SELECT ! Species(iSpec)%Init(iInit)%SpaceIC
      CASE('cell_local')
        SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution))
        CASE('constant','maxwell','maxwell_lpn','taylorgreenvortex')
          IF (LocalVolume.GT.0.) THEN
            Species(iSpec)%Init(iInit)%ParticleNumber = NINT(Species(iSpec)%Init(iInit)%PartDensity &
                                                                / Species(iSpec)%MacroParticleFactor * LocalVolume)
          ELSE
            CALL abort(__STAMP__,'Local mesh volume is zero!')
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__,'Given velocity distribution is not supported with the SpaceIC cell_local!')
        END SELECT  ! Species(iSpec)%Init(iInit)%velocityDistribution
      CASE('background')
        ! do nothing
      CASE DEFAULT
        SWRITE(*,*) 'SpaceIC is: ', TRIM(Species(iSpec)%Init(iInit)%SpaceIC)
        CALL abort(__STAMP__,'ERROR: Unknown SpaceIC for species: ', iSpec)
      END SELECT    ! Species(iSpec)%Init(iInit)%SpaceIC
    END IF
    ! Sum-up the number of particles to be inserted
    IF(Symmetry%Order.LE.2) THEN
      ! The radial scaling of the weighting factor has to be considered
      IF(RadialWeighting%DoRadialWeighting) Species(iSpec)%Init(iInit)%ParticleNumber = &
                                  INT(Species(iSpec)%Init(iInit)%ParticleNumber * 2. / (RadialWeighting%PartScaleFactor),8)
    END IF
#if USE_MPI
    insertParticles = insertParticles + INT(REAL(Species(iSpec)%Init(iInit)%ParticleNumber)/PartMPI%nProcs,8)
#else
    insertParticles = insertParticles + INT(Species(iSpec)%Init(iInit)%ParticleNumber,8)
#endif
  END DO
END DO

IF (insertParticles.GT.PDM%maxParticleNumber) THEN
  IPWRITE(UNIT_stdOut,*)' Maximum particle number : ',PDM%maxParticleNumber
  IPWRITE(UNIT_stdOut,*)' To be inserted particles: ',INT(insertParticles,4)
  CALL abort(__STAMP__,&
    'Number of to be inserted particles per init-proc exceeds max. particle number! ')
END IF

END SUBROUTINE DetermineInitialParticleNumber

END MODULE MOD_Particle_Emission_Init
