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
PUBLIC :: InitializeVariablesSpeciesBoundary
PUBLIC :: InitializeEmissionSpecificMPF
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

CALL prms%CreateIntOption(    'Part-Species[$]-nInits'  , 'Number of different initial particle placements for Species [$]'      , '0'      , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Species[$]-Reset'   , 'Flag for resetting species distribution with init during restart'     , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-ChargeIC', 'Particle charge of species [$], multiple of an elementary charge [C]' , '0.'     , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-MassIC'  , 'Atomic mass of species [$] [kg]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-MacroParticleFactor' ,'Particle weighting factor: number of simulation particles per real particle for species [$]' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-TimeStepFactor' ,'Species-specific time step factor, multiplied by the ManualTimeStep [$]' , '1.', numberedmulti=.TRUE.)
#if defined(IMPA)
CALL prms%CreateLogicalOption(  'Part-Species[$]-IsImplicit'  , 'Flag if specific particle species is implicit', '.FALSE.', numberedmulti=.TRUE.)
#endif
CALL prms%CreateLogicalOption(  'Part-Species[$]-IsIMDSpecies', 'Flag if particle species is used for IMD coupling', '.FALSE.', numberedmulti=.TRUE.)

CALL prms%SetSection("Particle Initialization (Inits)")

CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-SpaceIC'             , 'Keyword for particle space condition of species [$] in case of multiple inits' , 'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-velocityDistribution', 'Keyword for velocity distribution', 'constant', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-InflowRiseTime'      , 'Time to ramp the number of inflow particles linearly from zero to unity', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-RadiusIC'            , 'Radius for IC circle'                 , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-Radius2IC'           , 'Radius2 for IC cylinder (ring)' , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-RadiusICGyro'        ,'Radius for Gyrotron gyro radius','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-NormalIC'            ,'Normal / Orientation of circle', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-BasePointIC'         ,'Base point for IC cuboid and IC sphere ', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-BaseVector1IC'       ,'First base vector for IC cuboid', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-BaseVector2IC'       ,'Second base vector for IC cuboid', '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-CuboidHeightIC'      ,'Height of cuboid if SpaceIC = cuboid. (set 0 for flat rectangle),negative value = opposite direction', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-CylinderHeightIC'    ,'Third measure of cylinder  (set 0 for flat rectangle), negative value = opposite direction', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-VeloIC', 'Velocity magnitude [m/s]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-VeloVecIC', 'Normalized velocity vector', '0. , 0. , 0.', numberedmulti=.TRUE.)
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
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-BGG-Distribution-SpeciesIndex'  &
                                , 'Background gas with a distribution: Input the species index to use from the read-in '//&
                                  'distribution of the DSMCState file', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-BGG-Region'  &
                                , 'Number of the region in which the given conditions shall be applied to', numberedmulti=.TRUE.)

! === Cell local
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-MinimalLocation', 'Minimal location', '-999. , -999., -999.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-MaximalLocation', 'Maximal location', '999. , 999. , 999.', numberedmulti=.TRUE.)

! === Neutralization BC
CALL prms%CreateStringOption(   'Part-Species[$]-Init[$]-NeutralizationSource'  , 'Name of the boundary used for calculating the charged particle balance used for thruster neutralization (no default).' ,numberedmulti=.TRUE.)
! === Photoionization
CALL prms%CreateLogicalOption('Part-Species[$]-Init[$]-FirstQuadrantOnly','Only insert particles in the first quadrant that is'//&
                              ' spanned by the vectors x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R]',  '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-PulseDuration','Pulse duration tau for a Gaussian-type pulse with I~exp(-(t/tau)^2) [s]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-WaistRadius','Beam waist radius (in focal spot) w_b for Gaussian-type pulse with I~exp(-(r/w_b)^2) [m]',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-IntensityAmplitude','Beam intensity maximum I0 Gaussian-type pulse with I=I0*exp(-(t/tau)^2)exp(-(r/w_b)^2) [W/m^2]','-1.0',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-WaveLength','Beam wavelength [m]',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-YieldSEE','Secondary photoelectron yield [-]. Number of emitted electrons per incident photon',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-RepetitionRate','Pulse repetition rate (pulses per second) [Hz]',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-Power','Average pulse power (energy of a single pulse times repetition rate) [W]','-1.0',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-Energy','Single pulse energy [J]','-1.0',numberedmulti=.TRUE.)
CALL prms%CreateIntOption( 'Part-Species[$]-Init[$]-NbrOfPulses','Number of pulses [-]','1',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-WorkFunctionSEE','Photoelectron work function [eV]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-EffectiveIntensityFactor', 'Scaling factor that increases I0 [-]','1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Species[$]-Init[$]-TraceSpecies','Flag background species as trace element.'//&
                              ' Different weighting factor can be used',  '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption( 'Part-Species[$]-Init[$]-PartBCIndex','Associated particle boundary ID','-1',numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-Species[$]-Init[$]-MacroParticleFactor', 'Emission-specific particle weighting factor: number of simulation particles per real particle',numberedmulti=.TRUE.)
! === Emission distribution from an external input
CALL prms%CreateStringOption( 'Part-Species[$]-Init[$]-EmissionDistributionName' , 'Name of the species, e.g., "electron" or "ArIon" used for initial emission via interpolation of n, T and v from equidistant field data (no default).' ,numberedmulti=.TRUE.)
CALL prms%CreateStringOption( 'Part-EmissionDistributionFileName'                , 'H5 or CSV file containing the data for initial emission via interpolation of n, T and v from equidistant field data', 'none')
CALL prms%CreateIntOption(    'Part-EmissionDistributionN'                       , 'Polynomial degree for particle emission in each element. The default value is 2(N+1) with N being the polynomial degree of the solution.')
END SUBROUTINE DefineParametersParticleEmission


SUBROUTINE InitializeVariablesSpeciesInits()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_DSMC_Vars        ,ONLY: useDSMC, BGGas
USE MOD_DSMC_BGGas       ,ONLY: BGGas_Initialize
#if USE_MPI
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_MPI*/
USE MOD_Restart_Vars     ,ONLY: DoRestart
USE MOD_Analyze_Vars     ,ONLY: DoSurfModelAnalyze
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iInit
CHARACTER(32)         :: hilf, hilf2, DefStr
REAL                  :: MPFOld
!===================================================================================================================================
ALLOCATE(SpecReset(1:nSpecies))
SpecReset=.FALSE.
UseNeutralization = .FALSE.
UseEmissionDistribution = .FALSE.
! Species-specific time step
VarTimeStep%UseSpeciesSpecific = .FALSE.
! Background gas
BGGas%NumberOfSpecies = 0
ALLOCATE(BGGas%BackgroundSpecies(nSpecies))
BGGas%BackgroundSpecies = .FALSE.
ALLOCATE(BGGas%TraceSpecies(nSpecies))
BGGas%TraceSpecies = .FALSE.

BGGas%UseDistribution = GETLOGICAL('Particles-BGGas-UseDistribution')
BGGas%nRegions = GETINT('Particles-BGGas-nRegions')
IF(BGGas%UseDistribution.AND.(BGGas%nRegions.GT.0)) THEN
  CALL CollectiveStop(__STAMP__,'ERORR: Background gas can either be used with a distribution OR regions!')
ELSEIF (BGGas%UseDistribution) THEN
  ALLOCATE(BGGas%DistributionSpeciesIndex(nSpecies))
ELSEIF (BGGas%nRegions.GT.0) THEN
  BGGas%UseRegions = .TRUE.
ELSE
  ALLOCATE(BGGas%NumberDensity(nSpecies))
  BGGas%NumberDensity = 0.
END IF

DO iSpec = 1, nSpecies
  LBWRITE (UNIT_stdOut,'(66(". "))')
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
  IF((iSpec.GT.1).AND.UseDSMC.AND.(.NOT.UsevMPF))THEN
    IF(.NOT.ALMOSTEQUALRELATIVE(Species(iSpec)%MacroParticleFactor,MPFOld,1e-5)) CALL CollectiveStop(__STAMP__,&
        'Different MPFs only allowed when using Part-vMPF=T')
  END IF ! (iSpec.GT.1).AND.UseDSMC.AND.(.NOT.UsevMPF)
  MPFOld = Species(iSpec)%MacroParticleFactor
  ! Species-specific time step
  Species(iSpec)%TimeStepFactor              = GETREAL('Part-Species'//TRIM(hilf)//'-TimeStepFactor')
  IF(Species(iSpec)%TimeStepFactor.NE.1.) THEN
    VarTimeStep%UseSpeciesSpecific = .TRUE.
    IF(Species(iSpec)%TimeStepFactor.GT.1.) CALL CollectiveStop(__STAMP__,'ERROR: Species-specific time step only allows factors below 1!')
#if (USE_HDG) && !(PP_TimeDiscMethod==500) && !(PP_TimeDiscMethod==508) && !(PP_TimeDiscMethod==509)
    CALL CollectiveStop(__STAMP__,'ERROR: Species-specific time step is only implemented with Euler, Leapfrog & Boris-Leapfrog time discretization!')
#endif /*(USE_HDG)*/
#if !(USE_HDG) && !(PP_TimeDiscMethod==4)
    CALL CollectiveStop(__STAMP__,'ERROR: Species-specific time step is only implemented with HDG and/or DSMC')
#endif /*!(USE_HDG) && !(PP_TimeDiscMethod==4)*/
  END IF
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
    CASE('disc','circle','circle_equidistant','gyrotron_circle','cylinder','sphere','photon_cylinder','photon_SEE_disc',&
          'photon_honeycomb','photon_SEE_honeycomb')
      Species(iSpec)%Init(iInit)%RadiusIC               = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusIC')
      Species(iSpec)%Init(iInit)%Radius2IC              = GETREAL('Part-Species'//TRIM(hilf2)//'-Radius2IC')
      Species(iSpec)%Init(iInit)%CylinderHeightIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-CylinderHeightIC')
      IF(Species(iSpec)%Init(iInit)%Radius2IC.GE.Species(iSpec)%Init(iInit)%RadiusIC) CALL CollectiveStop(__STAMP__,&
          'For this emission type RadiusIC must be greater than Radius2IC!')
    CASE('cuboid','photon_rectangle')
      Species(iSpec)%Init(iInit)%CuboidHeightIC = GETREAL('Part-Species'//TRIM(hilf2)//'-CuboidHeightIC')
    CASE('cell_local')
      Species(iSpec)%Init(iInit)%MinLocation            = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-MinimalLocation',3)
      Species(iSpec)%Init(iInit)%MaxLocation            = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-MaximalLocation',3)
    END SELECT
    ! Space-ICs requiring basic geometry information and other options (all excluding cell_local and background)
    IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cell_local').AND. &
       (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'background')) THEN
      Species(iSpec)%Init(iInit)%NormalIC               = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-NormalIC',3)
      Species(iSpec)%Init(iInit)%BasePointIC            = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BasePointIC',3)
      !--- Get BaseVector1IC and normalize it
      IF(Symmetry%Order.GE.2) THEN
        Species(iSpec)%Init(iInit)%BaseVector1IC          = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector1IC',3)
        Species(iSpec)%Init(iInit)%NormalVector1IC        = UNITVECTOR(Species(iSpec)%Init(iInit)%BaseVector1IC)
      ELSE
        Species(iSpec)%Init(iInit)%BaseVector1IC          = (/0.,1.,0./)
        Species(iSpec)%Init(iInit)%NormalVector1IC        = (/0.,1.,0./)
      END IF
      !--- Get BaseVector2IC and normalize it
      IF(Symmetry%Order.GE.3) THEN
        Species(iSpec)%Init(iInit)%BaseVector2IC          = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector2IC',3)
        Species(iSpec)%Init(iInit)%NormalVector2IC        = UNITVECTOR(Species(iSpec)%Init(iInit)%BaseVector2IC)
      ELSE IF(Symmetry%Order.EQ.2.AND..NOT.Symmetry%Axisymmetric.AND.TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
        Species(iSpec)%Init(iInit)%BaseVector2IC          = (/0.,1.,0./)
        Species(iSpec)%Init(iInit)%NormalVector2IC        = (/0.,1.,0./)
      ELSE
        Species(iSpec)%Init(iInit)%BaseVector2IC          = (/0.,0.,1./)
        Species(iSpec)%Init(iInit)%NormalVector2IC        = (/0.,0.,1./)
      END IF
      !--- Normalize NormalIC (and BaseVector 1 & 3 IC for cylinder/sphere) for Inits
      Species(iSpec)%Init(iInit)%NormalIC = UNITVECTOR(Species(iSpec)%Init(iInit)%NormalIC)
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
    ELSE
      Species(iSpec)%Init(iInit)%InflowRiseTime         = 0.
    END IF
    ! Specifically for photon SEE on rectangle
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_SEE_rectangle')THEN
      ! Calculate the height using the base vectors. Use the length of the smaller one.
      Species(iSpec)%Init(iInit)%CuboidHeightIC = &
          MIN(VECNORM(Species(iSpec)%Init(iInit)%BaseVector1IC),VECNORM(Species(iSpec)%Init(iInit)%BaseVector2IC))
    END IF ! TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'gyrotron_circle'
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
    ! photon_cylinder, photon_SEE_disc, photon_honeycomb, photon_SEE_honeycomb
    IF(StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_'))THEN
      CALL InitializeVariablesPhotoIonization(iSpec,iInit,hilf2)
    END IF
    !--- Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
    ! for low-temperature partially magnetized plasmas (2019)
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('2D_landmark','2D_landmark_copy')
      Species(iSpec)%Init(iInit)%ParticleEmissionType = 8
      Species(iSpec)%Init(iInit)%NINT_Correction      = 0.0
    CASE('2D_landmark_neutralization','2D_Liu2010_neutralization','3D_Liu2010_neutralization','2D_Liu2010_neutralization_Szabo',&
         '3D_Liu2010_neutralization_Szabo')
      Species(iSpec)%Init(iInit)%ParticleEmissionType = 9
      NeutralizationSource = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-NeutralizationSource'))
      NeutralizationBalance = 0
      UseNeutralization = .TRUE.
      DoSurfModelAnalyze = .TRUE.
      IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'3D_Liu2010_neutralization').OR.&
         (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'3D_Liu2010_neutralization_Szabo'))THEN
        Species(iSpec)%Init(iInit)%FirstQuadrantOnly = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-FirstQuadrantOnly')
      END IF ! TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'3D_Liu2010_neutralization'
    END SELECT ! SpaceIC = '2D_landmark' or '2D_landmark_copy'
    !---Inserted Particles
    Species(iSpec)%Init(iInit)%InsertedParticle         = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus  = 0
    !----------- various checks/calculations after read-in of Species(iSpec)%Init(iInit)%-data ----------------------------------!
    !--- Check if initial ParticleInserting is really used
    IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0) THEN
      IF ( (Species(iSpec)%Init(iInit)%ParticleNumber.EQ.0) .AND. (Species(iSpec)%Init(iInit)%PartDensity.EQ.0.)) THEN
        LBWRITE(*,*) "WARNING: Initial ParticleInserting disabled as neither ParticleNumber"
        LBWRITE(*,*) "nor PartDensity detected for Species, Init ", iSpec, iInit
        Species(iSpec)%Init(iInit)%ParticleEmissionType = -1
      END IF
    END IF
    ! Check if number density and particle number have been defined
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.)) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleNumber.GT.0) THEN
        CALL CollectiveStop(__STAMP__, 'Either ParticleNumber or PartDensity can be defined for selected parameters, not both!')
      END IF
    END IF
    ! 2D simulation/variable time step only with cell_local and/or surface flux
    IF(UseVarTimeStep.OR.VarTimeStep%UseSpeciesSpecific) THEN
      IF(Species(iSpec)%Init(iInit)%ParticleEmissionType.GT.0) &
        CALL CollectiveStop(__STAMP__,'ERROR: Particle insertion/emission for variable time step '//&
            'only possible with initial particle insertion and/or surface flux!')
    END IF
    !--- integer check for ParticleEmissionType 2
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         (ABS(Species(iSpec)%Init(iInit)%ParticleNumber-INT(Species(iSpec)%Init(iInit)%ParticleNumber,8)).GT.0.0)) THEN
      CALL CollectiveStop(__STAMP__,' If ParticleEmissionType = 2 (parts per iteration), ParticleNumber has to be an integer number')
    END IF
    !--- Background gas
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'background') THEN
      IF(BGGas%BackgroundSpecies(iSpec)) THEN
        ! Only regions allows multiple background inits (additionally, avoid counting the same species multiple times)
        IF(.NOT.BGGas%UseRegions) CALL CollectiveStop(__STAMP__, 'ERROR: Only one background definition per species is allowed!')
      ELSE
        ! Count each species only once
        BGGas%NumberOfSpecies = BGGas%NumberOfSpecies + 1
      END IF
      BGGas%BackgroundSpecies(iSpec) = .TRUE.
      BGGas%TraceSpecies(iSpec)      = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-TraceSpecies')
      Species(iSpec)%Init(iInit)%ParticleEmissionType = -1
      ! Read-in the species index for background gas distribution
      IF(BGGas%UseDistribution) THEN
        BGGas%DistributionSpeciesIndex(iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-BGG-Distribution-SpeciesIndex')
      ELSE IF(BGGas%UseRegions) THEN
        Species(iSpec)%Init(iInit)%BGGRegion = GETINT('Part-Species'//TRIM(hilf2)//'-BGG-Region')
        IF(Species(iSpec)%Init(iInit)%BGGRegion.GT.BGGas%nRegions) THEN
          CALL CollectiveStop(__STAMP__, 'ERROR: Given background gas region number is greater than the defined number of regions!')
        END IF
      ELSE
        BGGas%NumberDensity(iSpec) = Species(iSpec)%Init(iInit)%PartDensity
      END IF
    END IF
    !--- InflowRise
    IF(Species(iSpec)%Init(iInit)%InflowRiseTime.GT.0.)THEN
      IF(.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow)  CALL CollectiveStop(__STAMP__, &
        ' Linearly ramping of inflow-number-of-particles is only possible with PoissonRounding or DoTimeDepInflow!')
    END IF
    !--- Emission distribution
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'EmissionDistribution')THEN
      UseEmissionDistribution                             = .TRUE.
      Species(iSpec)%Init(iInit)%EmissionDistributionName = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-EmissionDistributionName'))
      Species(iSpec)%Init(iInit)%ParticleEmissionType     = 0 ! initial particle insert
      Species(iSpec)%Init(iInit)%ParticleNumber           = 0  ! force 0
    END IF ! TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'EmissionDistribution'
  END DO ! iInit
END DO ! iSpec
IF(nSpecies.GT.0)THEN
  LBWRITE (UNIT_stdOut,'(66(". "))')
END IF ! nSpecies.GT.0

!-- reading BG Gas stuff
!   (moved here from dsmc_init for switching off the initial emission)
IF (useDSMC) THEN
  IF (BGGas%NumberOfSpecies.GT.0) THEN
    CALL BGGas_Initialize()
  ELSE
    IF(BGGas%UseDistribution) THEN
      DEALLOCATE(BGGas%DistributionSpeciesIndex)
    ELSE
      DEALLOCATE(BGGas%NumberDensity)
    END IF
  END IF ! BGGas%NumberOfSpecies.GT.0
ELSE
  IF((BGGas%NumberOfSpecies.GT.0).OR.BGGas%UseDistribution) CALL CollectiveStop(__STAMP__,'BGG requires UseDSMC=T')
END IF !useDSMC

!-- Read Emission Distribution stuff
IF(UseEmissionDistribution.AND.(.NOT.DoRestart)) THEN
  EmissionDistributionFileName = GETSTR('Part-EmissionDistributionFileName')
  WRITE(DefStr,'(i4)') 2*(PP_N+1)
  EmissionDistributionN   = GETINT('Part-EmissionDistributionN', DefStr)
  CALL ReadUseEmissionDistribution()
END IF

END SUBROUTINE InitializeVariablesSpeciesInits


SUBROUTINE InitialParticleInserting()
!===================================================================================================================================
!> Insert initial particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, Pi
USE MOD_TimeDisc_Vars           ,ONLY: ManualTimeStep, RKdtFrac
USE MOD_Mesh_Vars               ,ONLY: SideToElem
USE MOD_Dielectric_Vars         ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
USE MOD_DSMC_Vars               ,ONLY: useDSMC, DSMC, SpecDSMC, CollisMode
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SetParticleTimeStep
USE MOD_part_emission_tools     ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_Part_Pos_and_Velo       ,ONLY: SetParticlePosition,SetParticleVelocity,ParticleEmissionFromDistribution
USE MOD_Part_Pos_and_Velo       ,ONLY: ParticleEmissionCellLocal
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_SetInitElectronVelo
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition, IncreaseMaxParticleNumber, GetNextFreePosition
USE MOD_Restart_Vars            ,ONLY: DoRestart
USE MOD_Particle_Vars           ,ONLY: Species,nSpecies,PDM,PEM, usevMPF, SpecReset, UseVarTimeStep, VarTimeStep
USE MOD_Particle_Sampling_Vars  ,ONLY: UseAdaptiveBC, AdaptBCMacroVal, AdaptBCMapElemToSample, AdaptBCPartNumOut
USE MOD_Particle_Sampling_Adapt ,ONLY: AdaptiveBCSampling
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfFluxSideSize, SurfMeshSubSideData
USE MOD_DSMC_Init               ,ONLY: SetVarVibProb2Elems
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iSpec,NbrOfParticle,iInit,iPart,PositionNbr,iSF,iSide,ElemID,SampleElemID,currentBC,jSample,iSample,BCSideID
REAL                :: TimeStepOverWeight, v_thermal, dtVar
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INITIAL PARTICLE INSERTING...'

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
          ' Integer of initial particle number larger than max integer size: ',IntInfoOpt=HUGE(1))
      SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
      ! --------------------------------------------------------------------------------------------------
      ! Cell-local particle emission: every processors loops over its own elements
      CASE('cell_local')
        LBWRITE(UNIT_stdOut,'(A,I0,A)') ' Initial cell local particle emission for species ',iSpec,' ... '
        CALL ParticleEmissionCellLocal(iSpec,iInit,NbrOfParticle)
        ! TODO: MOVE EVERYTHING INTO THE EMISSION ROUTINE
        CALL SetParticleVelocity(iSpec,iInit,NbrOfParticle)
        ! SetParticleMPF cannot be performed anymore since the PartMPF was set based on the SplitThreshold
        ! SetParticleTimeStep is done during the emission as well
        ! SetParticleChargeAndMass is done as well
      ! --------------------------------------------------------------------------------------------------
      ! Cell-local particle emission from a given distribution
      CASE('EmissionDistribution')
        CALL ParticleEmissionFromDistribution(iSpec,iInit,NbrOfParticle)
        ! Particle velocity and species is set inside the above routine
        IF (usevMPF) CALL SetParticleMPF(iSpec,iInit,NbrOfParticle)
        IF (UseVarTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
      ! --------------------------------------------------------------------------------------------------
      ! Global particle emission
      CASE DEFAULT
        NbrOfParticle = INT(Species(iSpec)%Init(iInit)%ParticleNumber,4)
        LBWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle position for species ',iSpec,' ... '
        CALL SetParticlePosition(iSpec,iInit,NbrOfParticle)
        LBWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',iSpec,' ... '
        CALL SetParticleVelocity(iSpec,iInit,NbrOfParticle)
        IF (usevMPF) CALL SetParticleMPF(iSpec,iInit,NbrOfParticle)
        CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)
        IF (UseVarTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
      END SELECT
      IF (useDSMC) THEN
        IF (DSMC%DoAmbipolarDiff) CALL AD_SetInitElectronVelo(iSpec,iInit,NbrOfParticle)
        IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
          DO iPart = 1, NbrOfParticle
            PositionNbr = GetNextFreePosition(iPart)
            IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
              CALL DSMC_SetInternalEnr_Poly(iSpec,iInit,PositionNbr,1)
            ELSE
              CALL DSMC_SetInternalEnr_LauxVFD(iSpec,iInit,PositionNbr,1)
            END IF
          END DO
        END IF
      END IF
      ! Add new particles to particle vector length
      IF(NbrOfParticle.GT.0) PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(NbrOfParticle))
#ifdef CODE_ANALYZE
      IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
      DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
        IF (PDM%ParticleInside(iPart)) THEN
          IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
          CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
        END IF
      END DO
#endif
      ! Update
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

IF(UseAdaptiveBC.OR.(nPorousBC.GT.0)) THEN
!--- Sampling of near adaptive boundary element values after particle insertion to get initial distribution
  CALL AdaptiveBCSampling(initSampling_opt=.TRUE.)
  ! Adaptive BC Type = 4: Approximation of particles leaving the domain, assuming zero bulk velocity, using the fall-back values
  DO iSpec=1,nSpecies
    ! Species-specific time step
    IF(VarTimeStep%UseSpeciesSpecific) THEN
      dtVar = ManualTimeStep * RKdtFrac * Species(iSpec)%TimeStepFactor
    ELSE
      dtVar = ManualTimeStep * RKdtFrac
    END IF
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
      ! Skip processors without a surface flux
      IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) CYCLE
      ! Skip other regular surface flux and other types
      IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) CYCLE
      ! Calculate the velocity for the surface flux with the thermal velocity assuming a zero bulk velocity
      TimeStepOverWeight = dtVar / Species(iSpec)%MacroParticleFactor
      v_thermal = SQRT(2.*BoltzmannConst*Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) / (2.0*SQRT(PI))
      ! Loop over sides on the surface flux
      DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
        BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        SampleElemID = AdaptBCMapElemToSample(ElemID)
        IF(SampleElemID.GT.0) THEN
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            AdaptBCPartNumOut(iSpec,iSF) = AdaptBCPartNumOut(iSpec,iSF) + INT(AdaptBCMacroVal(4,SampleElemID,iSpec) &
              * TimeStepOverWeight * SurfMeshSubSideData(iSample,jSample,BCSideID)%area * v_thermal)
          END DO; END DO
        END IF  ! SampleElemID.GT.0
      END DO    ! iSide=1,BCdata_auxSF(currentBC)%SideNumber
    END DO      ! iSF=1,Species(iSpec)%nSurfacefluxBCs
  END DO        ! iSpec=1,nSpecies
END IF

IF((DSMC%VibRelaxProb.EQ.2).AND.(CollisMode.GE.2)) CALL SetVarVibProb2Elems()

LBWRITE(UNIT_stdOut,'(A)') ' INITIAL PARTICLE INSERTING DONE!'

END SUBROUTINE InitialParticleInserting


SUBROUTINE InitializeVariablesPhotoIonization(iSpec,iInit,hilf2)
!===================================================================================================================================
!> Read-in and initialize variables for the case of photo-ionization in a volume (cylinder) and SEE at a surface (disc)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars    ,ONLY: PI
USE MOD_ReadInTools
USE MOD_Particle_Vars   ,ONLY: Species
USE MOD_DSMC_Vars       ,ONLY: BGGas
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
! Abort if a background gas distribution is used (CalcPhotoIonizationNumber assumes a constant distribution)
IF(BGGas%UseDistribution) CALL abort(__STAMP__,&
  'ERROR: Photo-ionization and a background gas distribution is not implemented yet!')
! Check coordinate system of normal vector and two tangential vectors (they must form an orthogonal basis)
ASSOCIATE( n1 => UNITVECTOR(Species(iSpec)%Init(iInit)%NormalIC)      ,&
           n2 => UNITVECTOR(Species(iSpec)%Init(iInit)%BaseVector1IC) ,&
           n3 => UNITVECTOR(Species(iSpec)%Init(iInit)%BaseVector2IC) ,&
           v2 => Species(iSpec)%Init(iInit)%BaseVector1IC             ,&
           v3 => Species(iSpec)%Init(iInit)%BaseVector2IC             )
  !IF(DOT_PRODUCT(n1,n2).GT.1e-4) CALL abort(__STAMP__&
      !,TRIM(hilf2)//': NormalIC and BaseVector1IC are not perpendicular! Their dot product yields ',RealInfoOpt=DOT_PRODUCT(n1,n2))
  !IF(DOT_PRODUCT(n1,n3).GT.1e-4) CALL abort(__STAMP__&
      !,TRIM(hilf2)//': NormalIC and BaseVector2IC are not perpendicular! Their dot product yields ',RealInfoOpt=DOT_PRODUCT(n1,n3))
  IF(DOT_PRODUCT(n2,n3).GT.1e-4) CALL abort(__STAMP__&
      ,TRIM(hilf2)//': BaseVector1IC and BaseVector2IC are not perpendicular! Their dot product yields ',RealInfoOpt=DOT_PRODUCT(n2,n3))
  ! Settings only for rectangle emission
  SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
  CASE('photon_SEE_rectangle','photon_rectangle')
    ! Area is calculated from the cross product of the two base vectors
    Species(iSpec)%Init(iInit)%Area = VECNORM(CROSS(v2,v3))
    CALL PrintOption('Rectangular emission area for Part-Species'//TRIM(hilf2)//': A [m2]','CALCUL.',&
                    RealOpt=Species(iSpec)%Init(iInit)%Area)
  END SELECT
END ASSOCIATE
Species(iSpec)%Init(iInit)%FirstQuadrantOnly = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-FirstQuadrantOnly')
Species(iSpec)%Init(iInit)%PulseDuration      = GETREAL('Part-Species'//TRIM(hilf2)//'-PulseDuration')
Species(iSpec)%Init(iInit)%tShift = SQRT(8.0) * Species(iSpec)%Init(iInit)%PulseDuration
! Honeycomb and rectangle do not require a waist radius
SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
CASE('photon_SEE_honeycomb','photon_honeycomb','photon_SEE_rectangle','photon_rectangle')
  Species(iSpec)%Init(iInit)%WaistRadius        = GETREAL('Part-Species'//TRIM(hilf2)//'-WaistRadius','-1.0')
  IF(Species(iSpec)%Init(iInit)%FirstQuadrantOnly) CALL abort(__STAMP__,'FirstQuadrantOnly=T is not implemented for honeycombs')
CASE DEFAULT
  Species(iSpec)%Init(iInit)%WaistRadius        = GETREAL('Part-Species'//TRIM(hilf2)//'-WaistRadius')
END SELECT
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
  LBWRITE(*,*) 'Photoionization in cylindrical/rectangular/honeycomb volume: Selecting mode [RepetitionRate and Power]'

  Species(iSpec)%Init(iInit)%Energy = Species(iSpec)%Init(iInit)%Power / Species(iSpec)%Init(iInit)%RepetitionRate
  CALL PrintOption('Single pulse energy: Part-Species'//TRIM(hilf2)//'-Energy [J]','CALCUL.',&
                    RealOpt=Species(iSpec)%Init(iInit)%Energy)

  ASSOCIATE( E0   => Species(iSpec)%Init(iInit)%Energy             ,&
             wb   => Species(iSpec)%Init(iInit)%WaistRadius        ,&
             tau  => Species(iSpec)%Init(iInit)%PulseDuration      ,&
             Rout => Species(iSpec)%Init(iInit)%RadiusIC           ,&
             Rin  => Species(iSpec)%Init(iInit)%Radius2IC          ,&
             A    => Species(iSpec)%Init(iInit)%Area               )
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('photon_SEE_rectangle','photon_rectangle')
      Species(iSpec)%Init(iInit)%IntensityAmplitude = E0 / (SQRT(PI)*tau*A)
    CASE('photon_SEE_honeycomb','photon_honeycomb')
      Species(iSpec)%Init(iInit)%IntensityAmplitude = E0 / (SQRT(PI)*tau*(1.5*SQRT(3.0))*(Rout**2-Rin**2))
    CASE DEFAULT
      Species(iSpec)%Init(iInit)%IntensityAmplitude = E0 / (wb**2 * tau * PI**(3.0/2.0))
    END SELECT
  END ASSOCIATE

  CALL PrintOption('Intensity amplitude: I0 [W/m^2]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%IntensityAmplitude)
ELSEIF(Species(iSpec)%Init(iInit)%Energy.GT.0.0)THEN
  ! Check if more than one pulse is required
  IF(Species(iSpec)%Init(iInit)%NbrOfPulses.GT.1)THEN
    Species(iSpec)%Init(iInit)%RepetitionRate     = GETREAL('Part-Species'//TRIM(hilf2)//'-RepetitionRate')
    Species(iSpec)%Init(iInit)%Period = 1./Species(iSpec)%Init(iInit)%RepetitionRate
  ELSE
    Species(iSpec)%Init(iInit)%Period = 2.0 * Species(iSpec)%Init(iInit)%tShift
  END IF ! Species(iSpec)%Init(iInit)%NbrOfPulses
  LBWRITE(*,*) 'Photoionization in cylindrical/rectangular/honeycomb volume: Selecting mode [Energy]'

  ASSOCIATE( E0   => Species(iSpec)%Init(iInit)%Energy             ,&
             wb   => Species(iSpec)%Init(iInit)%WaistRadius        ,&
             tau  => Species(iSpec)%Init(iInit)%PulseDuration      ,&
             Rout => Species(iSpec)%Init(iInit)%RadiusIC           ,&
             Rin  => Species(iSpec)%Init(iInit)%Radius2IC          ,&
             A    => Species(iSpec)%Init(iInit)%Area               )
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('photon_SEE_rectangle','photon_rectangle')
      Species(iSpec)%Init(iInit)%IntensityAmplitude = E0 / (SQRT(PI)*tau*A)
    CASE('photon_SEE_honeycomb','photon_honeycomb')
      Species(iSpec)%Init(iInit)%IntensityAmplitude = E0 / (SQRT(PI)*tau*(1.5*SQRT(3.0))*(Rout**2-Rin**2))
    CASE DEFAULT
      Species(iSpec)%Init(iInit)%IntensityAmplitude = E0 / (wb**2 * tau * PI**(3.0/2.0))
    END SELECT
  END ASSOCIATE

  CALL PrintOption('Intensity amplitude: I0 [W/m^2]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%IntensityAmplitude)
ELSEIF(Species(iSpec)%Init(iInit)%IntensityAmplitude.GT.0.0)THEN
  ! Check if more than one pulse is required
  IF(Species(iSpec)%Init(iInit)%NbrOfPulses.GT.1)THEN
    Species(iSpec)%Init(iInit)%RepetitionRate     = GETREAL('Part-Species'//TRIM(hilf2)//'-RepetitionRate')
    Species(iSpec)%Init(iInit)%Period = 1./Species(iSpec)%Init(iInit)%RepetitionRate
  ELSE
    Species(iSpec)%Init(iInit)%Period = 2.0 * Species(iSpec)%Init(iInit)%tShift
  END IF ! Species(iSpec)%Init(iInit)%NbrOfPulses
  LBWRITE(*,*) 'Photoionization in cylindrical/rectangular/honeycomb volume: Selecting mode [IntensityAmplitude]'

  ! Calculate energy: E = I0*w_b**2*tau*PI**(3.0/2.0)
  ASSOCIATE( I0   => Species(iSpec)%Init(iInit)%IntensityAmplitude ,&
             wb   => Species(iSpec)%Init(iInit)%WaistRadius        ,&
             tau  => Species(iSpec)%Init(iInit)%PulseDuration      ,&
             Rout => Species(iSpec)%Init(iInit)%RadiusIC           ,&
             Rin  => Species(iSpec)%Init(iInit)%Radius2IC          ,&
             A    => Species(iSpec)%Init(iInit)%Area               )
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('photon_SEE_rectangle','photon_rectangle')
      Species(iSpec)%Init(iInit)%Energy = I0 * SQRT(PI)*tau*A
    CASE('photon_SEE_honeycomb','photon_honeycomb')
      Species(iSpec)%Init(iInit)%Energy = I0 * SQRT(PI)*tau*(1.5*SQRT(3.0))*(Rout**2-Rin**2)
    CASE DEFAULT
      Species(iSpec)%Init(iInit)%Energy = I0 * wb**2*tau*PI**(3.0/2.0)
    END SELECT
  END ASSOCIATE

  CALL PrintOption('Single pulse energy: Part-Species'//TRIM(hilf2)//'-Energy [J]','CALCUL.',&
                    RealOpt=Species(iSpec)%Init(iInit)%Energy)
ELSE
  CALL abort(__STAMP__,'Photoionization in cylinder: Either power P and repetition rate f, or energy E or intensity maximum I0!')
END IF ! use RepetitionRate and Power

! Sanity check: overlapping of pulses is not implemented (use multiple emissions for this)
IF(2.0*Species(iSpec)%Init(iInit)%tShift.GT.Species(iSpec)%Init(iInit)%Period) CALL abort(__STAMP__&
  ,'Pulse length (2*tShift) is greater than the pulse period. This is not implemented!')

! Calculate the corrected intensity amplitude I0_corrected (due to temporal "-tShift to tShift" and spatial cut-off "0 to R")
ASSOCIATE( tShift => Species(iSpec)%Init(iInit)%tShift             ,&
           wb     => Species(iSpec)%Init(iInit)%WaistRadius        ,&
           tau    => Species(iSpec)%Init(iInit)%PulseDuration      ,&
           Rout   => Species(iSpec)%Init(iInit)%RadiusIC           ,&
           Rin    => Species(iSpec)%Init(iInit)%Radius2IC          ,&
             A    => Species(iSpec)%Init(iInit)%Area               )
  SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
  CASE('photon_SEE_rectangle','photon_rectangle')
   ! no need for correction in space because the function is not cut-off in space
   ! just consider the temporal cut-off for the rectangle
    factor = ERF(tShift/tau)
    factor = SQRT(PI)*tau*A
  CASE('photon_SEE_honeycomb','photon_honeycomb')
   ! no need for correction in space because the function is not cut-off in space
   ! just consider the temporal cut-off for the hexagon
    factor = ERF(tShift/tau)
    factor = SQRT(PI)*tau*(1.5*SQRT(3.0))*(Rout**2-Rin**2) * factor
  CASE DEFAULT
    factor = (1.0-EXP(-Rout**2/(wb**2)))*ERF(tShift/tau)
    factor = PI**(3.0/2.0) * wb**2 * tau * factor
  END SELECT
END ASSOCIATE

Species(iSpec)%Init(iInit)%IntensityAmplitude = Species(iSpec)%Init(iInit)%Energy / factor
CALL PrintOption('Corrected Intensity amplitude: I0_corr [W/m^2]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%IntensityAmplitude)

CALL PrintOption('Pulse period (Time between maximum of two pulses) [s]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%Period)

CALL PrintOption('Temporal pulse width (pulse time 2x tShift) [s]','CALCUL.',RealOpt=2.0*Species(iSpec)%Init(iInit)%tShift)
Species(iSpec)%Init(iInit)%tActive = REAL(Species(iSpec)%Init(iInit)%NbrOfPulses - 1)*Species(iSpec)%Init(iInit)%Period &
                                        + 2.0*Species(iSpec)%Init(iInit)%tShift
CALL PrintOption('Pulse will end at tActive (pulse final time) [s]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%tActive)

! Read effective intensity for volume emission or yield for SEE (surface)
SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
CASE('photon_cylinder','photon_honeycomb','photon_rectangle')
  Species(iSpec)%Init(iInit)%EffectiveIntensityFactor = GETREAL('Part-Species'//TRIM(hilf2)//'-EffectiveIntensityFactor')
CASE DEFAULT
  Species(iSpec)%Init(iInit)%YieldSEE = GETREAL('Part-Species'//TRIM(hilf2)//'-YieldSEE')
END SELECT

END SUBROUTINE InitializeVariablesPhotoIonization


SUBROUTINE DetermineInitialParticleNumber()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: PI
USE MOD_DSMC_Vars           ,ONLY: RadialWeighting, DSMC
USE MOD_Particle_Mesh_Vars  ,ONLY: LocalVolume
USE MOD_Particle_Vars       ,ONLY: Species,nSpecies,SpecReset,Symmetry
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
            CALL abort(__STAMP__,'BaseVectors are parallel or zero!')
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__,'Given velocity distribution is not supported with the SpaceIC cuboid/sphere/cylinder!')
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
        IF(Symmetry%Order.LE.2) THEN
          ! The radial scaling of the weighting factor has to be considered
          IF(RadialWeighting%DoRadialWeighting) Species(iSpec)%Init(iInit)%ParticleNumber = &
                                      INT(Species(iSpec)%Init(iInit)%ParticleNumber * 2. / (RadialWeighting%PartScaleFactor),8)
        END IF
      CASE('background')
        ! do nothing
      CASE DEFAULT
        SWRITE(*,*) 'SpaceIC is: ', TRIM(Species(iSpec)%Init(iInit)%SpaceIC)
        CALL abort(__STAMP__,'ERROR: Unknown SpaceIC for species: ', iSpec)
      END SELECT    ! Species(iSpec)%Init(iInit)%SpaceIC
    END IF
    ! Sum-up the number of particles to be inserted
#if USE_MPI
    insertParticles = insertParticles + INT(REAL(Species(iSpec)%Init(iInit)%ParticleNumber)/REAL(nProcessors),8)
#else
    insertParticles = insertParticles + INT(Species(iSpec)%Init(iInit)%ParticleNumber,8)
#endif
  END DO ! iInit = 1, Species(iSpec)%NumberOfInits
END DO ! iSpec=1,nSpecies

END SUBROUTINE DetermineInitialParticleNumber


!===================================================================================================================================
!> Initialize the particle boundary IDs for all emission inits
!> This routine is only called for DoBoundaryParticleOutputHDF5=T, which is determined in particle boundary init that comes after
!> the general variable species emission initialization
!===================================================================================================================================
SUBROUTINE InitializeVariablesSpeciesBoundary()
! MODULES
USE MOD_ReadInTools   ,ONLY: GETINT
USE MOD_Particle_Vars ,ONLY: Species,nSpecies
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iSpec, iInit
CHARACTER(32) :: hilf, hilf2
!===================================================================================================================================

! Loop all species
DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  ! Loop all inits
  DO iInit = 1, Species(iSpec)%NumberOfInits
    WRITE(UNIT=hilf2,FMT='(I0)') iInit
    hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
    ! Read-in of type and particle number for emission per iteration
    Species(iSpec)%Init(iInit)%PartBCIndex = GETINT('Part-Species'//TRIM(hilf2)//'-PartBCIndex')
  END DO ! iInit = 1, Species(iSpec)%NumberOfInits
END DO ! iSpec = 1, nSpecies

END SUBROUTINE InitializeVariablesSpeciesBoundary


!===================================================================================================================================
!> Initialize emission-specific macro particle factors when Part-vMPF=T. The default value is the species MPF.
!===================================================================================================================================
SUBROUTINE InitializeEmissionSpecificMPF()
! MODULES
USE MOD_ReadInTools   ,ONLY: GETREAL,GETINT
USE MOD_Particle_Vars ,ONLY: Species,nSpecies
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iSpec, iInit
CHARACTER(32) :: hilf, hilf2, hilf3
!===================================================================================================================================

! Loop all species
DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  ! Loop all inits
  DO iInit = 1, Species(iSpec)%NumberOfInits
    WRITE(UNIT=hilf2,FMT='(I0)') iInit
    hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
    ! Read-in of type and particle number for emission per iteration
    WRITE(UNIT=hilf3,FMT='(G0)') Species(iSpec)%MacroParticleFactor
    Species(iSpec)%Init(iInit)%MacroParticleFactor = GETREAL('Part-Species'//TRIM(hilf2)//'-MacroParticleFactor',TRIM(hilf3))
  END DO ! iInit = 1, Species(iSpec)%NumberOfInits
END DO ! iSpec = 1, nSpecies

END SUBROUTINE InitializeEmissionSpecificMPF


SUBROUTINE ReadUseEmissionDistribution()
!===================================================================================================================================
! ATTENTION: The fields (density, temperature and velocity) need to be defined on equidistant data-points as either .csv or .h5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionFileName
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionDim,EmissionDistributionAxisSym
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionDelta,EmissionDistributionDim
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionMin
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionMax,EmissionDistributionNum
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionRadInd,EmissionDistributionAxisDir
USE MOD_HDF5_Input_Field       ,ONLY: ReadExternalFieldFromHDF5
USE MOD_Particle_Vars          ,ONLY: Species,nSpecies
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER     :: lenmin=4
INTEGER               :: lenstr
INTEGER               :: iSpec,iInit
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A,3X,A,65X,A)') ' INITIALIZATION OF EMISSION DISTRIBUTION FOR PARTICLES '

! Check if file exists
IF(.NOT.FILEEXISTS(EmissionDistributionFileName)) CALL abort(__STAMP__,"File not found: "//TRIM(EmissionDistributionFileName))

! Check length of file name
lenstr=LEN(TRIM(EmissionDistributionFileName))
IF(lenstr.LT.lenmin) CALL abort(__STAMP__,"File name too short: "//TRIM(EmissionDistributionFileName))

! Check file ending, either .csv or .h5
IF(TRIM(EmissionDistributionFileName(lenstr-lenmin+2:lenstr)).EQ.'.h5')THEN
  DO iSpec = 1,nSpecies
    DO iInit = 1, Species(iSpec)%NumberOfInits
      IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'EmissionDistribution')THEN
        CALL ReadExternalFieldFromHDF5(TRIM(Species(iSpec)%Init(iInit)%EmissionDistributionName)                                ,&
                                            Species(iSpec)%Init(iInit)%EmissionDistribution        , EmissionDistributionDelta  ,&
            EmissionDistributionFileName , EmissionDistributionDim , EmissionDistributionAxisSym   , EmissionDistributionRadInd ,&
            EmissionDistributionAxisDir  , EmissionDistributionMin , EmissionDistributionMax       , EmissionDistributionNum     )
        IF(.NOT.ALLOCATED(Species(iSpec)%Init(iInit)%EmissionDistribution)) CALL abort(__STAMP__,&
            "Failed to load data from: "//TRIM(EmissionDistributionFileName))
      END IF ! TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'EmissionDistribution'
    END DO ! iInit = 1, Species(iSpec)%NumberOfInits
  END DO ! iSpec = 1,nSpecies
ELSEIF(TRIM(EmissionDistributionFileName(lenstr-lenmin+1:lenstr)).EQ.'.csv')THEN
  CALL abort(__STAMP__,'ReadUseEmissionDistribution(): Read-in from .csv is not implemented')
ELSE
  CALL abort(__STAMP__,"Unrecognised file format for : "//TRIM(EmissionDistributionFileName))
END IF

LBWRITE(UNIT_stdOut,'(A)')' ...EMISSION DISTRIBUTION INITIALIZATION DONE'
END SUBROUTINE ReadUseEmissionDistribution


END MODULE MOD_Particle_Emission_Init
