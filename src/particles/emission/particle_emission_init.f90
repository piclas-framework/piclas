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
PUBLIC :: DefineParametersParticleEmission, InitializeVariablesSpeciesInits, InitializeParticleEmission
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
CALL prms%CreateRealOption(     'Part-Species[$]-ChargeIC'  &
                                , '[TODO-DEFINE-PARAMETER]\n'//&
                                  'Particle Charge (without MPF) of species[$] dim' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MassIC'  &
                                , 'Particle Mass (without MPF) of species [$] [kg]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MacroParticleFactor' &
                                , 'Number of Microparticle per Macroparticle for species [$]', '1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-IsImplicit'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Flag if specific particle is implicit', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-IsIMDSpecies' &
                                , 'TODO-DEFINE-PARAMETER', '.FALSE.', numberedmulti=.TRUE.)

CALL prms%SetSection("Particle Initialization (Inits)")

CALL prms%CreateLogicalOption(  'Part-Species[$]-Init[$]-UseForInit' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Flag to use Init/Emission for init', '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Init[$]-UseForEmission' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Flag to use Init/Emission for emission', '.FALSE.', numberedmulti=.TRUE.)
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
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-initialParticleNumber' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Number of Particles at time 0.0', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-RadiusIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Radius for IC circle', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-Radius2IC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Radius2 for IC cylinder (ring)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-RadiusICGyro' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Radius for Gyrotron gyro radius', '1.', numberedmulti=.TRUE.)
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
CALL prms%CreateLogicalOption(  'Part-Species[$]-Init[$]-CalcHeightFromDt'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Calculate cuboid/cylinder height from v and dt?', '.FALSE.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-VeloIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Velocity for inital Data', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Init[$]-VeloVecIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Normalized velocity vector', '0. , 0. , 0.', numberedmulti=.TRUE.)
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
                                  'WaveNumber for sin-deviation initiation.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-MWTemperatureIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Temperature for Maxwell Distribution', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-PartDensity' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'PartDensity (real particles per m^3)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Init[$]-ParticleEmissionType'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Emission Type \n'//&
                                  '1 = emission rate in 1/s,\n'//&
                                  '2 = emission rate 1/iteration\n', '2', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-ParticleEmission' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Emission in [1/s] or [1/Iteration]', '0.', numberedmulti=.TRUE.)
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
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
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
INTEGER               :: iSpec, iInit, iExclude
CHARACTER(32)         :: hilf , hilf2, hilf3
LOGICAL               :: PartDens_OnlyInit
REAL                  :: lineVector(3), v_drift_line, A_ins, factor
!===================================================================================================================================
BGGas%NumberOfSpecies = 0
ALLOCATE(BGGas%BackgroundSpecies(nSpecies))
BGGas%BackgroundSpecies = .FALSE.
ALLOCATE(BGGas%NumberDensity(nSpecies))
BGGas%NumberDensity = 0.
ALLOCATE(SpecReset(1:nSpecies))
SpecReset=.FALSE.

DO iSpec = 1, nSpecies
  SWRITE (UNIT_stdOut,'(66(". "))')
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(hilf)//'-nInits','0')
#if USE_MPI
  IF(.NOT.PerformLoadBalance) THEN
#endif /*USE_MPI*/
    SpecReset(iSpec)                     = GETLOGICAL('Part-Species'//TRIM(hilf)//'-Reset','.FALSE.')
#if USE_MPI
  END IF
#endif /*USE_MPI*/
  Species(iSpec)%ChargeIC              = GETREAL('Part-Species'//TRIM(hilf)//'-ChargeIC','0.')
  Species(iSpec)%MassIC                = GETREAL('Part-Species'//TRIM(hilf)//'-MassIC')
  Species(iSpec)%MacroParticleFactor   = GETREAL('Part-Species'//TRIM(hilf)//'-MacroParticleFactor','1.')
#if defined(IMPA)
  Species(iSpec)%IsImplicit            = GETLOGICAL('Part-Species'//TRIM(hilf)//'-IsImplicit','.FALSE.')
#endif
  ALLOCATE(Species(iSpec)%Init(1:Species(iSpec)%NumberOfInits))
  DO iInit = 1, Species(iSpec)%NumberOfInits
    WRITE(UNIT=hilf2,FMT='(I0)') iInit
    hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-SpaceIC','cell_local'))
    !-------------------------------------------------------------------------------------------------------------------------------
    IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
      ! get emission and init data
      Species(iSpec)%Init(iInit)%UseForInit            = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForInit')
      Species(iSpec)%Init(iInit)%UseForEmission        = .FALSE.
    ELSE ! SpaceIC not cell_local
      Species(iSpec)%Init(iInit)%UseForInit            = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForInit')
      Species(iSpec)%Init(iInit)%UseForEmission        = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForEmission')
      IF((Symmetry%Order.EQ.2).OR.VarTimeStep%UseVariableTimeStep) THEN
        CALL abort(__STAMP__&
            ,'ERROR: Particle insertion/emission for 2D/axisymmetric or variable time step only possible with'//&
             'cell_local-SpaceIC and/or surface flux!')
      END IF
    END IF
    !-------------------------------------------------------------------------------------------------------------------------------
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution'&
      ,'maxwell_lpn'))
   IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'photon_SEE_energy')THEN
      Species(iSpec)%Init(iInit)%WorkFunctionSEE  = GETREAL('Part-Species'//TRIM(hilf2)//'-WorkFunctionSEE')
    END IF
    Species(iSpec)%Init(iInit)%InflowRiseTime        = GETREAL('Part-Species'//TRIM(hilf2)//'-InflowRiseTime','0.')
    Species(iSpec)%Init(iInit)%initialParticleNumber = GETINT('Part-Species'//TRIM(hilf2)//'-initialParticleNumber','0')
    Species(iSpec)%Init(iInit)%RadiusIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusIC','1.')
    Species(iSpec)%Init(iInit)%Radius2IC             = GETREAL('Part-Species'//TRIM(hilf2)//'-Radius2IC','0.')
    Species(iSpec)%Init(iInit)%RadiusICGyro          = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusICGyro','1.')
    Species(iSpec)%Init(iInit)%NormalIC              = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-NormalIC',3,'0. , 0. , 1.')
    Species(iSpec)%Init(iInit)%BasePointIC           = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BasePointIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector1IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector1IC',3,'1. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector2IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector2IC',3,'0. , 1. , 0.')
    Species(iSpec)%Init(iInit)%CuboidHeightIC        = GETREAL('Part-Species'//TRIM(hilf2)//'-CuboidHeightIC','1.')
    Species(iSpec)%Init(iInit)%CylinderHeightIC      = GETREAL('Part-Species'//TRIM(hilf2)//'-CylinderHeightIC','1.')
    Species(iSpec)%Init(iInit)%CalcHeightFromDt      = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CalcHeightFromDt','.FALSE.')
    Species(iSpec)%Init(iInit)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC','0.')
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%Amplitude             = GETREAL('Part-Species'//TRIM(hilf2)//'-Amplitude','0.01')
    Species(iSpec)%Init(iInit)%WaveNumber            = GETREAL('Part-Species'//TRIM(hilf2)//'-WaveNumber','2.')
    Species(iSpec)%Init(iInit)%maxParticleNumberX    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-x','0')
    Species(iSpec)%Init(iInit)%maxParticleNumberY    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-y','0')
    Species(iSpec)%Init(iInit)%maxParticleNumberZ    = GETINT('Part-Species'//TRIM(hilf2)//'-maxParticleNumber-z','0')
    Species(iSpec)%Init(iInit)%Alpha                 = GETREAL('Part-Species'//TRIM(hilf2)//'-Alpha','0.')
    Species(iSpec)%Init(iInit)%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC','0.')
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity','0.')
    ! Background gas definition
    IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'background') THEN
      IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN
        BGGas%NumberOfSpecies = BGGas%NumberOfSpecies + 1
        BGGas%BackgroundSpecies(iSpec)  = .TRUE.
        BGGas%NumberDensity(iSpec)      = Species(iSpec)%Init(iInit)%PartDensity
      ELSE
        CALL abort(__STAMP__&
            ,'ERROR: Only one background definition per species is allowed!')
      END IF
    END IF
    IF (Species(iSpec)%Init(iInit)%UseForEmission) THEN
      Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(hilf2)//'-ParticleEmissionType','2')
      Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL('Part-Species'//TRIM(hilf2)//'-ParticleEmission','0.')
    ELSE
      Species(iSpec)%Init(iInit)%ParticleEmissionType  = 0 !dummy
      Species(iSpec)%Init(iInit)%ParticleEmission      = 0. !dummy
    END IF
    ! Photoionization in cylinderical volume (modelling a laser pulse) and SEE based on photon impact on a surface
    IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_SEE_disc')          &
   .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_cylinder')) THEN
      Species(iSpec)%Init(iInit)%ParticleEmissionType = 7
      Species(iSpec)%Init(iInit)%UseForEmission = .TRUE.
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
        SWRITE(*,*) 'Photoionization in cylinderical volume: Selecting mode [RepetitionRate and Power]'

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
        SWRITE(*,*) 'Photoionization in cylinderical volume: Selecting mode [Energy]'

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
        SWRITE(*,*) 'Photoionization in cylinderical volume: Selecting mode [IntensityAmplitude]'

        ! Calculate energy: E = I0*w_b**2*tau*PI**(3.0/2.0)
        Species(iSpec)%Init(iInit)%Energy = Species(iSpec)%Init(iInit)%IntensityAmplitude*Species(iSpec)%Init(iInit)%WaistRadius**2&
                                            *Species(iSpec)%Init(iInit)%PulseDuration*PI**(3.0/2.0)
        CALL PrintOption('Single pulse energy: Part-Species'//TRIM(hilf2)//'-Energy [J]','CALCUL.',&
                         RealOpt=Species(iSpec)%Init(iInit)%Energy)
      ELSE
        CALL abort(&
          __STAMP__&
          ,'Photoionization in cylinderical volume: Supply either power P and repetition rate f, or energy E or intensity maximum I0!')
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

      CALL PrintOption('Temporal pulse width (pulse time) [s]','CALCUL.',RealOpt=2.0*Species(iSpec)%Init(iInit)%tShift)
      Species(iSpec)%Init(iInit)%tActive = REAL(Species(iSpec)%Init(iInit)%NbrOfPulses - 1)*Species(iSpec)%Init(iInit)%Period &
                                             + 2.0*Species(iSpec)%Init(iInit)%tShift
      CALL PrintOption('Pulse will end at tActive (pulse final time) [s]','CALCUL.',RealOpt=Species(iSpec)%Init(iInit)%tActive)

      IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_cylinder') THEN
        Species(iSpec)%Init(iInit)%EffectiveIntensityFactor = GETREAL('Part-Species'//TRIM(hilf2)//'-EffectiveIntensityFactor')
      ELSE
        Species(iSpec)%Init(iInit)%YieldSEE           = GETREAL('Part-Species'//TRIM(hilf2)//'-YieldSEE')
      END IF
    END IF
    Species(iSpec)%Init(iInit)%NumberOfExcludeRegions= GETINT('Part-Species'//TRIM(hilf2)//'-NumberOfExcludeRegions','0')
    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus = 0

    !----------- various checks/calculations after read-in of Species(i)%Init(iInit)%-data ----------------------------------!
    !--- Check if Initial ParticleInserting is really used
    !IF ( ((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.1).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2)) &
    !  .AND.
    IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
      IF ( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0) .AND. (Species(iSpec)%Init(iInit)%PartDensity.EQ.0.)) THEN
        Species(iSpec)%Init(iInit)%UseForInit=.FALSE.
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled as neither ParticleNumber"
        SWRITE(*,*) "nor PartDensity detected for Species, Init ", iSpec, iInit
      END IF
    END IF
    !--- cuboid-/cylinder-height calculation from v and dt
    IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid')) THEN
        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CuboidHeightIC,-1.)) THEN ! flag is initialized with -1, compatibility issue
          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.
          SWRITE(*,*) "WARNING: Cuboid height will be calculated from v and dt!"
        END IF
      ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CylinderHeightIC,-1.)) THEN !flag is initialized with -1, compatibility issue
          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.
          SWRITE(*,*) "WARNING: Cylinder height will be calculated from v and dt!"
        END IF
      END IF
    END IF
    IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) .AND. (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2) ) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is only supported for EmiType1 or EmiType2(=default)!')
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cuboid')        .AND.&
          (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cylinder'))          &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is only supported for cuboid or cylinder!')
      IF (Species(iSpec)%Init(iInit)%UseForInit) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is not supported for initial ParticleInserting!')
    END IF
    ! 1D Simulation with cuboid-SpaceIC
    IF((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').AND.Symmetry%Order.EQ.1) THEN
      IF(Species(iSpec)%Init(iInit)%BasePointIC(2).NE.-0.5 &
         .AND.Species(iSpec)%Init(iInit)%BasePointIC(3).NE.-0.5 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector1IC(2).NE.1 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector1IC(3).NE.0 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector2IC(1).NE.0 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector2IC(2).NE.0 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector1IC(1).NE.0 &
         .AND.Species(iSpec)%Init(iInit)%BaseVector2IC(3).NE.1 ) THEN
        SWRITE(*,*) 'For 1D Simulation with SpaceIC cuboid, the vectors have to be defined in the following from:'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-BasePointIC=(/x,-0.5,-0.5/), with x as the basepoint in x direction'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-BaseVector1IC=(/0.,1.,0/)'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-BaseVector2IC=(/0.,0.,1/)'
        SWRITE(*,*) 'Part-Species[$]-Init[$]-CuboidHeightIC is the extension of the insertion region and has to be positive'
        CALL abort(__STAMP__&
        ,'See above')
      END IF
    END IF
    !--- integer check for ParticleEmissionType 2
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         (ABS(Species(iSpec)%Init(iInit)%ParticleEmission-INT(Species(iSpec)%Init(iInit)%ParticleEmission,8)).GT.0.0)) THEN
       CALL abort(&
__STAMP__&
       ,' If ParticleEmissionType = 2 (parts per iteration), ParticleEmission has to be an integer number')
    END IF
    !--- normalize VeloVecIC and NormalIC (and BaseVector 1 & 2 IC for cylinder) for Inits
    IF (.NOT. ALL(Species(iSpec)%Init(iInit)%VeloVecIC(:).eq.0.)) THEN
      Species(iSpec)%Init(iInit)%VeloVecIC = Species(iSpec)%Init(iInit)%VeloVecIC            / &
        SQRT(Species(iSpec)%Init(iInit)%VeloVecIC(1)*Species(iSpec)%Init(iInit)%VeloVecIC(1) + &
        Species(iSpec)%Init(iInit)%VeloVecIC(2)*Species(iSpec)%Init(iInit)%VeloVecIC(2)      + &
        Species(iSpec)%Init(iInit)%VeloVecIC(3)*Species(iSpec)%Init(iInit)%VeloVecIC(3))
    END IF
    Species(iSpec)%Init(iInit)%NormalIC = Species(iSpec)%Init(iInit)%NormalIC /                 &
      SQRT(Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1) + &
      Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2) + &
      Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')&
        .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere')) THEN
        Species(iSpec)%Init(iInit)%BaseVector1IC =&
                  Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector1IC /     &
        SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)*Species(iSpec)%Init(iInit)%BaseVector1IC(1) + &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2)*Species(iSpec)%Init(iInit)%BaseVector1IC(2) + &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3)*Species(iSpec)%Init(iInit)%BaseVector1IC(3))
        Species(iSpec)%Init(iInit)%BaseVector2IC =&
                   Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector2IC /    &
        SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)*Species(iSpec)%Init(iInit)%BaseVector2IC(1) + &
        Species(iSpec)%Init(iInit)%BaseVector2IC(2)*Species(iSpec)%Init(iInit)%BaseVector2IC(2)      + &
        Species(iSpec)%Init(iInit)%BaseVector2IC(3)*Species(iSpec)%Init(iInit)%BaseVector2IC(3))
    END IF
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere')) THEN
        Species(iSpec)%Init(iInit)%NormalIC =&
                  Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%NormalIC /     &
        SQRT(Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1) + &
        Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2) + &
        Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
    END IF
    !--- read stuff for ExcludeRegions and normalize/calculate corresponding vectors
    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      ALLOCATE(Species(iSpec)%Init(iInit)%ExcludeRegion(1:Species(iSpec)%Init(iInit)%NumberOfExcludeRegions))
      IF (((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') &
      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') &
      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere'))) THEN
        DO iExclude=1,Species(iSpec)%Init(iInit)%NumberOfExcludeRegions
          WRITE(UNIT=hilf3,FMT='(I0)') iExclude
          hilf3=TRIM(hilf2)//'-ExcludeRegion'//TRIM(hilf3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC             &
            = TRIM(GETSTR('Part-Species'//TRIM(hilf3)//'-SpaceIC','cuboid'))
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC             &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-RadiusIC','1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC            &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-Radius2IC','0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC             &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-NormalIC',3,'0. , 0. , 1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC          &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BasePointIC',3,'0. , 0. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC        &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector1IC',3,'1. , 0. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC        &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector2IC',3,'0. , 1. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC       &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-CuboidHeightIC','1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC     &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-CylinderHeightIC','1.')
          !--normalize and stuff
          IF((TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).EQ.'cuboid') .OR. &
               ((((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1),1.) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2),0.)) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3),0.)) &
            .OR. ((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1),0.) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2),1.)) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3),0.))) &
            .AND. (((ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1),0.)) &
              .AND. (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2),0.))) &
              .AND. (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3),1.))))) THEN
            !-- cuboid; or BV are non-default and NormalIC is default: calc. NormalIC for ExcludeRegions from BV1/2
            !   (for def. BV and non-def. NormalIC; or all def. or non-def.: Use User-defined NormalIC when ExclRegion is cylinder)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3) &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)
          ELSE IF ( (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cuboid')        .AND. &
                    (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cylinder') )THEN
            CALL abort(&
__STAMP__&
,'Error in ParticleInit, ExcludeRegions must be cuboid or cylinder!')
          END IF
          IF (Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 + &
              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 + &
              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2 .GT. 0.) THEN
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC &
              / SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1) &
              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)**2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2) &
              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)**2)
          ELSE
            CALL abort(&
__STAMP__&
,'Error in ParticleInit, NormalIC Vector must not be zero!')
          END IF
        END DO !iExclude
      ELSE
        CALL abort(&
__STAMP__&
,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid, sphere or cylinder!')
      END IF
    END IF
    !--- stuff for calculating ParticleEmission/InitialParticleNumber from PartDensity
    PartDens_OnlyInit=.FALSE.
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.)) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) THEN
        IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2 .OR. Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0) &
            .AND. (Species(iSpec)%Init(iInit)%UseForInit) ) THEN
          PartDens_OnlyInit=.TRUE.
        ELSE
          CALL abort(&
            __STAMP__&
            , 'PartDensity is only supported for EmiType1 or initial ParticleInserting with EmiType1/2!')
        END IF
      END IF
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere')) THEN
        IF  ((((TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') ) ) THEN
          IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
            CALL abort(&
              __STAMP__&
              ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
          END IF
          !---calculation of Base-Area and corresponding component of VeloVecIC
          lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
            Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
          lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
            Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
          lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
            Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
          A_ins = lineVector(1)*lineVector(1) + lineVector(2)*lineVector(2) + lineVector(3)*lineVector(3)
          IF (A_ins .GT. 0.) THEN
            A_ins = SQRT(A_ins)
            lineVector = lineVector / A_ins
            IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
              v_drift_line = Species(iSpec)%Init(iInit)%VeloIC * &
                ( Species(iSpec)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(iSpec)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
                + Species(iSpec)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
            ELSE
              v_drift_line = 0.
              IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
                PartDens_OnlyInit=.TRUE.
              ELSE
                CALL abort(&
                  __STAMP__&
                  ,'PartDensity is only supported for CalcHeightFromDt, or initial ParticleInserting!')
              END IF
            END IF
            IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .EQ. 'cylinder' ) THEN
              A_ins = Pi * (Species(iSpec)%Init(iInit)%RadiusIC**2-Species(iSpec)%Init(iInit)%Radius2IC**2)
            END IF
            !---calculation of particle flow (macroparticles/s) through boundary
            IF (.NOT.PartDens_OnlyInit) THEN
              Species(iSpec)%Init(iInit)%ParticleEmission &
                = Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor * v_drift_line * A_ins
            END IF
            !---calculation of initial (macro)particle number
            IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
              IF (Species(iSpec)%Init(iInit)%initialParticleNumber .GT. 0) THEN
                CALL abort(&
                  __STAMP__&
                  ,'Either initialParticleNumber or PartDensity can be defined for selected parameters, not both!')
              END IF
              IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * Species(iSpec)%Init(iInit)%CuboidHeightIC * A_ins)
              ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * Species(iSpec)%Init(iInit)%CylinderHeightIC * A_ins)
              ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'sphere') THEN
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * (Species(iSpec)%Init(iInit)%RadiusIC**3 * 4./3. * PI))
              END IF
            END IF
          ELSE
            CALL abort(&
              __STAMP__&
              ,'BaseVectors are parallel or zero!')
          END IF
        ELSE
          CALL abort(&
            __STAMP__&
            ,'Only const. or maxwell(_lpn) is supported as velocityDistr. for PartDensity!')
        END IF
      ELSE IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cell_local')) THEN
        IF( (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
           .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') &
           .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') &
           .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'taylorgreenvortex') )THEN
          IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
            CALL abort(&
              __STAMP__&
              ,'Either ParticleEmission or PartDensity can be defined for cell_local emission parameters, not both!')
          END IF
          IF (LocalVolume.GT.0.) THEN
            IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
              IF (Species(iSpec)%Init(iInit)%initialParticleNumber .GT. 0) THEN
                CALL abort(&
                  __STAMP__&
                  ,'Either initialParticleNumber or PartDensity can be defined for selected parameters, not both!')
              END IF
              Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = NINT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor * LocalVolume)
            END IF
          ELSE
            CALL abort(&
              __STAMP__&
              ,'Local mesh volume is zero!')
          END IF
        ELSE
          CALL abort(&
            __STAMP__&
            ,'Only const. or maxwell_lpn is supported as velocityDistr. using cell_local inserting with PartDensity!')
        END IF
      ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'background') THEN
        ! do nothing
      ELSE
        SWRITE(*,*) 'SpaceIC is: ', TRIM(Species(iSpec)%Init(iInit)%SpaceIC)
        CALL abort(&
          __STAMP__&
          ,'ERROR: Unknown SpaceIC for species: ', iSpec)
      END IF
    END IF
    IF(Species(iSpec)%Init(iInit)%InflowRiseTime.GT.0.)THEN
      IF(.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow)  CALL CollectiveStop(&
__STAMP__, &
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


SUBROUTINE InitializeParticleEmission()
!===================================================================================================================================
! Initialize particles / Insert initial particles
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Vars   ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Dielectric_Vars     ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
USE MOD_DSMC_Vars           ,ONLY: useDSMC, RadialWeighting, DSMC
USE MOD_Part_Emission_Tools ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SetParticleTimeStep
USE MOD_Part_Pos_and_Velo   ,ONLY: SetParticlePosition,SetParticleVelocity, AD_SetInitElectronVelo
USE MOD_Part_Tools          ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Mesh_Vars  ,ONLY: LocalVolume
USE MOD_Particle_Vars       ,ONLY: Species,nSpecies,PDM,PEM, usevMPF, SpecReset, Symmetry, VarTimeStep
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
INTEGER               :: i, NbrOfParticle,iInit,iPart,PositionNbr
INTEGER(KIND=8)       :: insertParticles
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' Initial particle inserting... '

CALL UpdateNextFreePosition()

! Do sanity check of max. particle number compared to the number that is to be inserted for certain insertion types
insertParticles = 0
DO i=1,nSpecies
  IF (DSMC%DoAmbipolarDiff) THEN
    IF (i.EQ.DSMC%AmbiDiffElecSpec) CYCLE
  END IF
  IF (DoRestart .AND. .NOT.SpecReset(i)) CYCLE
  DO iInit = 1, Species(i)%NumberOfInits
    IF (TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
      IF(Symmetry%Order.LE.2) THEN
        ! The correct 2D/axisymmetric LocalVolume could only be calculated after the symmetry axis was defined (through the boundary
        ! conditions). However, the initialParticleNumber was already determined before the 2D volume calculation was performed.
        ! This can lead to initialParticleNumbers of 0, thus skipping the insertion entirely.
        Species(i)%Init(iInit)%initialParticleNumber &
                  = NINT(Species(i)%Init(iInit)%PartDensity / Species(i)%MacroParticleFactor * LocalVolume)
        ! The radial scaling of the weighting factor has to be considered
        IF(RadialWeighting%DoRadialWeighting) Species(i)%Init(iInit)%initialParticleNumber = &
                                    INT(Species(i)%Init(iInit)%initialParticleNumber * 2. / (RadialWeighting%PartScaleFactor),8)
      END IF
#if USE_MPI
      insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs,8)
#else
      insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber,8)
#endif
    ELSE IF ((TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cuboid') &
         .OR.(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
#if USE_MPI
      insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs,8)
#else
      insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber,8)
#endif
    END IF
  END DO
END DO
IF (insertParticles.GT.PDM%maxParticleNumber) THEN
  IPWRITE(UNIT_stdOut,*)' Maximum particle number : ',PDM%maxParticleNumber
  IPWRITE(UNIT_stdOut,*)' To be inserted particles: ',INT(insertParticles,4)
  CALL abort(&
__STAMP__&
,'Number of to be inserted particles per init-proc exceeds max. particle number! ')
END IF
DO i = 1,nSpecies
  IF (DSMC%DoAmbipolarDiff) THEN
    IF (i.EQ.DSMC%AmbiDiffElecSpec) CYCLE
  END IF
  IF (DoRestart .AND. .NOT.SpecReset(i)) CYCLE
  DO iInit = 1, Species(i)%NumberOfInits
    IF (Species(i)%Init(iInit)%UseForInit) THEN ! no special emissiontype to be used
      IF(Species(i)%Init(iInit)%initialParticleNumber.GT.HUGE(1)) CALL abort(&
__STAMP__&
,' Integer of initial particle number larger than max integer size: ',HUGE(1))
      NbrOfParticle = INT(Species(i)%Init(iInit)%initialParticleNumber,4)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle position for species ',i,' ... '
      CALL SetParticlePosition(i,iInit,NbrOfParticle)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',i,' ... '
      CALL SetParticleVelocity(i,iInit,NbrOfParticle)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle charge and mass for species ',i,' ... '
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
      IF (VarTimeStep%UseVariableTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
      IF (useDSMC) THEN
        IF (DSMC%DoAmbipolarDiff) CALL AD_SetInitElectronVelo(i,iInit,NbrOfParticle)
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            PDM%PartInit(PositionNbr) = iInit
          ELSE
            CALL abort(&
            __STAMP__&
            ,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
          END IF
          iPart = iPart + 1
        END DO
      END IF
      !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      CALL UpdateNextFreePosition()
    END IF ! not Emissiontype 4
  END DO !inits
END DO ! species

!--- set last element to current element (needed when ParticlePush is not executed, e.g. "delay")
DO i = 1,PDM%ParticleVecLength
  PEM%LastGlobalElemID(i) = PEM%GlobalElemID(i)
END DO

!--- Remove particles from dielectric regions if DielectricNoParticles=.TRUE.
IF(DoDielectric)THEN
  IF(DielectricNoParticles)THEN
    DO i = 1,PDM%ParticleVecLength
      ! Remove particles in dielectric elements
      IF(isDielectricElem(PEM%LocalElemID(i)))THEN
        PDM%ParticleInside(i) = .FALSE.
      END IF
    END DO
  END IF
END IF

SWRITE(UNIT_stdOut,'(A)') ' ...DONE '

END SUBROUTINE InitializeParticleEmission

END MODULE MOD_Particle_Emission_Init