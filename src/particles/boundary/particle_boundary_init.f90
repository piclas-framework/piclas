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

MODULE MOD_Particle_Boundary_Init
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
PUBLIC :: DefineParametersParticleBoundary, InitializeVariablesPartBoundary, FinalizeParticleBoundary
PUBLIC :: InitParticleBoundaryRotPeriodic, InitAdaptiveWallTemp
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle boundaries
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBoundary()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Boundaries")

CALL prms%CreateIntOption(      'Part-RotPeriodicAxi' , 'Axis of rotational periodicity: x = 1, y = 2, z = 3')
CALL prms%CreateIntOption(      'Part-nBounds', 'Number of particle boundaries.', '1')
CALL prms%CreateStringOption(   'Part-Boundary[$]-SourceName', &
                                  'No Default. Source Name of Boundary[i]. Has to be selected for all'//&
                                  'nBounds. Has to be same name as defined in preproc tool', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Boundary[$]-Condition', &
                                  'Possible conditions for boundary[$] are:\n'//&
                                  '- open\n'//&
                                  '- reflective\n'//&
                                  '- periodic\n'//&
                                  '- simple_anode\n'//&
                                  '- simple_cathode.\n'//&
                                  '- rot_periodic.\n'//&
                                  'If condition=reflective (Part-Boundary[$]-=PB): PB-MomentumACC,PB-WallTemp,PB-TransACC,PB-VibACC,PB-RotACC,'//&
                                  'PB-WallVelo,SpeciesSwaps.\nIf condition=periodic:Part-nPeriodicVectors', 'open', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Boundary[$]-Dielectric' , 'Define if particle boundary [$] is a '//&
                              'dielectric interface, i.e. an interface between a dielectric and a non-dielectric or a between two'//&
                              ' different dielectrics [.TRUE.] or not [.FALSE.] (requires reflective BC and species swap for nSpecies)'&
                              , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Boundary[$]-BoundaryParticleOutput' , 'Define if the properties of particles impacting on '//&
                              'boundary [$] are to be stored in an additional .h5 file for post-processing analysis [.TRUE.] '//&
                              'or not [.FALSE.].'&
                              , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-UseAdaptedWallTemp', &
                                'Use an adapted wall temperature, calculated if Part-AdaptWallTemp = T, otherwise read-in '//&
                                'from the restart file','.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-RadiativeEmissivity',  &
                                'Radiative emissivity of the boundary used to determine the wall temperature assuming a '//&
                                'radiative equilibrium at the wall, Part-Boundary1-UseAdaptedWallTemp = T and ' //&
                                'Part-AdaptWallTemp = T must be set', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-WallTemp'  &
                                , 'Wall temperature (in [K]) of reflective particle boundary [$].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-MomentumACC'  &
                                , 'Momentum accommodation coefficient of reflective particle boundary [$].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-TransACC'  &
                                , 'Translation accommodation coefficient of reflective particle boundary [$].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-VibACC'  &
                                , 'Vibrational accommodation coefficient of reflective particle boundary [$].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-RotACC'  &
                                , 'Rotational accommodation coefficient of reflective particle boundary [$].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-ElecACC '  &
                                , 'Electronic accommodation coefficient of reflective particle boundary [$].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-Resample', &
                                  'Sample particle properties from equilibrium distribution after reflection', '.FALSE.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-WallVelo'  &
                                , 'Velocity (global x,y,z in [m/s]) of reflective particle boundary [$].' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-RotVelo'  &
                                , 'Flag for rotating walls:'//&
                                  ' Particles will be accelerated additionally to the boundary interaction'//&
                                  ' through the rotating wall depending on their POI, rotation frequency and rotation axis.'//&
                                  ' In that case Part-Boundary[$]-WallVelo will be overwritten.' &
                                , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-RotFreq'  &
                                , 'Rotation frequency of the wall in [Hz]. Note: Rotation direction based on right-hand rule!' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Boundary[$]-RotAxis'  &
                                , 'Definition of rotation axis, only major axis: x=1,y=2,z=3.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(      'Part-Boundary[$]-RotPeriodicAngle' , 'Angle and Direction of rotation periodicity, either + or - '//&
                                'Note: Rotation direction based on right-hand rule!', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(      'Part-Boundary[$]-RotPeriodicMin' , 'Minimum coordinate at rotational axis for this segment', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(      'Part-Boundary[$]-RotPeriodicMax' , 'Maximum coordinate at rotational axis for this segment', numberedmulti=.TRUE.)
!CALL prms%CreateLogicalOption(  'Part-RotPeriodicReBuild', 'Force re-creation of rotational periodic mapping (which might already exist in the mesh file).', '.FALSE.')
CALL prms%CreateRealOption(     'Part-Boundary[$]-WallTemp2'  &
                                , 'Second wall temperature (in [K]) of reflective particle boundary for a temperature gradient.' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-TempGradStart'  &
                                , 'Impose a temperature gradient by supplying a start/end vector and a second wall temperature.' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-TempGradEnd'  &
                                , 'Impose a temperature gradient by supplying a start/end vector and a second wall temperature.' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Boundary[$]-TempGradDir', 'Optional definition of the temperature '//&
                                'gradient direction along a major axis: x = 1, y = 2, z = 3. Default = 0: Gradient is along '//&
                                'the vector defined by the start and end values', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Boundary[$]-SurfaceModel'  &
                                , 'Defining surface to be treated reactively by defining Model used for particle surface interaction. If any >0 then look in section SurfaceModel.\n'//&
                                '0: Maxwell scattering\n'//&
                                '5: SEE-E and SEE-I (secondary e- emission due to e- or i+ bombardment) by Levko2015 for copper electrodes\n'//&
                                '6: SEE-E (secondary e- emission due to e- bombardment) by Pagonakis2016 for molybdenum, originally from Harrower1956. Currently not available\n'//&
                                '7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for '//&
                                'secondary e- emission with 0.13 probability) by D. Depla, "Magnetron sputter deposition: Linking discharge voltage with target properties", 2009\n'// &
                                '8: SEE-E (e- on dielectric materials is considered for SEE and three different outcomes) '//&
                                'by A.I. Morozov, "Structure of Steady-State Debye Layers in a Low-Density Plasma near a Dielectric Surface", 2004\n'//&
                                '9: SEE-I when Ar+ ion bombards surface with 0.01 probability and fixed SEE electron energy of 6.8 eV\n'//&
                                '10: SEE-I when Ar+ bombards copper by J.G. Theis "Computing the Paschen curve for argon with speed-limited particle-in-cell simulation", 2021 (originates from Phelps1999)\n'// &
                                '11: SEE-E when e- bombard quartz (SiO2) by A. Dunaevsky, "Secondary electron emission from dielectric materials of a Hall thruster with segmented electrodes", 2003'&
                                , '0', numberedmulti=.TRUE.)
CALL prms%SetSection('Particle Boundaries: Species Swap')
CALL prms%CreateIntOption(      'Part-Boundary[$]-NbrOfSpeciesSwaps'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Number of Species to be changed at wall.', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-ProbOfSpeciesSwaps'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Probability of SpeciesSwaps at wall', '1.', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Part-Boundary[$]-SpeciesSwaps[$]'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Species to be changed at wall (out=: delete)', '0 , 0'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-AdaptWallTemp','Perform wall temperature adaptation at every macroscopic output.', '.FALSE.')

END SUBROUTINE DefineParametersParticleBoundary


SUBROUTINE InitializeVariablesPartBoundary()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_ReadInTools
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas
USE MOD_Mesh_Vars              ,ONLY: BoundaryName,BoundaryType, nBCs
USE MOD_Particle_Vars          ,ONLY: PDM, nSpecies, PartMeshHasPeriodicBCs, RotRefFrameAxis, SpeciesDatabase, Species
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound,DoBoundaryParticleOutputHDF5,PartStateBoundary, AdaptWallTemp
USE MOD_Particle_Boundary_Vars ,ONLY: nVarPartStateBoundary
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Emission_Init ,ONLY: InitializeVariablesSpeciesBoundary
USE MOD_PICDepo_Vars           ,ONLY: DepositionType,DoHaloDepo
USE MOD_HDF5_input             ,ONLY: OpenDataFile, ReadArray, DatasetExists, GetDataSize, nDims, HSize, CloseDataFile
USE MOD_SurfaceModel_Vars      ,ONLY: StickingCoefficientData
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Vars          ,ONLY: PartMeshHasReflectiveBCs
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPartBound, iBC, iPBC, iSwaps, MaxNbrOfSpeciesSwaps, RotAxis, nRotPeriodicBCs, TempGradDir
INTEGER               :: ALLOCSTAT, dummy_int
REAL                  :: omegaTemp, RotFreq
CHARACTER(32)         :: hilf , hilf2
CHARACTER(200)        :: tmpString
CHARACTER(LEN=64)     :: dsetname
LOGICAL               :: StickingCoefficientExists,FoundPartBoundSEE
INTEGER               :: iInit,iSpec
!===================================================================================================================================
! Read in boundary parameters
dummy_int  = CountOption('Part-nBounds') ! check if Part-nBounds is present in .ini file
nPartBound = GETINT('Part-nBounds')      ! get number of particle boundaries
! Read-in number of porous boundaries
nPorousBC  = GETINT('Surf-nPorousBC')
IF ((nPartBound.LE.0).OR.(dummy_int.LT.0)) CALL abort(__STAMP__  ,'ERROR: nPartBound .LE. 0:', nPartBound)

ALLOCATE(PartBound%SourceBoundName(  1:nPartBound))
PartBound%SourceBoundName = ''
ALLOCATE(PartBound%TargetBoundCond(  1:nPartBound))
PartBound%TargetBoundCond = -1
ALLOCATE(PartBound%MomentumACC(      1:nPartBound))
PartBound%MomentumACC = -1
ALLOCATE(PartBound%OnlySpecular(     1:nPartBound))
PartBound%OnlySpecular = .FALSE.
ALLOCATE(PartBound%OnlyDiffuse(      1:nPartBound))
PartBound%OnlyDiffuse = .FALSE.
ALLOCATE(PartBound%WallTemp(         1:nPartBound))
PartBound%WallTemp = -1.
ALLOCATE(PartBound%WallTemp2(        1:nPartBound))
PartBound%WallTemp2 = -1.
ALLOCATE(PartBound%WallTempDelta(    1:nPartBound))
PartBound%WallTempDelta = 0.
ALLOCATE(PartBound%TransACC(         1:nPartBound))
PartBound%TransACC = -1.
ALLOCATE(PartBound%VibACC(           1:nPartBound))
PartBound%VibACC = -1.
ALLOCATE(PartBound%RotACC(           1:nPartBound))
PartBound%RotACC = -1.
ALLOCATE(PartBound%ElecACC(          1:nPartBound))
PartBound%ElecACC = -1.
ALLOCATE(PartBound%Resample(         1:nPartBound))
PartBound%Resample = .FALSE.
ALLOCATE(PartBound%WallVelo(     1:3,1:nPartBound))
PartBound%WallVelo = 0.
ALLOCATE(PartBound%RotVelo(          1:nPartBound))
PartBound%RotVelo = .FALSE.
ALLOCATE(PartBound%RotOmega(       1:3,1:nPartBound))
PartBound%RotOmega = 0.
ALLOCATE(PartBound%RotPeriodicAngle(  1:nPartBound))
PartBound%RotPeriodicAngle = 0
ALLOCATE(PartBound%RotPeriodicMin(  1:nPartBound))
PartBound%RotPeriodicMin = -HUGE(1.)
ALLOCATE(PartBound%RotPeriodicMax(  1:nPartBound))
PartBound%RotPeriodicMax = HUGE(1.)
ALLOCATE(PartBound%TempGradStart(1:3,1:nPartBound))
PartBound%TempGradStart = 0.
ALLOCATE(PartBound%TempGradEnd(  1:3,1:nPartBound))
PartBound%TempGradEnd = 0.
ALLOCATE(PartBound%TempGradVec(  1:3,1:nPartBound))
PartBound%TempGradVec = 0.
ALLOCATE(PartBound%TempGradDir(1:nPartBound))
PartBound%TempGradDir = 0
ALLOCATE(PartBound%SurfaceModel(     1:nPartBound))
PartBound%SurfaceModel = 0
ALLOCATE(PartBound%Reactive(         1:nPartBound))
PartBound%Reactive = .FALSE.
ALLOCATE(PartBound%NbrOfSpeciesSwaps(1:nPartBound))
PartBound%NbrOfSpeciesSwaps = 0
ALLOCATE(PartBound%UseAdaptedWallTemp(1:nPartBound))
PartBound%UseAdaptedWallTemp = .FALSE.
ALLOCATE(PartBound%RadiativeEmissivity(1:nPartBound))
PartBound%RadiativeEmissivity = 1.

! Output of wall temperature per default off
PartBound%OutputWallTemp = .FALSE.

!--determine MaxNbrOfSpeciesSwaps for correct allocation
MaxNbrOfSpeciesSwaps=0
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I0)') iPartBound
  PartBound%NbrOfSpeciesSwaps(iPartBound)= GETINT('Part-Boundary'//TRIM(hilf)//'-NbrOfSpeciesSwaps','0')
  MaxNbrOfSpeciesSwaps=max(PartBound%NbrOfSpeciesSwaps(iPartBound),MaxNbrOfSpeciesSwaps)
END DO
IF (MaxNbrOfSpeciesSwaps.gt.0) THEN
  ALLOCATE(PartBound%ProbOfSpeciesSwaps(1:nPartBound))
  PartBound%ProbOfSpeciesSwaps = -1
  ALLOCATE(PartBound%SpeciesSwaps(1:2,1:MaxNbrOfSpeciesSwaps,1:nPartBound))
  PartBound%SpeciesSwaps = -1
END IF
! Dielectric Surfaces
ALLOCATE(PartBound%Dielectric(1:nPartBound))
PartBound%Dielectric      = .FALSE.
DoDielectricSurfaceCharge = .FALSE.
DoHaloDepo                = .FALSE. ! dielectric surfaces or implicit particle deposition
! Surface particle output to .h5
ALLOCATE(PartBound%BoundaryParticleOutputHDF5(1:nPartBound))
PartBound%BoundaryParticleOutputHDF5=.FALSE.
DoBoundaryParticleOutputHDF5=.FALSE.

PartMeshHasPeriodicBCs=.FALSE.
GEO%RotPeriodicBC =.FALSE.
nRotPeriodicBCs  = 0
! TODO: REMOVE THIS CALL WHEN MERGED WITH UNIFIED SPECIES DATABASE BRANCH
SpeciesDatabase = GETSTR('Particles-Species-Database', 'none')

#if defined(IMPA) || defined(ROS)
PartMeshHasReflectiveBCs=.FALSE.
#endif
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I0)') iPartBound
  tmpString = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-Condition','open'))

  SELECT CASE (TRIM(tmpString))
  CASE('open')
    PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
  CASE('reflective')
#if defined(IMPA) || defined(ROS)
    PartMeshHasReflectiveBCs=.TRUE.
#endif
    PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
    PartBound%MomentumACC(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-MomentumACC')
    IF(PartBound%MomentumACC(iPartBound).EQ.0.0) THEN
      PartBound%OnlySpecular(iPartBound) = .TRUE.
    ELSE IF(PartBound%MomentumACC(iPartBound).EQ.1.0) THEN
      PartBound%OnlyDiffuse(iPartBound)  = .TRUE.
    END IF
    PartBound%WallTemp(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp')
    PartBound%TransACC(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-TransACC')
    PartBound%VibACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-VibACC')
    PartBound%RotACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotACC')
    PartBound%ElecACC(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-ElecACC')
    PartBound%Resample(iPartBound)        = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Resample')
    PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-WallVelo',3)
    PartBound%RotVelo(iPartBound)         = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-RotVelo')
    IF(PartBound%RotVelo(iPartBound)) THEN
      RotFreq                             = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotFreq')
      RotAxis                             = GETINT('Part-Boundary'//TRIM(hilf)//'-RotAxis')
      omegaTemp = 2. * PI * RotFreq
      SELECT CASE(RotAxis)
        CASE(1)
          PartBound%RotOmega(1:3,iPartBound) = (/omegaTemp,0.,0./)
        CASE(2)
          PartBound%RotOmega(1:3,iPartBound) = (/0.,omegaTemp,0./)
        CASE(3)
          PartBound%RotOmega(1:3,iPartBound) = (/0.,0.,omegaTemp/)
        CASE DEFAULT
          CALL abort(__STAMP__,'ERROR Rotational Wall Velocity: Axis must be between 1 and 3. Selected axis: ',IntInfoOpt=RotRefFrameAxis)
      END SELECT
    END IF
    ! Utilize an adaptive wall temparature
    PartBound%UseAdaptedWallTemp(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-UseAdaptedWallTemp')
    ! Activate wall temperature output in the DSMCSurfState, required for the initialization of the array as well
    IF(PartBound%UseAdaptedWallTemp(iPartBound)) PartBound%OutputWallTemp = .TRUE.
    PartBound%RadiativeEmissivity(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-RadiativeEmissivity')
    ! Selection of the surface model (e.q. SEE, sticking, etc.)
    PartBound%SurfaceModel(iPartBound)    = GETINT('Part-Boundary'//TRIM(hilf)//'-SurfaceModel')
    ! Impose a wall temperature gradient
    PartBound%WallTemp2(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp2')
    IF(PartBound%WallTemp2(iPartBound).GT.0.) THEN
      ! Activate wall temperature output in the DSMCSurfState
      PartBound%OutputWallTemp = .TRUE.
      PartBound%TempGradStart(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-TempGradStart',3)
      PartBound%TempGradEnd(1:3,iPartBound)   = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-TempGradEnd',3)
      PartBound%TempGradDir(iPartBound)   = GETINT('Part-Boundary'//TRIM(hilf)//'-TempGradDir')
      TempGradDir = PartBound%TempGradDir(iPartBound)
      IF((TempGradDir.LT.0).AND.(TempGradDir.GT.3)) THEN
        CALL abort(__STAMP__,'ERROR Wall Temperature Gradient: Input must be between 1 (=x), 2 (=y), 3 (=z) or 0. Input: ', &
                    IntInfoOpt=TempGradDir)
      END IF
      ! Calculate the magnitude of the temperature gradient
      PartBound%WallTempDelta(iPartBound)   = PartBound%WallTemp2(iPartBound) - PartBound%WallTemp(iPartBound)
      ! Determine the direction of the temperature gradient
      SELECT CASE(TempGradDir)
      CASE(0)
        PartBound%TempGradVec(1:3,iPartBound) = PartBound%TempGradEnd(1:3,iPartBound) - PartBound%TempGradStart(1:3,iPartBound)
      CASE(1,2,3)
        PartBound%TempGradVec(TempGradDir,iPartBound) = PartBound%TempGradEnd(TempGradDir,iPartBound) &
                                                        - PartBound%TempGradStart(TempGradDir,iPartBound)
      END SELECT
      ! Sanity check: defined vector shall be above zero
      IF(ALL(PartBound%TempGradVec(1:3,iPartBound).EQ.0.)) THEN
        SWRITE(*,*) 'ERROR Temperature gradient vector: ', PartBound%TempGradVec(1:3,iPartBound)
        CALL abort(__STAMP__,'ERROR Wall Temperature Gradient: gradient vector appears to be zero!')
      END IF
    END IF
    ! check for correct surfacemodel input
    IF (PartBound%SurfaceModel(iPartBound).GT.0)THEN
      IF (.NOT.useDSMC) CALL abort(__STAMP__,'Cannot use surfacemodel>0 with useDSMC=F for particle boundary: ',iPartBound)
      SELECT CASE (PartBound%SurfaceModel(iPartBound))
      CASE (0)
        PartBound%Reactive(iPartBound)        = .FALSE.
      CASE (1)
        PartBound%Reactive(iPartBound)        = .FALSE.
        IF(TRIM(SpeciesDatabase).EQ.'none') &
          CALL abort(__STAMP__,'ERROR in InitializeVariablesPartBoundary: SpeciesDatabase is required for the boundary #', iPartBound)
      CASE (SEE_MODELS_ID)
        ! SEE models require reactive BC
        PartBound%Reactive(iPartBound)        = .TRUE.
      CASE DEFAULT
        CALL abort(__STAMP__,'Error in particle init: only allowed SurfaceModels: 0,SEE_MODELS_ID! SurfaceModel=',&
        IntInfoOpt=PartBound%SurfaceModel(iPartBound))
      END SELECT
    END IF
    IF (PartBound%NbrOfSpeciesSwaps(iPartBound).GT.0) THEN
      !read Species to be changed at wall (in, out), out=0: delete
      PartBound%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-Boundary'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
      DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(iPartBound)
        WRITE(UNIT=hilf2,FMT='(I0)') iSwaps
        PartBound%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
            GETINTARRAY('Part-Boundary'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
      END DO
      IF(PartBound%Reactive(iPartBound)) CALL abort(__STAMP__&
          ,'ERROR: Species swap is only supported in combination with Maxwell scattering (SurfModel = 0). PartBound: ',iPartBound)
    END IF
    ! Dielectric Surfaces
    PartBound%Dielectric(iPartBound)      = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Dielectric')
    ! Sanity check: PartBound%Dielectric=T requires supplying species swap for every species
    IF(PartBound%Dielectric(iPartBound))THEN
      IF((PartBound%NbrOfSpeciesSwaps(iPartBound).LT.(nSpecies-BGGas%NumberOfSpecies)).AND.&
          (.NOT.PartBound%Reactive(iPartBound)))THEN
        CALL abort(__STAMP__,'PartBound%Dielectric=T requires\n   a) supplying a species swap (Part-BoundaryX-NbrOfSpeciesSwaps)'//&
            ' for every species (except background gas species) or\n   '//&
            'b) surface model that is reactive (Part-BoundaryX-SurfaceModel)!')
      ELSE
        DoDielectricSurfaceCharge = .TRUE.
        DoHaloDepo                = .TRUE.
        IF(TRIM(DepositionType).NE.'cell_volweight_mean') CALL CollectiveStop(__STAMP__,&
            'PartBound%Dielectric=T requires cell_volweight_mean (12) as deposition method')
      END IF ! PartBound%NbrOfSpeciesSwaps(iPartBound).NE.nSpecies
    END IF ! PartBound%Dielectric(iPartBound)
  CASE('periodic')
    PartBound%TargetBoundCond(iPartBound) = PartBound%PeriodicBC
    PartMeshHasPeriodicBCs = .TRUE.
  CASE('simple_anode')
    PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleAnodeBC
  CASE('simple_cathode')
    PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleCathodeBC
  CASE('symmetric')
#if defined(IMPA) || defined(ROS)
    PartMeshHasReflectiveBCs=.TRUE.
#endif
    PartBound%TargetBoundCond(iPartBound) = PartBound%SymmetryBC
    PartBound%WallVelo(1:3,iPartBound)    = (/0.,0.,0./)
  CASE('symmetric_axis')
    PartBound%TargetBoundCond(iPartBound) = PartBound%SymmetryAxis
    PartBound%WallVelo(1:3,iPartBound)    = (/0.,0.,0./)
  CASE('rot_periodic')
    GEO%RotPeriodicBC = .TRUE.
    nRotPeriodicBCs  = nRotPeriodicBCs + 1
    PartBound%TargetBoundCond(iPartBound)  = PartBound%RotPeriodicBC
    PartBound%RotPeriodicAngle(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotPeriodicAngle')
    IF(ALMOSTZERO(PartBound%RotPeriodicAngle(iPartBound))) THEN
      CALL abort(__STAMP__,'Angle for rotational periodicity can not be zero (using the right-hand rule!)')
    END IF
   ! Rotate the particle slightly inside the domain
    PartBound%RotPeriodicAngle(iPartBound) = PartBound%RotPeriodicAngle(iPartBound) / 180. * PI
  CASE DEFAULT
    SWRITE(*,*) ' Boundary does not exist: ', TRIM(tmpString)
    CALL abort(__STAMP__,'Particle Boundary Condition does not exist')
  END SELECT
  PartBound%SourceBoundName(iPartBound) = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-SourceName'))

  ! Surface particle output to .h5
  PartBound%BoundaryParticleOutputHDF5(iPartBound)      = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-BoundaryParticleOutput')
  IF(PartBound%BoundaryParticleOutputHDF5(iPartBound)) DoBoundaryParticleOutputHDF5=.TRUE.
END DO

! Check if there is an particle init with photon SEE
FoundPartBoundSEE=.FALSE.
SpecLoop: DO iSpec=1,nSpecies
  DO iInit=1, Species(iSpec)%NumberOfInits
    IF(StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_SEE'))THEN
      FoundPartBoundSEE = .TRUE.
      EXIT SpecLoop
    END IF ! StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_SEE')
  END DO
END DO SpecLoop

! Connect emission inits to particle boundaries for output
IF(DoBoundaryParticleOutputHDF5.OR.FoundPartBoundSEE) CALL InitializeVariablesSpeciesBoundary()

AdaptWallTemp = GETLOGICAL('Part-AdaptWallTemp')

IF(GEO%RotPeriodicBC) THEN
  GEO%RotPeriodicAxi   = GETINT('Part-RotPeriodicAxi')
! Check whether two corresponding RotPeriodic BCs are always set
  IF(MOD(nRotPeriodicBCs,2).NE.0) THEN
    CALL abort(__STAMP__,'ERROR: Uneven number of rot_periodic BCs. Check whether two corresponding RotPeriodic BCs are set!')
  ELSE IF(nRotPeriodicBCs.NE.2) THEN
    DO iPartBound=1,nPartBound
      WRITE(UNIT=hilf,FMT='(I0)') iPartBound
      IF(PartBound%TargetBoundCond(iPartBound).EQ.PartBound%RotPeriodicBC) THEN
        PartBound%RotPeriodicMin(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotPeriodicMin')
        PartBound%RotPeriodicMax(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotPeriodicMax')
        IF(PartBound%RotPeriodicMin(iPartBound).GE.PartBound%RotPeriodicMax(iPartBound)) THEN
          CALL abort(__STAMP__,'ERROR: Minimum coordinate at rotational axis is higher than maximum coordinate at rotational axis')
        END IF
      END IF
    END DO
  END IF
END IF

! Surface particle output to .h5
IF(DoBoundaryParticleOutputHDF5)THEN
  ! This array is not de-allocated during load balance as it is only written to .h5 during WriteStateToHDF5()
  IF(.NOT.ALLOCATED(PartStateBoundary))THEN
    ALLOCATE(PartStateBoundary(1:nVarPartStateBoundary,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary=0.
  END IF ! .NOT.ALLOCATED(PartStateBoundary)
END IF

! Set mapping from field boundary to particle boundary index and vice versa
ALLOCATE(PartBound%MapToPartBC(1:nBCs))
PartBound%MapToPartBC(:)=-10
ALLOCATE(PartBound%MapToFieldBC(1:nPartBound))
PartBound%MapToFieldBC=-1
DO iPBC=1,nPartBound
  DO iBC = 1, nBCs
    IF (BoundaryType(iBC,BC_TYPE).EQ.0) THEN
      PartBound%MapToPartBC(iBC) = -1 !there are no internal BCs in the mesh, they are just in the name list!
      LBWRITE(*,*)"... PartBound",iPBC,"is internal bound, no mapping needed"
    ELSEIF(BoundaryType(iBC,BC_TYPE).EQ.100)THEN
      IF(TrackingMethod.EQ.REFMAPPING)THEN
        SWRITE(UNIT_STDOUT,'(A)') ' Analyze sides are not implemented for RefMapping=T, because '//  &
                                  ' orientation of SideNormVec is unknown.'
        CALL abort(__STAMP__,' Analyze-BCs cannot be used for internal reflection in general cases! ')
      END IF
    END IF
    IF (TRIM(BoundaryName(iBC)).EQ.TRIM(PartBound%SourceBoundName(iPBC))) THEN
      PartBound%MapToPartBC(iBC) = iPBC !PartBound%TargetBoundCond(iPBC)
      PartBound%MapToFieldBC(iPBC) = iBC ! part BC to field BC
      LBWRITE(*,*)"... Mapped PartBound",iPBC,"on FieldBound", iBC,",i.e.:",TRIM(BoundaryName(iBC))
    END IF
  END DO
END DO
! Errorhandler for PartBound-Types that could not be mapped to the
! FieldBound-Types.
DO iBC = 1,nBCs
  IF (PartBound%MapToPartBC(iBC).EQ.-10) CALL abort(__STAMP__,' PartBound%MapToPartBC for Boundary is not set. iBC: :',iBC)
END DO

!-- Floating Potential
ALLOCATE(BCdata_auxSF(1:nPartBound))
DO iPartBound=1,nPartBound
  BCdata_auxSF(iPartBound)%SideNumber=-1 !init value when not used
  BCdata_auxSF(iPartBound)%GlobalArea=0.
  BCdata_auxSF(iPartBound)%LocalArea=0.
END DO

IF(ANY(PartBound%SurfaceModel.EQ.1)) THEN
  ! Open the species database
  CALL OpenDataFile(TRIM(SpeciesDatabase),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  ! Check if the correct dataset exists
  StickingCoefficientExists = .FALSE.
  dsetname = TRIM('/Surface-Chemistry/StickingCoefficient')
  CALL DatasetExists(File_ID,TRIM(dsetname),StickingCoefficientExists)
  IF(.NOT.StickingCoefficientExists) CALL abort(__STAMP__,'ERROR in InitializeVariablesPartBoundary: '//  &
                                      'No /Surface-Chemistry/StickingCoefficient dataset found in SpeciesDatabase!')
  ! Get dimensions
  CALL GetDataSize(File_ID,dsetname,nDims,HSize,attrib=.FALSE.)
  ! Allocate the data array
  ALLOCATE(StickingCoefficientData(INT(HSize(1),4),INT(HSize(2),4)))
  StickingCoefficientData = 0.
  ! Read-in array
  CALL ReadArray(TRIM(dsetname),2,INT(HSize,IK),0_IK,1,RealArray=StickingCoefficientData)
  CALL CloseDataFile()
END IF

END SUBROUTINE InitializeVariablesPartBoundary


!===================================================================================================================================
!> Read mapping for rotational periodicity from mesh file. If it does not yet exist, build the mapping and store in mesh.h5 for
!> faster initialization later on.
!===================================================================================================================================
SUBROUTINE InitParticleBoundaryRotPeriodic()
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input
USE MOD_HDF5_output       ,ONLY: WriteAttributeToHDF5
USE MOD_Mesh_Vars         ,ONLY: MeshFile!,RotPeriodicReBuild
USE MOD_Particle_MPI_Vars ,ONLY: halo_eps_velo
USE MOD_TimeDisc_Vars     ,ONLY: ManualTimeStep
USE MOD_ReadInTools       ,ONLY: GETLOGICAL
USE MOD_Globals_Vars      ,ONLY: ProjectName
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL           :: DatasetFound
CHARACTER(LEN=64) :: DatasetName,hilf
REAL              :: StartT,EndT
INTEGER           :: notMappedTotal
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT ROTATIONAL PERIODIC BOUNDARY...'
!RotPeriodicReBuild = GETLOGICAL('Part-RotPeriodicReBuild')
GETTIME(StartT)

DatasetFound = .FALSE.
notMappedTotal = 0
WRITE(UNIT=hilf,FMT='(ES10.4)') halo_eps_velo
DatasetName = 'RotPeriodicMap-v'//TRIM(hilf)
WRITE(UNIT=hilf,FMT='(ES10.4)') ManualTimeStep
DatasetName = TRIM(DatasetName)//'-dt'//TRIM(hilf)
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
CALL DatasetExists(File_ID,TRIM(DatasetName),DatasetFound)
CALL CloseDataFile()

IF(DatasetFound)THEN!.AND.(.NOT.RotPeriodicReBuild))THEN
  ! Read mapping from mesh file
  LBWRITE(Unit_StdOut,'(A)')" Reading ["//TRIM(DatasetName)//"] from mesh file ["//TRIM(MeshFile)//"]"
  ! This feature is not implemented yet. If more performance is required, e.g., during load balancing, store the mapping in the
  ! mesh file and simply read it during initialization without the need for searching the connections again
  !CALL ReadArray(TRIM(DatasetName),2,(/x,x,x/),0_IK,1,RealArray=mapping)
  CALL ReadAttribute(File_ID,'notMappedTotal',1,IntScalar=notMappedTotal)
  CALL CloseDataFile()
ELSE
  ! Create mapping and store in mesh file
  CALL BuildParticleBoundaryRotPeriodic(notMappedTotal)

  ! Store number of unmapped elements in mesh file
  IF(MPIRoot)THEN
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'notMappedTotal',1,IntegerScalar=notMappedTotal,Overwrite=.TRUE.)
    CALL CloseDataFile()
  END IF ! MPIRoot
END IF ! DatasetFound


IF(notMappedTotal.GT.0)THEN
  LBWRITE(Unit_StdOut,'(A,I0,A)') ' | Warning: Found ',notMappedTotal,' rot periodic sides that did not find a corresponding side.'
  LBWRITE(Unit_StdOut,'(A)')" | The halo region (halo flag 2) merely reaches a rot periodic BC side but not any further."
  LBWRITE(Unit_StdOut,'(A)')" | See ElemData container 'LostRotPeriodicSides' for more information on where sides were unmatched."
  LBWRITE(Unit_StdOut,'(A)')" | This information is written to "//TRIM(ProjectName)//"_LostRotPeriodicSides.h5 (only when CalcMeshInfo=T)"
  !LBWRITE(Unit_StdOut,'(A)')" | If the file does not exist, it can be re-created with RotPeriodicReBuild=T"
END IF ! notMappedTotal.GT.0

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'INIT ROTATIONAL PERIODIC BOUNDARY DONE!')

END SUBROUTINE InitParticleBoundaryRotPeriodic


!===================================================================================================================================
!> Build Mapping for rotational periodicity: RotPeriodicSide -> SideID2 (Side on corresponding BC).
!> In RotPeriodicBC (particle_boundary_condition.f90): SideID -> SurfSideID -> RotPeriodicSide
!>                                                     RotPeriodicSide -> SideID2
!> (1) counting rotational periodic sides and build mapping from SurfSideID -> RotPeriodicSide
!> (2) Build bounding boxes (in 2D reference system) for all nRotPeriodicSides
!> (3) find Side on corresponding BC and build mapping RotPeriodicSide -> SideID2 (and vice versa)
!>     counting potential rotational periodic sides (for not conform meshes)
!> (4) reallocate array due to number of potential rotational periodic sides
!===================================================================================================================================
SUBROUTINE BuildParticleBoundaryRotPeriodic(notMappedTotal)
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfTotalSides,SurfSide2GlobalSide,PartBound
USE MOD_Particle_Boundary_Vars  ,ONLY: RotPeriodicSideMapping, NumRotPeriodicNeigh, SurfSide2RotPeriodicSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared, NodeCoords_Shared, ElemSideNodeID_Shared, GEO, ElemInfo_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared, NodeInfo_Shared, NodeToElemInfo, NodeToElemMapping
USE MOD_Mesh_Vars               ,ONLY: LostRotPeriodicSides,nElems
USE MOD_Analyze_Vars            ,ONLY: CalcMeshInfo
USE MOD_IO_HDF5                 ,ONLY: AddToElemData,ElementOut
USE MOD_HDF5_Output_ElemData    ,ONLY: WriteLostRotPeriodicSidesToHDF5
#if USE_MPI
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2RotPeriodicSide_Shared,SurfSide2RotPeriodicSide_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: Rot2Glob_temp_Shared,Rot2Glob_temp_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: RotPeriodicSideMapping_temp_Shared,RotPeriodicSideMapping_temp_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: RotPeriodicSideMapping_Shared,RotPeriodicSideMapping_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: BoundingBox_Shared,BoundingBox_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: NumRotPeriodicNeigh_Shared,NumRotPeriodicNeigh_Shared_Win
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared
#else
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)  :: notMappedTotal
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iSide, jSide, nRotPeriodicSides, SideID,SideID2, MaxNumRotPeriodicNeigh, iNode, iNeigh, jNeigh
INTEGER                           :: NodeID, CNElemID, LocSideID, k, l, m, CNElemID2, LocSideID2, TestElemID, UniqueNodeID
INTEGER                           :: iElem, jElem, NewNeighNumber, kNeigh
LOGICAL                           :: mySide, FoundConnection
REAL                              :: iNodeVec(1:3), jNodeVec(1:3)
REAL                              :: iNodeR, iNodeH, jNodeR, jNodeH, Node2Rmin, Node2Rmax, Node2Hmin, Node2Hmax, dh, dr
INTEGER,PARAMETER                 :: NbrOfRotConnections=1000
INTEGER                           :: notMapped,GlobalElemID
INTEGER                           :: firstSide,lastSide,offsetSide
INTEGER,ALLOCPOINT,DIMENSION(:)   :: Rot2Glob_temp
INTEGER,ALLOCPOINT,DIMENSION(:,:) :: RotPeriodicSideMapping_temp
REAL,ALLOCPOINT,DIMENSION(:,:)    :: BoundingBox
!===================================================================================================================================

nRotPeriodicSides = 0
offsetSide        = 0

! Surf sides are shared, array calculation can be distributed
#if USE_MPI
CALL Allocate_Shared((/nComputeNodeSurfTotalSides/) , SurfSide2RotPeriodicSide_Shared_Win , SurfSide2RotPeriodicSide_Shared)
CALL MPI_WIN_LOCK_ALL(0 , SurfSide2RotPeriodicSide_Shared_Win , IERROR)
CALL Allocate_Shared((/nComputeNodeSurfTotalSides/) , Rot2Glob_temp_Shared_Win            , Rot2Glob_temp_Shared)
CALL MPI_WIN_LOCK_ALL(0 , Rot2Glob_temp_Shared_Win            , IERROR)
SurfSide2RotPeriodicSide => SurfSide2RotPeriodicSide_Shared
Rot2Glob_temp            => Rot2Glob_temp_Shared
IF (myComputeNodeRank.EQ.0) SurfSide2RotPeriodicSide = -1.
IF (myComputeNodeRank.EQ.0) Rot2Glob_temp            = -1.
CALL BARRIER_AND_SYNC(SurfSide2RotPeriodicSide_Shared_Win , MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Rot2Glob_temp_Shared_Win            , MPI_COMM_SHARED)
firstSide = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nSurfTotalSides
ALLOCATE(SurfSide2RotPeriodicSide(firstSide:lastSide))
SurfSide2RotPeriodicSide(:) = -1
ALLOCATE(Rot2Glob_temp(firstSide:lastSide))
#endif /*USE_MPI*/

! (1) Count rotational periodic sides and build mapping from SurfSideID -> RotPeriodicSide
#if USE_MPI
! Only when more than 1 core per node
IF(nComputeNodeProcessors.GT.1)THEN
  ! Loop once and get offset
  DO iSide=firstSide, lastSide
    SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
    IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%RotPeriodicBC) &
        nRotPeriodicSides = nRotPeriodicSides + 1
  END DO

  CALL MPI_EXSCAN(nRotPeriodicSides,offsetSide,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
  nRotPeriodicSides = offsetSide
END IF ! nComputeNodeProcessors.GT.1
#endif /*USE_MPI*/

DO iSide=firstSide, lastSide
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%RotPeriodicBC) THEN
    nRotPeriodicSides = nRotPeriodicSides + 1
    Rot2Glob_temp(nRotPeriodicSides) = SideID
    SurfSide2RotPeriodicSide(iSide) = nRotPeriodicSides ! Store RotSideID
  END IF
END DO

#if USE_MPI
! last proc knows CN total number of rot periodic sides
CALL MPI_BCAST(nRotPeriodicSides,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
CALL BARRIER_AND_SYNC(SurfSide2RotPeriodicSide_Shared_Win , MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Rot2Glob_temp_Shared_Win            , MPI_COMM_SHARED)
CALL Allocate_Shared((/nRotPeriodicSides/) , NumRotPeriodicNeigh_Shared_Win , NumRotPeriodicNeigh_Shared)
CALL MPI_WIN_LOCK_ALL(0 , NumRotPeriodicNeigh_Shared_Win , IERROR)
NumRotPeriodicNeigh => NumRotPeriodicNeigh_Shared
CALL Allocate_Shared((/nRotPeriodicSides,NbrOfRotConnections/) , RotPeriodicSideMapping_temp_Shared_Win , RotPeriodicSideMapping_temp_Shared)
CALL MPI_WIN_LOCK_ALL(0 , RotPeriodicSideMapping_temp_Shared_Win , IERROR)
RotPeriodicSideMapping_temp => RotPeriodicSideMapping_temp_Shared
IF (myComputeNodeRank.EQ.0) NumRotPeriodicNeigh = 0
IF (myComputeNodeRank.EQ.0) RotPeriodicSideMapping_temp = 0
CALL BARRIER_AND_SYNC(NumRotPeriodicNeigh_Shared_Win         , MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(RotPeriodicSideMapping_temp_Shared_Win , MPI_COMM_SHARED)
firstSide = INT(REAL( myComputeNodeRank   )*REAL(nRotPeriodicSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nRotPeriodicSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nRotPeriodicSides
ALLOCATE(NumRotPeriodicNeigh(nRotPeriodicSides))
NumRotPeriodicNeigh = 0
! number of potential rotational periodic sides is unknown => allocate mapping array with fixed number of NbrOfRotConnections
! and reallocate at the end of subroutine
ALLOCATE(RotPeriodicSideMapping_temp(nRotPeriodicSides,NbrOfRotConnections))
RotPeriodicSideMapping_temp = 0
#endif /*USE_MPI*/

notMapped=0
MaxNumRotPeriodicNeigh = 0
! Defining rotation matrix
SELECT CASE(GEO%RotPeriodicAxi)
  CASE(1) ! x-rotation axis
    k = 1
    l = 2
    m = 3
  CASE(2) ! y-rotation axis
    k = 2
    l = 3
    m = 1
  CASE(3) ! z-rotation axis
    k = 3
    l = 1
    m = 2
END SELECT

! (2) Build bounding boxes (in 2D reference system) for all nRotPeriodicSides
#if USE_MPI
CALL Allocate_Shared((/4,nRotPeriodicSides/) , BoundingBox_Shared_Win , BoundingBox_Shared)
CALL MPI_WIN_LOCK_ALL(0 , BoundingBox_Shared_Win , IERROR)
BoundingBox => BoundingBox_Shared
CALL BARRIER_AND_SYNC(BoundingBox_Shared_Win , MPI_COMM_SHARED)
#else
ALLOCATE(BoundingBox(4,nRotPeriodicSides))
#endif /*USE_MPI*/

DO iSide = firstSide, lastSide

  ! Get side information
  SideID    = Rot2Glob_temp(iSide)
  CNElemID  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
  LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)

  ! Loop over all 4 nodes
  DO iNode=1, 4
    NodeID        = ElemSideNodeID_Shared(iNode,LocSideID,CNElemID) + 1
    jNodeVec(1:3) = NodeCoords_Shared(1:3,NodeID)
    jNodeH        = jNodeVec(k)
    jNodeR        = SQRT(jNodeVec(l)**2+jNodeVec(m)**2)
    IF(iNode.EQ. 1) THEN
      Node2Hmin = jNodeH
      Node2Hmax = jNodeH
      Node2Rmin = jNodeR
      Node2Rmax = jNodeR
    ELSE
      Node2Hmin = MIN(Node2Hmin,jNodeH)
      Node2Hmax = MAX(Node2Hmax,jNodeH)
      Node2Rmin = MIN(Node2Rmin,jNodeR)
      Node2Rmax = MAX(Node2Rmax,jNodeR)
    END IF
  END DO
  ! Add tolerance by increasing the bounding box size by 1 percent of the length in h and r
  dh = Node2Hmax-Node2Hmin
  dr = Node2Rmax-Node2Rmin
  BoundingBox(1,iSide) = Node2Hmin-dh*0.01
  BoundingBox(2,iSide) = Node2Hmax+dh*0.01
  BoundingBox(3,iSide) = Node2Rmin-dr*0.01
  BoundingBox(4,iSide) = Node2Rmax+dr*0.01

END DO ! iSide = firstSide, lastSide
#if USE_MPI
CALL BARRIER_AND_SYNC(BoundingBox_Shared_Win , MPI_COMM_SHARED)
#endif /*USE_MPI*/

! (3) find Side on corresponding BC and build mapping RotPeriodicSide -> SideID2 (and vice versa)
!     counting potential rotational periodic sides (for non-conforming meshes)
! Use named loops: Loop over the assigned iSides and compare against all nRotPeriodicSides
iSideLoop: DO iSide = firstSide, lastSide
  SideID    = Rot2Glob_temp(iSide)
  CNElemID  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
  LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)

  jSideLoop: DO jSide = 1, nRotPeriodicSides
    SideID2 = Rot2Glob_temp(jSide)
    FoundConnection = .FALSE.

    ! Check if both sides are on the same boundary, i.e., they cannot be connected
    IF(PartBound%RotPeriodicAngle(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ. &
        PartBound%RotPeriodicAngle(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID2)))) CYCLE jSideLoop

    ! Check if jSide is assigned to the proc
    mySide = (jSide.GE.firstSide).AND.(jSide.LE.lastSide)

    ! Check if the side was added in a previous step
    IF(mySide)THEN
      DO iNeigh = 1, NumRotPeriodicNeigh(jSide)
        IF(RotPeriodicSideMapping_temp(jSide,iNeigh).EQ.SideID) CYCLE jSideLoop
      END DO ! iNeigh = 1, NumRotPeriodicNeigh(jSide)
    END IF ! mySide

    ! Loop the 4 nodes of iSide and test against the bounding box of jSide
    iNodeLoop: DO iNode = 1, 4
      NodeID = ElemSideNodeID_Shared(iNode,LocSideID,CNElemID) + 1
      ! Calculate node coordinates in reference system
      iNodeVec(1:3) = NodeCoords_Shared(1:3,NodeID)
      iNodeH        = iNodeVec(k)
      iNodeR        = SQRT(iNodeVec(l)**2+iNodeVec(m)**2)

      ! Cycle if outside of bounding box of jSide
      IF(BoundingBox(1,jSide).GT.iNodeH) CYCLE iNodeLoop
      IF(BoundingBox(2,jSide).LT.iNodeH) CYCLE iNodeLoop
      IF(BoundingBox(3,jSide).GT.iNodeR) CYCLE iNodeLoop
      IF(BoundingBox(4,jSide).LT.iNodeR) CYCLE iNodeLoop

      FoundConnection = .TRUE.

      ! Check if the side was added in a previous step
      DO iNeigh = 1, NumRotPeriodicNeigh(iSide)
        IF(RotPeriodicSideMapping_temp(iSide,iNeigh).EQ.SideID2) CYCLE iNodeLoop
      END DO ! iNeigh = 1, NumRotPeriodicNeigh(iSide)

      ! Found connection
      NumRotPeriodicNeigh(iSide) = NumRotPeriodicNeigh(iSide) + 1
      IF(NumRotPeriodicNeigh(iSide).GT.NbrOfRotConnections) CALL abort(__STAMP__,&
            ' ERROR: Number of rotational periodic side exceed fixed number of ',IntInfoOpt=NbrOfRotConnections)
      RotPeriodicSideMapping_temp(iSide,NumRotPeriodicNeigh(iSide)) = SideID2

      ! Only do vice versa mapping if the processor has been assigned this side
      IF(mySide)THEN
        ! Check if the side was added in a previous step
        DO iNeigh = 1, NumRotPeriodicNeigh(jSide)
          IF(RotPeriodicSideMapping_temp(jSide,iNeigh).EQ.SideID) CYCLE iNodeLoop
        END DO ! iNeigh = 1, NumRotPeriodicNeigh(iSide)
        NumRotPeriodicNeigh(jSide) = NumRotPeriodicNeigh(jSide) + 1
        IF(NumRotPeriodicNeigh(jSide).GT.NbrOfRotConnections) CALL abort(__STAMP__,&
            ' jSide: Number of rotational periodic side exceed fixed number of ',IntInfoOpt=NbrOfRotConnections)
        RotPeriodicSideMapping_temp(jSide,NumRotPeriodicNeigh(jSide)) = SideID
      END IF ! mySide

      ! Exit loop over nodes when side has been assigned
      EXIT iNodeLoop

    END DO iNodeLoop ! iNode = 1, 4

    IF(.NOT.FoundConnection) THEN
      ! Double check is needed if bounding box of jSide is within the original side (iSide) -> check whether jSide is within
      ! bounding box of iSide. Need to get the NodeIDs of jSide (= SideID2).
      CNElemID2  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID2))
      LocSideID2 = SideInfo_Shared(SIDE_LOCALID,SideID2)

      ! Loop the 4 nodes of jSide and test against the bounding box of iSide
      jNodeLoop: DO iNode = 1, 4
        NodeID = ElemSideNodeID_Shared(iNode,LocSideID2,CNElemID2) + 1
        ! Calculate node coordinates in reference system
        iNodeVec(1:3) = NodeCoords_Shared(1:3,NodeID)
        iNodeH        = iNodeVec(k)
        iNodeR        = SQRT(iNodeVec(l)**2+iNodeVec(m)**2)

        ! Cycle if outside of bounding box of iSide
        IF(BoundingBox(1,iSide).GT.iNodeH) CYCLE jNodeLoop
        IF(BoundingBox(2,iSide).LT.iNodeH) CYCLE jNodeLoop
        IF(BoundingBox(3,iSide).GT.iNodeR) CYCLE jNodeLoop
        IF(BoundingBox(4,iSide).LT.iNodeR) CYCLE jNodeLoop

        ! Check if the side was added in a previous step
        DO iNeigh = 1, NumRotPeriodicNeigh(iSide)
          IF(RotPeriodicSideMapping_temp(iSide,iNeigh).EQ.SideID2) CYCLE jNodeLoop
        END DO ! iNeigh = 1, NumRotPeriodicNeigh(iSide)

        ! Found connection
        NumRotPeriodicNeigh(iSide) = NumRotPeriodicNeigh(iSide) + 1
        IF(NumRotPeriodicNeigh(iSide).GT.NbrOfRotConnections) CALL abort(__STAMP__,&
              ' ERROR: Number of rotational periodic side exceed fixed number of ',IntInfoOpt=NbrOfRotConnections)
        RotPeriodicSideMapping_temp(iSide,NumRotPeriodicNeigh(iSide)) = SideID2

        ! Only do vice versa mapping if the processor has been assigned this side
        IF(mySide)THEN
          ! Check if the side was added in a previous step
          DO iNeigh = 1, NumRotPeriodicNeigh(jSide)
            IF(RotPeriodicSideMapping_temp(jSide,iNeigh).EQ.SideID) CYCLE jNodeLoop
          END DO ! iNeigh = 1, NumRotPeriodicNeigh(iSide)
          NumRotPeriodicNeigh(jSide) = NumRotPeriodicNeigh(jSide) + 1
          IF(NumRotPeriodicNeigh(jSide).GT.NbrOfRotConnections) CALL abort(__STAMP__,&
              ' jSide: Number of rotational periodic side exceed fixed number of ',IntInfoOpt=NbrOfRotConnections)
          RotPeriodicSideMapping_temp(jSide,NumRotPeriodicNeigh(jSide)) = SideID
        END IF ! mySide

        ! Exit loop over nodes when side has been assigned
        EXIT jNodeLoop

      END DO jNodeLoop ! iNode = 1, 4
    END IF ! .NOT.FoundConnection

  END DO jSideLoop ! jSide = 1, nRotPeriodicSides

  ! Check if iSide could not be mapped to any other side. There should be at least one
  IF(NumRotPeriodicNeigh(iSide).EQ.0) THEN
    IF(ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,SideID)).NE.3) THEN
      ! Count number of sides that could not be mapped (warning output + info in h5 file when CalcMeshInfo=T)
      notMapped = notMapped + 1
    ELSE
      ! Found side on element that is a neighbor element in rot halo region (they have halo flag 3)
      NumRotPeriodicNeigh(iSide) = 1
      RotPeriodicSideMapping_temp(iSide,NumRotPeriodicNeigh(iSide)) = 0
    END IF ! ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,SideID)).NE.3
  END IF ! NumRotPeriodicNeigh(iSide).EQ.0

END DO iSideLoop ! iSide = firstSide, lastSide

! Addition of neighbour elements for each node of a mapped side (especially a problem for tetrahedron-based meshes)
DO iSide = firstSide, lastSide
  NewNeighNumber = NumRotPeriodicNeigh(iSide)
  DO iNeigh=1, NumRotPeriodicNeigh(iSide)
    SideID = RotPeriodicSideMapping_temp(iSide,iNeigh)
    CNElemID  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
    LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
    kNodeLoop: DO iNode = 1, 4
      NodeID = ElemSideNodeID_Shared(iNode,LocSideID,CNElemID) + 1
      UniqueNodeID = NodeInfo_Shared(NodeID)
      ElemLoop: DO iElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
        TestElemID = NodeToElemInfo(iElem)
        ! Check if its the same element
        IF(CNElemID.EQ.TestElemID) CYCLE ElemLoop
        ! Check if element is already in the list of the OLD neighbours
        NeighLoop: DO jNeigh=1, NumRotPeriodicNeigh(iSide)
          ! Skip yourself
          IF(iNeigh.EQ.jNeigh) CYCLE NeighLoop
          CNElemID2  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,RotPeriodicSideMapping_temp(iSide,jNeigh)))
          IF(CNElemID2.EQ.TestElemID) CYCLE ElemLoop
        END DO NeighLoop
        ! Check if element is already in the list of the NEW neighbours
        NewNeighLoop: DO kNeigh = NumRotPeriodicNeigh(iSide),NewNeighNumber
          jElem = ABS(RotPeriodicSideMapping_temp(iSide,kNeigh))
          IF(jElem.EQ.GetGlobalElemID(TestElemID)) CYCLE ElemLoop
        END DO NewNeighLoop
        ! Add element to the neighbour list
        NewNeighNumber = NewNeighNumber + 1
        IF(NewNeighNumber.GT.NbrOfRotConnections) CALL abort(__STAMP__,&
            ' NewNeighNumber: Number of rotational periodic side exceed fixed number of ',IntInfoOpt=NbrOfRotConnections)
        RotPeriodicSideMapping_temp(iSide,NewNeighNumber) = -GetGlobalElemID(TestElemID)
      END DO ElemLoop
    END DO kNodeLoop
  END DO
  NumRotPeriodicNeigh(iSide) = NewNeighNumber
END DO

! (4) reallocate array due to number of potential rotational periodic sides
#if USE_MPI
CALL BARRIER_AND_SYNC(NumRotPeriodicNeigh_Shared_Win, MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(RotPeriodicSideMapping_temp_Shared_Win, MPI_COMM_SHARED)
! The allreduce is only required when a global array for writing to .h5 is to be used
!CALL MPI_ALLREDUCE(MAXVAL(NumRotPeriodicNeigh) , MaxNumRotPeriodicNeigh , 1 , MPI_INTEGER , MPI_MAX , MPI_COMM_WORLD , iError)
#endif /*USE_MPI*/
MaxNumRotPeriodicNeigh = MAXVAL(NumRotPeriodicNeigh)

#if USE_MPI
CALL Allocate_Shared((/nRotPeriodicSides,MaxNumRotPeriodicNeigh/) , RotPeriodicSideMapping_Shared_Win , RotPeriodicSideMapping_Shared)
CALL MPI_WIN_LOCK_ALL(0 , RotPeriodicSideMapping_Shared_Win , IERROR)
RotPeriodicSideMapping => RotPeriodicSideMapping_Shared
IF (myComputeNodeRank.EQ.0) RotPeriodicSideMapping = -1
CALL BARRIER_AND_SYNC(RotPeriodicSideMapping_Shared_Win , MPI_COMM_SHARED)
#else
ALLOCATE(RotPeriodicSideMapping(nRotPeriodicSides,MaxNumRotPeriodicNeigh))
RotPeriodicSideMapping = -1
#endif /*USE_MPI*/

DO iSide=1, nRotPeriodicSides
  DO iNeigh=1, MaxNumRotPeriodicNeigh
    SideID = RotPeriodicSideMapping_temp(iSide,iNeigh)
    IF(SideID.GT.0)THEN
      GlobalElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      RotPeriodicSideMapping(iSide,iNeigh) = GlobalElemID
    ELSE IF(SideID.LT.0) THEN
      RotPeriodicSideMapping(iSide,iNeigh) = ABS(SideID)
    END IF ! SideID.NE.0
  END DO
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(notMapped , notMappedTotal , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_WORLD , IERROR)
#else
notMappedTotal = notMapped
#endif /*USE_MPI*/

IF(notMappedTotal.GT.0)THEN
  IF(CalcMeshInfo)THEN
    ALLOCATE(LostRotPeriodicSides(1:nElems))
    LostRotPeriodicSides=notMapped
    CALL AddToElemData(ElementOut,'LostRotPeriodicSides',LongIntArray=LostRotPeriodicSides)
    CALL WriteLostRotPeriodicSidesToHDF5()
  END IF ! CalcMeshInfo
  !IF(MPIroot) CALL abort(__STAMP__,'At least one rot periodic side did not find a corresponding side.')
END IF ! notMappedTotal.GT.0

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
CALL UNLOCK_AND_FREE(Rot2Glob_temp_Shared_Win)
ADEALLOCATE(Rot2Glob_temp_Shared)
ADEALLOCATE(Rot2Glob_temp)
CALL UNLOCK_AND_FREE(RotPeriodicSideMapping_temp_Shared_Win)
ADEALLOCATE(RotPeriodicSideMapping_temp_Shared)
ADEALLOCATE(RotPeriodicSideMapping_temp)
CALL UNLOCK_AND_FREE(BoundingBox_Shared_Win)
ADEALLOCATE(BoundingBox_Shared)
ADEALLOCATE(BoundingBox)
#else
DEALLOCATE(Rot2Glob_temp)
DEALLOCATE(RotPeriodicSideMapping_temp)
DEALLOCATE(BoundingBox)
#endif /*USE_MPI*/

END SUBROUTINE BuildParticleBoundaryRotPeriodic


!===================================================================================================================================
!> Allocate shared array for the side-local wall temperature: BoundaryWallTemp_Shared and initialize with read-in wall temperature
!===================================================================================================================================
SUBROUTINE InitAdaptiveWallTemp()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, nComputeNodeSurfTotalSides, BoundaryWallTemp, SurfSide2GlobalSide,nSurfSample
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared, NodeCoords_Shared, ElemSideNodeID_Shared
USE MOD_SurfaceModel_Tools        ,ONLY: CalcWallTempGradient
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars           ,ONLY: MPI_COMM_SHARED, myComputeNodeRank, nComputeNodeProcessors
USE MOD_Particle_Boundary_Vars    ,ONLY: BoundaryWallTemp_Shared, BoundaryWallTemp_Shared_Win
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: firstSide, lastSide, iSide, SideID, iBC, ElemID, CNElemID, LocSideID, iNode
REAL                              :: SideMidpoint(1:3)
!===================================================================================================================================
IF (.NOT.PartBound%OutputWallTemp) RETURN

#if USE_MPI
!> Then shared arrays for boundary sampling
CALL Allocate_Shared((/nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),BoundaryWallTemp_Shared_Win,BoundaryWallTemp_Shared)
CALL MPI_WIN_LOCK_ALL(0,BoundaryWallTemp_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  BoundaryWallTemp_Shared = 0.
END IF
BoundaryWallTemp => BoundaryWallTemp_Shared
CALL BARRIER_AND_SYNC(BoundaryWallTemp_Shared_Win,MPI_COMM_SHARED)
firstSide = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(BoundaryWallTemp(nSurfSample,nSurfSample,1:nComputeNodeSurfTotalSides))
BoundaryWallTemp = 0.
firstSide = 1
lastSide  = nComputeNodeSurfTotalSides
#endif /*USE_MPI*/

DO iSide = firstSide,LastSide
  ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  IF (PartBound%MomentumACC(iBC).GT.0.0) THEN
    IF(PartBound%WallTemp2(iBC).GT.0.) THEN
      ElemID    = SideInfo_Shared(SIDE_ELEMID ,SideID)
      CNElemID  = GetCNElemID(ElemID)
      LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
      SideMidpoint(1:3) = 0.
      DO iNode = 1,4
        SideMidpoint(1:3) = SideMidpoint(1:3) + NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode,LocSideID,CNElemID)+1) / 4
      END DO
      BoundaryWallTemp(:,:,iSide) = CalcWallTempGradient(SideMidpoint,iBC)
    ELSE
      BoundaryWallTemp(:,:,iSide) = PartBound%WallTemp(iBC)
    END IF
  END IF
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(BoundaryWallTemp_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE InitAdaptiveWallTemp


SUBROUTINE FinalizeParticleBoundary()
!----------------------------------------------------------------------------------------------------------------------------------!
!> Finalize particle boundary variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars
#if USE_MPI
USE MOD_MPI_Shared_vars    ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%MomentumACC)
SDEALLOCATE(PartBound%OnlySpecular)
SDEALLOCATE(PartBound%OnlyDiffuse)
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%WallTemp2)
SDEALLOCATE(PartBound%WallTempDelta)
SDEALLOCATE(PartBound%TempGradStart)
SDEALLOCATE(PartBound%TempGradEnd)
SDEALLOCATE(PartBound%TempGradVec)
SDEALLOCATE(PartBound%TempGradDir)
SDEALLOCATE(PartBound%TransACC)
SDEALLOCATE(PartBound%VibACC)
SDEALLOCATE(PartBound%RotACC)
SDEALLOCATE(PartBound%ElecACC)
SDEALLOCATE(PartBound%Resample)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%RotVelo)
SDEALLOCATE(PartBound%RotOmega)
SDEALLOCATE(PartBound%RotPeriodicAngle)
SDEALLOCATE(PartBound%RotPeriodicMin)
SDEALLOCATE(PartBound%RotPeriodicMax)
SDEALLOCATE(PartBound%NbrOfSpeciesSwaps)
SDEALLOCATE(PartBound%ProbOfSpeciesSwaps)
SDEALLOCATE(PartBound%SpeciesSwaps)
SDEALLOCATE(PartBound%MapToPartBC)
SDEALLOCATE(PartBound%MapToFieldBC)
SDEALLOCATE(PartBound%SurfaceModel)
SDEALLOCATE(PartBound%Reactive)
SDEALLOCATE(PartBound%Dielectric)
SDEALLOCATE(PartBound%BoundaryParticleOutputHDF5)
SDEALLOCATE(PartBound%RadiativeEmissivity)

! Rotational periodic boundary condition
IF(GEO%RotPeriodicBC)THEN
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  CALL UNLOCK_AND_FREE(SurfSide2RotPeriodicSide_Shared_Win)
  CALL UNLOCK_AND_FREE(NumRotPeriodicNeigh_Shared_Win)
  ADEALLOCATE(SurfSide2RotPeriodicSide_Shared)
  ADEALLOCATE(SurfSide2RotPeriodicSide)
  ADEALLOCATE(NumRotPeriodicNeigh_Shared)
  ADEALLOCATE(NumRotPeriodicNeigh)
  CALL UNLOCK_AND_FREE(RotPeriodicSideMapping_Shared_Win)
  ADEALLOCATE(RotPeriodicSideMapping_Shared)
  ADEALLOCATE(RotPeriodicSideMapping)
#else
  SDEALLOCATE(SurfSide2RotPeriodicSide)
  SDEALLOCATE(NumRotPeriodicNeigh)
  SDEALLOCATE(RotPeriodicSideMapping)
#endif
END IF ! GEO%RotPeriodicBC

! Adaptive wall temperature (e.g. calculate from sampled heat flux)
IF (PartBound%OutputWallTemp) THEN
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  CALL UNLOCK_AND_FREE(BoundaryWallTemp_Shared_Win)
  ADEALLOCATE(BoundaryWallTemp_Shared)
  ADEALLOCATE(BoundaryWallTemp)
#else
  SDEALLOCATE(BoundaryWallTemp)
#endif
END IF
SDEALLOCATE(PartBound%UseAdaptedWallTemp)

! Do not deallocate during load balance (either communicate to new processor or simply keep on current processor)
#if USE_LOADBALANCE
IF(.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))THEN
#endif /*USE_LOADBALANCE*/
  SDEALLOCATE(PartStateBoundary)
#if USE_LOADBALANCE
END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/

END SUBROUTINE FinalizeParticleBoundary

END MODULE MOD_Particle_Boundary_Init
