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
PUBLIC :: DefineParametersParticleBoundary, InitializeVariablesPartBoundary, InitializeVariablesAuxBC, FinalizeParticleBoundary
PUBLIC :: InitParticleBoundaryRotPeriodic, InitAdaptiveWallTemp
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle boundaries and AuxBC
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBoundary()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Boundaries")

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
                                  'PB-WallVelo,Voltage,SpeciesSwaps.If condition=periodic:Part-nPeriodicVectors,'//&
                                  'Part-PeriodicVector[$]', 'open', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Boundary[$]-Dielectric' , 'Define if particle boundary [$] is a '//&
                              'dielectric interface, i.e. an interface between a dielectric and a non-dielectric or a between two'//&
                              ' different dielectrics [.TRUE.] or not [.FALSE.] (requires reflective BC and species swap for nSpecies)'&
                              , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Boundary[$]-BoundaryParticleOutput' , 'Define if the properties of particles impacting on '//&
                              'boundary [$] are to be stored in an additional .h5 file for post-processing analysis [.TRUE.] '//&
                              'or not [.FALSE.].'&
                              , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-Voltage', &
                                  'Deprecated parameter, do not utilize!', '1.0e20', numberedmulti=.TRUE.)
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
                                  ' Particles will be accelerated additionaly to the boundary interaction'//&
                                  ' through the rotating wall depending on their POI, rotation frequency and rotation axis.'//&
                                  ' In that case Part-Boundary[$]-WallVelo will be overwritten.' &
                                , '.FALSE.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-RotFreq'  &
                                , 'Rotation frequency of the wall in [Hz].' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-RotOrg'  &
                                , 'Origin of rotation axis (global x,y,z).' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-RotAxi'  &
                                , 'Direction of rotation axis (global x,y,z). Note: Rotation direction based on Right-hand rule!' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Boundary[$]-RotPeriodicDir'  &
                                , 'Angular degree of rotation periodicity in [deg].' &
                                , '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-WallTemp2'  &
                                , 'Second wall temperature (in [K]) of reflective particle boundary for a temperature gradient.' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-TemperatureGradientStart'  &
                                , 'Impose a temperature gradient by supplying a start/end vector and a second wall temperature.' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-TemperatureGradientEnd'  &
                                , 'Impose a temperature gradient by supplying a start/end vector and a second wall temperature.' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Boundary[$]-SurfaceModel'  &
                                , 'Defining surface to be treated reactively by defining Model used '//&
                                'for particle surface interaction. If any >0 then look in section SurfaceModel.\n'//&
                                '0: Maxwell scattering\n'//&
                                '5: SEE-E and SEE-I (secondary e- emission due to e- or i+ bombardment) '//&
                                    'by Levko2015 for copper electrodes\n'//&
                                '6: SEE-E (secondary e- emission due to e- bombardment) '//&
                                    'by Pagonakis2016 for molybdenum (originally from Harrower1956)'//&
                                '7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for '//&
                                'secondary e- emission with 0.13 probability) by Depla2009\n' &
                                , '0', numberedmulti=.TRUE.)
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
CALL prms%CreateIntOption(      'Part-nAuxBCs'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Number of auxillary BCs that are checked during tracing',  '0')
CALL prms%CreateIntOption(      'Part-AuxBC[$]-NbrOfSpeciesSwaps'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Number of Species to be changed at wall.',  '0', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-AuxBC[$]-Condition'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Used auxillary boundary condition for boundary[$].'//&
                                  '- open'//&
                                  '- reflective'//&
                                  '- periodic)'//&
                                  '-> more details see also Part-Boundary[$]-Condition',  'open', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-MomentumACC'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Momentum accommodation',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-WallTemp'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Wall temperature of boundary[$]',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-TransACC'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Translation accommodation on boundary [$]',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-VibACC'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Vibrational accommodation on boundary [$]',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-RotACC'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Rotational accommodation on boundary [$]',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-ElecACC'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Electronic accommodation on boundary [$]',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-AuxBC[$]-Resample'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Resample Equilibirum Distribution with reflection',  '.FALSE.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-AuxBC[$]-WallVelo'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Emitted velocity on boundary [$]', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-ProbOfSpeciesSwaps'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Probability of SpeciesSwaps at wall',  '1.', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Part-AuxBC[$]-SpeciesSwaps[$]'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Species to be changed at wall (out=: delete)', '0 , 0'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-AuxBC[$]-Type'  &
                                , 'TODO-DEFINE-PARAMETER'//&
                                  'Type of BC (plane, ...)',  'plane', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-AuxBC[$]-r_vec'  &
                                , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-radius'  &
                                , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.) !def. might be calculated!!!
CALL prms%CreateRealArrayOption('Part-AuxBC[$]-n_vec'  &
                                , 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-AuxBC[$]-axis'  &
                                , 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-lmin'  &
                                , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.) !def. might be calculated!!!
CALL prms%CreateRealOption(     'Part-AuxBC[$]-lmax'  &
                                , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.) !def. is calculated!!!
CALL prms%CreateLogicalOption(  'Part-AuxBC[$]-inwards'  &
                                , 'TODO-DEFINE-PARAMETER',  '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-rmax'  &
                                , 'TODO-DEFINE-PARAMETER',  '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-halfangle'  &
                                , 'TODO-DEFINE-PARAMETER',  '45.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-AuxBC[$]-zfac'  &
                                , 'TODO-DEFINE-PARAMETER',  '1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-AdaptWallTemp','Perform wall temperature adaptation at every macroscopic output.', '.TRUE.')

END SUBROUTINE DefineParametersParticleBoundary

SUBROUTINE InitializeVariablesPartBoundary()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_ReadInTools
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas
USE MOD_Mesh_Vars              ,ONLY: BoundaryName,BoundaryType, nBCs
USE MOD_Particle_Vars
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound,DoBoundaryParticleOutputHDF5,PartStateBoundary, AdaptWallTemp
USE MOD_Particle_Boundary_Vars ,ONLY: nVarPartStateBoundary
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPartBound, iBC, iPBC, iSwaps, MaxNbrOfSpeciesSwaps
INTEGER               :: ALLOCSTAT, dummy_int
CHARACTER(32)         :: hilf , hilf2
CHARACTER(200)        :: tmpString
LOGICAL               :: DeprecatedVoltage
!===================================================================================================================================
! Read in boundary parameters
dummy_int = CountOption('Part-nBounds')       ! check if Part-nBounds is present in .ini file
nPartBound = GETINT('Part-nBounds','1.') ! get number of particle boundaries
! Read-in number of porous boundaries
nPorousBC = GETINT('Surf-nPorousBC', '0')
IF ((nPartBound.LE.0).OR.(dummy_int.LT.0)) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nPartBound .LE. 0:', nPartBound)
END IF

ALLOCATE(PartBound%SourceBoundName(  1:nPartBound))
PartBound%SourceBoundName = ''
ALLOCATE(PartBound%TargetBoundCond(  1:nPartBound))
PartBound%TargetBoundCond = -1
ALLOCATE(PartBound%MomentumACC(      1:nPartBound))
PartBound%MomentumACC = -1
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
ALLOCATE(PartBound%RotFreq(          1:nPartBound))
PartBound%RotFreq = -1.
ALLOCATE(PartBound%RotOrg(       1:3,1:nPartBound))
PartBound%RotOrg = 0.
ALLOCATE(PartBound%RotAxi(       1:3,1:nPartBound))
PartBound%RotAxi = 0.
ALLOCATE(PartBound%RotPeriodicDir(  1:nPartBound))
PartBound%RotPeriodicDir = 0
ALLOCATE(PartBound%TempGradStart(1:3,1:nPartBound))
PartBound%TempGradStart = 0.
ALLOCATE(PartBound%TempGradEnd(  1:3,1:nPartBound))
PartBound%TempGradEnd = 0.
ALLOCATE(PartBound%TempGradVec(  1:3,1:nPartBound))
PartBound%TempGradVec = 0.
ALLOCATE(PartBound%SurfaceModel(     1:nPartBound))
PartBound%SurfaceModel = 0
ALLOCATE(PartBound%Reactive(         1:nPartBound))
PartBound%Reactive = .FALSE.
ALLOCATE(PartBound%Voltage(1:nPartBound))
PartBound%Voltage = 0.
DeprecatedVoltage = .FALSE.
ALLOCATE(PartBound%NbrOfSpeciesSwaps(1:nPartBound))
PartBound%NbrOfSpeciesSwaps = 0
ALLOCATE(PartBound%UseAdaptedWallTemp(1:nPartBound))
PartBound%UseAdaptedWallTemp = .FALSE.
ALLOCATE(PartBound%RadiativeEmissivity(1:nPartBound))
PartBound%RadiativeEmissivity = 1.

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
PartBound%Dielectric=.FALSE.
DoDielectricSurfaceCharge=.FALSE.
! Surface particle output to .h5
ALLOCATE(PartBound%BoundaryParticleOutputHDF5(1:nPartBound))
PartBound%BoundaryParticleOutputHDF5=.FALSE.
DoBoundaryParticleOutputHDF5=.FALSE.

PartMeshHasPeriodicBCs=.FALSE.
GEO%RotPeriodicBC =.FALSE.
#if defined(IMPA) || defined(ROS)
PartMeshHasReflectiveBCs=.FALSE.
#endif
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I0)') iPartBound
  tmpString = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-Condition','open'))

  SELECT CASE (TRIM(tmpString))
  CASE('open')
    PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
    PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage')
    IF(ABS(PartBound%Voltage(iPartBound)).LT.1e20) DeprecatedVoltage=.TRUE.
  CASE('reflective')
#if defined(IMPA) || defined(ROS)
    PartMeshHasReflectiveBCs=.TRUE.
#endif
    PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
    PartBound%MomentumACC(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-MomentumACC')
    PartBound%WallTemp(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp')
    PartBound%TransACC(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-TransACC')
    PartBound%VibACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-VibACC')
    PartBound%RotACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotACC')
    PartBound%ElecACC(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-ElecACC')
    PartBound%Resample(iPartBound)        = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Resample')
    PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-WallVelo',3)
    PartBound%RotVelo(iPartBound)         = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-RotVelo')
    IF(PartBound%RotVelo(iPartBound)) THEN
      PartBound%RotFreq(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotFreq')
      PartBound%RotOrg(1:3,iPartBound)      = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-RotOrg',3)
      PartBound%RotAxi(1:3,iPartBound)      = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-RotAxi',3)
    END IF
    PartBound%UseAdaptedWallTemp(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-UseAdaptedWallTemp')
    PartBound%RadiativeEmissivity(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-RadiativeEmissivity')
    PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage')
    IF(ABS(PartBound%Voltage(iPartBound)).LT.1e20) DeprecatedVoltage=.TRUE.
    PartBound%SurfaceModel(iPartBound)    = GETINT('Part-Boundary'//TRIM(hilf)//'-SurfaceModel')
    PartBound%WallTemp2(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp2')
    IF(PartBound%WallTemp2(iPartBound).GT.0.) THEN
      PartBound%TempGradStart(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-TemperatureGradientStart',3)
      PartBound%TempGradEnd(1:3,iPartBound)   = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-TemperatureGradientEnd',3)
      PartBound%WallTempDelta(iPartBound)   = PartBound%WallTemp2(iPartBound) - PartBound%WallTemp(iPartBound)
      PartBound%TempGradVec(1:3,iPartBound) = PartBound%TempGradEnd(1:3,iPartBound) - PartBound%TempGradStart(1:3,iPartBound)
    END IF
    ! check for correct surfacemodel input
    IF (PartBound%SurfaceModel(iPartBound).GT.0)THEN
      IF (.NOT.useDSMC) CALL abort(&
          __STAMP__&
          ,'Cannot use surfacemodel>0 with useDSMC=F for particle boundary: ',iPartBound)
      SELECT CASE (PartBound%SurfaceModel(iPartBound))
      CASE (0)
        PartBound%Reactive(iPartBound)        = .FALSE.
      CASE (5,6,7)
        PartBound%Reactive(iPartBound)        = .TRUE.
      CASE DEFAULT
        CALL abort(&
            __STAMP__&
            ,'Error in particle init: only allowed SurfaceModels: 0,5,6,7!')
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
      IF(PartBound%Reactive(iPartBound)) THEN
        CALL abort(&
          __STAMP__&
          ,'ERROR: Species swap is only supported in combination with Maxwell scattering (SurfModel = 0). PartBound: ',iPartBound)
      END IF
    END IF
    ! Dielectric Surfaces
    PartBound%Dielectric(iPartBound)      = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Dielectric')
    ! Sanity check: PartBound%Dielectric=T requires supplying species swap for every species
    IF(PartBound%Dielectric(iPartBound))THEN
      IF(PartBound%NbrOfSpeciesSwaps(iPartBound).LT.(nSpecies-BGGas%NumberOfSpecies))THEN
        CALL abort(__STAMP__,&
          'PartBound%Dielectric=T requires supplying a species swap (Part-BoundaryX-NbrOfSpeciesSwaps) for every species (except background gas species)!')
      ELSE
        DoDielectricSurfaceCharge=.TRUE.
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
    PartBound%TargetBoundCond(iPartBound)  = PartBound%RotPeriodicBC
    PartBound%RotPeriodicDir(iPartBound) = GETINT('Part-Boundary'//TRIM(hilf)//'-RotPeriodicDir','0.')
    IF(ABS(PartBound%RotPeriodicDir(iPartBound)).NE.1) THEN
      CALL abort(&
          __STAMP__&
          ,'Angle for for rotational periodicity must not be zero!')
    END IF
  CASE DEFAULT
    SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
    CALL abort(&
        __STAMP__&
        ,'Particle Boundary Condition does not exist')
  END SELECT
  PartBound%SourceBoundName(iPartBound) = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-SourceName'))

  ! Surface particle output to .h5
  PartBound%BoundaryParticleOutputHDF5(iPartBound)      = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-BoundaryParticleOutput')
  IF(PartBound%BoundaryParticleOutputHDF5(iPartBound)) DoBoundaryParticleOutputHDF5=.TRUE.
END DO
AdaptWallTemp = GETLOGICAL('Part-AdaptWallTemp')

IF(GEO%RotPeriodicBC) THEN
  GEO%RotPeriodicAxi   = GETINT('Part-RotPeriodicAxi')
  GEO%RotPeriodicAngle = GETREAL('Part-RotPeriodicAngle')
  IF(ALMOSTZERO(GEO%RotPeriodicAngle)) THEN
    CALL abort(&
        __STAMP__&
        ,'Angle for for rotational periodicity must not be zero!')
  END IF
  GEO%RotPeriodicAngle = GEO%RotPeriodicAngle / 180. * PI
END IF

! Surface particle output to .h5
IF(DoBoundaryParticleOutputHDF5)THEN
  ALLOCATE(PartStateBoundary(1:nVarPartStateBoundary,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
        __STAMP__&
        ,'ERROR in particle_init.f90: Cannot allocate PartStateBoundary array!')
  END IF
  PartStateBoundary=0.
END IF

! Set mapping from field boundary to particle boundary index
ALLOCATE(PartBound%MapToPartBC(1:nBCs))
PartBound%MapToPartBC(:)=-10
DO iPBC=1,nPartBound
  DO iBC = 1, nBCs
    IF (BoundaryType(iBC,BC_TYPE).EQ.0) THEN
      PartBound%MapToPartBC(iBC) = -1 !there are no internal BCs in the mesh, they are just in the name list!
      SWRITE(*,*)"... PartBound",iPBC,"is internal bound, no mapping needed"
    ELSEIF(BoundaryType(iBC,BC_TYPE).EQ.100)THEN
      IF(TrackingMethod.EQ.REFMAPPING)THEN
        SWRITE(UNIT_STDOUT,'(A)') ' Analyze sides are not implemented for RefMapping=T, because '//  &
                                  ' orientation of SideNormVec is unknown.'
     CALL abort(&
                __STAMP__&
                ,' Analyze-BCs cannot be used for internal reflection in general cases! ')
      END IF
    END IF
    IF (TRIM(BoundaryName(iBC)).EQ.TRIM(PartBound%SourceBoundName(iPBC))) THEN
      PartBound%MapToPartBC(iBC) = iPBC !PartBound%TargetBoundCond(iPBC)
      SWRITE(*,*)"... Mapped PartBound",iPBC,"on FieldBound",BoundaryType(iBC,1),",i.e.:",TRIM(BoundaryName(iBC))
    END IF
  END DO
END DO
! Errorhandler for PartBound-Types that could not be mapped to the
! FieldBound-Types.
DO iBC = 1,nBCs
  IF (PartBound%MapToPartBC(iBC).EQ.-10) THEN
    CALL abort(&
__STAMP__&
    ,' PartBound%MapToPartBC for Boundary is not set. iBC: :',iBC)
  END IF
END DO

!-- Floating Potential
ALLOCATE(BCdata_auxSF(1:nPartBound))
DO iPartBound=1,nPartBound
  BCdata_auxSF(iPartBound)%SideNumber=-1 !init value when not used
  BCdata_auxSF(iPartBound)%GlobalArea=0.
  BCdata_auxSF(iPartBound)%LocalArea=0.
END DO

!-- Sanity check: Deprecated voltage parameter
IF(DeprecatedVoltage) CALL abort(&
  __STAMP__&
  ,'Part-Boundary-Voltage is no longer supported. Use corresponding RefState parameter as described in the user guide.')

END SUBROUTINE InitializeVariablesPartBoundary


SUBROUTINE InitParticleBoundaryRotPeriodic()
!----------------------------------------------------------------------------------------------------------------------------------!
! Build Mapping for rotational periodicity: RotPeriodicSide -> SideID2 (Side on corresponding BC).
! In RotPeriodicBC (particle_boundary_condition.f90): SideID -> SurfSideID -> RotPeriodicSide
!                                                     RotPeriodicSide -> SideID2
! (1) counting rotational periodic sides and build mapping from SurfSideID -> RotPeriodicSide
! (2) find Side on corresponding BC and build mapping RotPeriodicSide -> SideID2 (and vice versa)
!     counting potential rotational periodic sides (for not conform meshes)
! (3) reallocate array due to number of potential rotational periodic sides
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: RotPeriodicSide2GlobalSide,nComputeNodeSurfTotalSides,SurfSide2GlobalSide,PartBound
USE MOD_Particle_Boundary_Vars  ,ONLY: RotPeriodicSideMapping, NumRotPeriodicNeigh, SurfSide2RotPeriodicSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared, NodeCoords_Shared, ElemSideNodeID_Shared, GEO, ElemInfo_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iSide, jSide, nRotPeriodicSides, SideID,SideID2, MaxNumRotPeriodicNeigh, iNode, jNode, iNeigh
INTEGER              :: NodeID, CNElemID, CNElemID2, LocSideID, k, l, m, NodeID2, LocSideID2
INTEGER,ALLOCATABLE  :: Rot2Glob_temp(:)
INTEGER,ALLOCATABLE  :: RotPeriodicSideMapping_temp(:,:)
LOGICAL              :: isMapped,SideIsMapped
REAL                 :: iNodeVec(1:3), jNodeVec(1:3)
REAL                 :: iNodeR, iNodeH, jNodeR, jNodeH, Node2Rmin, Node2Rmax, Node2Hmin, Node2Hmax
!===================================================================================================================================

ALLOCATE(Rot2Glob_temp(nComputeNodeSurfTotalSides))
ALLOCATE(SurfSide2RotPeriodicSide(nComputeNodeSurfTotalSides))
SurfSide2RotPeriodicSide(:) = -1
nRotPeriodicSides=0

! (1) counting rotational periodic sides and build mapping from SurfSideID -> RotPeriodicSide
DO iSide=1, nComputeNodeSurfTotalSides
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%RotPeriodicBC) THEN
    nRotPeriodicSides = nRotPeriodicSides + 1
    Rot2Glob_temp(nRotPeriodicSides) = SideID
    SurfSide2RotPeriodicSide(iSide) = nRotPeriodicSides
  END IF
END DO

ALLOCATE(RotPeriodicSide2GlobalSide(nRotPeriodicSides))
ALLOCATE(NumRotPeriodicNeigh(nRotPeriodicSides))
! number of potential rotational periodic sides is unknown => allocate mapping array with fixed number of 1000
! and reallocate at the end of subroutine
ALLOCATE(RotPeriodicSideMapping_temp(nRotPeriodicSides,1000))

DO iSide=1, nRotPeriodicSides
  RotPeriodicSide2GlobalSide(iSide) = Rot2Glob_temp(iSide)
  NumRotPeriodicNeigh(iSide) = 0
  RotPeriodicSideMapping_temp(iSide,1:1000) = -1
END DO

MaxNumRotPeriodicNeigh = 0
! Defining rotation matrix
SELECT CASE(GEO%RotPeriodicAxi)
  CASE(1) ! x-rotation axis
    k = 1
    l = 2
    m = 3
  CASE(2) ! x-rotation axis
    k = 2
    l = 3
    m = 1
  CASE(3) ! x-rotation axis
    k = 3
    l = 1
    m = 2
END SELECT
! (2) find Side on corresponding BC and build mapping RotPeriodicSide -> SideID2 (and vice versa)
!     counting potential rotational periodic sides (for not conform meshes)
DO iSide=1, nRotPeriodicSides
  SideID    = RotPeriodicSide2GlobalSide(iSide)
  CNElemID  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
  LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
  SideIsMapped = .FALSE.
  ! check if at least one node of iSide is inside bounding box of a Side on corresponding BC
  DO iNode=1, 4
    NodeID = ElemSideNodeID_Shared(iNode,LocSideID,CNElemID) + 1
    iNodeVec(1:3) = NodeCoords_Shared(1:3,NodeID)
    iNodeH = iNodeVec(k)
    iNodeR = SQRT(iNodeVec(l)*iNodeVec(l)+iNodeVec(m)*iNodeVec(m))
    DO jSide=1, nRotPeriodicSides
      SideID2 = RotPeriodicSide2GlobalSide(jSide)
      ! is on same RotPeriodicBC?
      IF(PartBound%RotPeriodicDir(SideInfo_Shared(SIDE_BCID,SideID)).EQ. &
        PartBound%RotPeriodicDir(SideInfo_Shared(SIDE_BCID,SideID2))) CYCLE
      isMapped = .FALSE.
      ! check wether RotPeriodicSides is already mapped
      IF(NumRotPeriodicNeigh(jSide).GT. 0) THEN
        DO iNeigh=1, NumRotPeriodicNeigh(jSide)
          IF(RotPeriodicSideMapping_temp(jSide,iNeigh).EQ.SideID) THEN
            isMapped = .TRUE.
            SideIsMapped = .TRUE.
            EXIT
          END IF
        END DO
      END IF
      IF(isMapped) CYCLE
      ! get ElemID for node mapping
      CNElemID2    = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID2))
      LocSideID2 = SideInfo_Shared(SIDE_LOCALID,SideID2)
      ! calc bounding box of jSide
      DO jNode=1, 4
        NodeID2       = ElemSideNodeID_Shared(jNode,LocSideID2,CNElemID2) + 1
        jNodeVec(1:3) = NodeCoords_Shared(1:3,NodeID2)
        jNodeH        = jNodeVec(k)
        jNodeR        = SQRT(jNodeVec(l)*jNodeVec(l)+jNodeVec(m)*jNodeVec(m))
        IF(jNode.EQ. 1) THEN
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
      IF( ( (Node2Hmin.LE.iNodeH).AND.(iNodeH.LE.Node2Hmax) ) .AND. &
          ( (Node2Rmin.LE.iNodeR).AND.(iNodeR.LE.Node2Rmax) )       ) THEN
      ! at least one node of iSide is inside bounding box of jSide =>
      !                                                     1. increase NumRotPeriodicNeigh
      !                                                     2. map:     iSide => SideID2 and
      !                                                     vise versa: jSide => SideID  s.o.
        NumRotPeriodicNeigh(iSide) = NumRotPeriodicNeigh(iSide) + 1
        IF(NumRotPeriodicNeigh(iSide).GT. 1000) THEN
          CALL abort(&
              __STAMP__&
              ,' ERROR: Number of rotational periodic side exceed fixed number of 1000!.')
        END IF
        RotPeriodicSideMapping_temp(iSide,NumRotPeriodicNeigh(iSide)) = SideID2
        NumRotPeriodicNeigh(jSide) = NumRotPeriodicNeigh(jSide) + 1
        IF(NumRotPeriodicNeigh(jSide).GT. 1000) THEN
          CALL abort(&
              __STAMP__&
              ,' ERROR: Number of rotational periodic side exceed fixed number of 1000!.')
        END IF
        RotPeriodicSideMapping_temp(jSide,NumRotPeriodicNeigh(jSide)) = SideID
        SideIsMapped = .TRUE.
      END IF
    END DO
  END DO
  IF(.NOT.SideIsMapped) THEN
    IF(ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,SideID)).NE.3) THEN
      CALL abort(&
          __STAMP__&
          ,' ERROR: One rot periodic side did not find a corresponding side.')
    ELSE
      NumRotPeriodicNeigh(iSide) = 1
      RotPeriodicSideMapping_temp(iSide,NumRotPeriodicNeigh(iSide)) = -1
    END IF
  END IF
END DO
! (3) reallocate array due to number of potential rotational periodic sides
MaxNumRotPeriodicNeigh = MAXVAL(NumRotPeriodicNeigh)
ALLOCATE(RotPeriodicSideMapping(nRotPeriodicSides,MaxNumRotPeriodicNeigh))
DO iSide=1, nRotPeriodicSides
  DO iNeigh=1, MaxNumRotPeriodicNeigh
    RotPeriodicSideMapping(iSide,iNeigh) = RotPeriodicSideMapping_temp(iSide,iNeigh)
  END DO
END DO

END SUBROUTINE InitParticleBoundaryRotPeriodic


SUBROUTINE InitAdaptiveWallTemp()
!===================================================================================================================================
!> Allocate shared array for the side-local wall temperature: BoundaryWallTemp_Shared and initialize with read-in wall temperature
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, nComputeNodeSurfTotalSides, BoundaryWallTemp, SurfSide2GlobalSide,nSurfSample
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared
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
INTEGER                           :: firstSide, lastSide, iSide, SideID, iBC
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND)    :: MPISharedSize
#endif
!===================================================================================================================================
IF (.NOT.(ANY(PartBound%UseAdaptedWallTemp))) RETURN

#if USE_MPI
!> Then shared arrays for boundary sampling
MPISharedSize = INT((nSurfSample*nSurfSample*nComputeNodeSurfTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),BoundaryWallTemp_Shared_Win,BoundaryWallTemp_Shared)
CALL MPI_WIN_LOCK_ALL(0,BoundaryWallTemp_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  BoundaryWallTemp_Shared = 0.
END IF
BoundaryWallTemp => BoundaryWallTemp_Shared
CALL MPI_WIN_SYNC(BoundaryWallTemp_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
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
  IF (PartBound%MomentumACC(iBC).GT.0.0) BoundaryWallTemp(:,:,iSide) = PartBound%WallTemp(iBC)
END DO 

#if USE_MPI
CALL MPI_WIN_SYNC(BoundaryWallTemp_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

END SUBROUTINE InitAdaptiveWallTemp

SUBROUTINE InitializeVariablesAuxBC()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Boundary_Vars ,ONLY: PartAuxBC
USE MOD_Particle_Boundary_Vars ,ONLY: nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol,UseAuxBCs
USE MOD_Particle_Boundary_Tools,ONLY: MarkAuxBCElems
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPartBound, iSwaps, MaxNbrOfSpeciesSwaps
INTEGER               :: iAuxBC, nAuxBCplanes, nAuxBCcylinders, nAuxBCcones, nAuxBCparabols
CHARACTER(32)         :: hilf , hilf2
CHARACTER(200)        :: tmpString
REAL                  :: n_vec(3), cos2, rmax
REAL, DIMENSION(3,1)  :: norm,norm1,norm2
REAL, DIMENSION(3,3)  :: rot1, rot2
REAL                  :: alpha1, alpha2
!===================================================================================================================================
nAuxBCs=GETINT('Part-nAuxBCs','0')
IF (nAuxBCs.GT.0) THEN
  UseAuxBCs=.TRUE.
  ALLOCATE (AuxBCType(1:nAuxBCs) &
            ,AuxBCMap(1:nAuxBCs) )
  AuxBCMap=0
  !- Read in BC parameters
  ALLOCATE(PartAuxBC%TargetBoundCond(1:nAuxBCs))
  ALLOCATE(PartAuxBC%MomentumACC(1:nAuxBCs))
  ALLOCATE(PartAuxBC%WallTemp(1:nAuxBCs))
  ALLOCATE(PartAuxBC%TransACC(1:nAuxBCs))
  ALLOCATE(PartAuxBC%VibACC(1:nAuxBCs))
  ALLOCATE(PartAuxBC%RotACC(1:nAuxBCs))
  ALLOCATE(PartAuxBC%ElecACC(1:nAuxBCs))
  ALLOCATE(PartAuxBC%Resample(1:nAuxBCs))
  ALLOCATE(PartAuxBC%WallVelo(1:3,1:nAuxBCs))
  ALLOCATE(PartAuxBC%NbrOfSpeciesSwaps(1:nAuxBCs))
  !--determine MaxNbrOfSpeciesSwaps for correct allocation
  MaxNbrOfSpeciesSwaps=0
  DO iPartBound=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iPartBound
    PartAuxBC%NbrOfSpeciesSwaps(iPartBound)= GETINT('Part-AuxBC'//TRIM(hilf)//'-NbrOfSpeciesSwaps','0')
    MaxNbrOfSpeciesSwaps=max(PartAuxBC%NbrOfSpeciesSwaps(iPartBound),MaxNbrOfSpeciesSwaps)
  END DO
  IF (MaxNbrOfSpeciesSwaps.gt.0) THEN
    ALLOCATE(PartAuxBC%ProbOfSpeciesSwaps(1:nAuxBCs))
    ALLOCATE(PartAuxBC%SpeciesSwaps(1:2,1:MaxNbrOfSpeciesSwaps,1:nAuxBCs))
  END IF
  !--
  DO iPartBound=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iPartBound
    tmpString = TRIM(GETSTR('Part-AuxBC'//TRIM(hilf)//'-Condition','open'))
    SELECT CASE (TRIM(tmpString))
    CASE('open')
      PartAuxBC%TargetBoundCond(iPartBound) = PartAuxBC%OpenBC          ! definitions see typesdef_pic
    CASE('reflective')
      PartAuxBC%TargetBoundCond(iPartBound) = PartAuxBC%ReflectiveBC
      PartAuxBC%MomentumACC(iPartBound)     = GETREAL('Part-AuxBC'//TRIM(hilf)//'-MomentumACC')
      PartAuxBC%WallTemp(iPartBound)        = GETREAL('Part-AuxBC'//TRIM(hilf)//'-WallTemp')
      PartAuxBC%TransACC(iPartBound)        = GETREAL('Part-AuxBC'//TRIM(hilf)//'-TransACC')
      PartAuxBC%VibACC(iPartBound)          = GETREAL('Part-AuxBC'//TRIM(hilf)//'-VibACC')
      PartAuxBC%RotACC(iPartBound)          = GETREAL('Part-AuxBC'//TRIM(hilf)//'-RotACC')
      PartAuxBC%ElecACC(iPartBound)         = GETREAL('Part-AuxBC'//TRIM(hilf)//'-ElecACC')
      PartAuxBC%Resample(iPartBound)        = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-Resample')
      PartAuxBC%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-WallVelo',3)
      IF (PartAuxBC%NbrOfSpeciesSwaps(iPartBound).gt.0) THEN
        !read Species to be changed at wall (in, out), out=0: delete
        PartAuxBC%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-AuxBC'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
        DO iSwaps=1,PartAuxBC%NbrOfSpeciesSwaps(iPartBound)
          WRITE(UNIT=hilf2,FMT='(I0)') iSwaps
          PartAuxBC%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
            GETINTARRAY('Part-AuxBC'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
        END DO
      END IF
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC Condition does not exists: ', TRIM(tmpString)
      CALL abort(&
        __STAMP__&
        ,'AuxBC Condition does not exist')
    END SELECT
  END DO
  !- read and count types
  nAuxBCplanes = 0
  nAuxBCcylinders = 0
  nAuxBCcones = 0
  nAuxBCparabols = 0
  DO iAuxBC=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iAuxBC
    AuxBCType(iAuxBC) = TRIM(GETSTR('Part-AuxBC'//TRIM(hilf)//'-Type','plane'))
    SELECT CASE (TRIM(AuxBCType(iAuxBC)))
    CASE ('plane')
      nAuxBCplanes = nAuxBCplanes + 1
      AuxBCMap(iAuxBC) = nAuxBCplanes
    CASE ('cylinder')
      nAuxBCcylinders = nAuxBCcylinders + 1
      AuxBCMap(iAuxBC) = nAuxBCcylinders
    CASE ('cone')
      nAuxBCcones = nAuxBCcones + 1
      AuxBCMap(iAuxBC) = nAuxBCcones
    CASE ('parabol')
      nAuxBCparabols = nAuxBCparabols + 1
      AuxBCMap(iAuxBC) = nAuxBCparabols
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist')
    END SELECT
  END DO
  !- allocate type-specifics
  IF (nAuxBCplanes.GT.0) THEN
    ALLOCATE (AuxBC_plane(1:nAuxBCplanes))
  END IF
  IF (nAuxBCcylinders.GT.0) THEN
    ALLOCATE (AuxBC_cylinder(1:nAuxBCcylinders))
  END IF
  IF (nAuxBCcones.GT.0) THEN
    ALLOCATE (AuxBC_cone(1:nAuxBCcones))
  END IF
  IF (nAuxBCparabols.GT.0) THEN
    ALLOCATE (AuxBC_parabol(1:nAuxBCparabols))
  END IF
  !- read type-specifics
  DO iAuxBC=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iAuxBC
    SELECT CASE (TRIM(AuxBCType(iAuxBC)))
    CASE ('plane')
      AuxBC_plane(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_plane(AuxBCMap(iAuxBC))%radius)
      AuxBC_plane(AuxBCMap(iAuxBC))%radius= GETREAL('Part-AuxBC'//TRIM(hilf)//'-radius',TRIM(hilf2))
      n_vec                               = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-n_vec',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-n_vec is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_plane(AuxBCMap(iAuxBC))%n_vec = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
    CASE ('cylinder')
      AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                                  = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-axis',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_cylinder(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
      AuxBC_cylinder(AuxBCMap(iAuxBC))%radius  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-radius','1.')
      WRITE(UNIT=hilf2,FMT='(G0)') -HUGE(AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmin',TRIM(hilf2))
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmax',TRIM(hilf2))
      AuxBC_cylinder(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-inwards','.TRUE.')
    CASE ('cone')
      AuxBC_cone(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                              = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-axis',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_cone(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
      AuxBC_cone(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmin','0.')
      IF (AuxBC_cone(AuxBCMap(iAuxBC))%lmin.LT.0.) CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-lminis .lt. zero for AuxBC',iAuxBC)
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_cone(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cone(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmax',TRIM(hilf2))
      rmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-rmax','0.')
      ! either define rmax at lmax or the halfangle
      IF (rmax.EQ.0.) THEN
        AuxBC_cone(AuxBCMap(iAuxBC))%halfangle  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-halfangle','45.')*PI/180.
      ELSE
        AuxBC_cone(AuxBCMap(iAuxBC))%halfangle  = ATAN(rmax/AuxBC_cone(AuxBCMap(iAuxBC))%lmax)
      END IF
      IF (AuxBC_cone(AuxBCMap(iAuxBC))%halfangle.LE.0.) CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-halfangle is .le. zero for AuxBC',iAuxBC)
      AuxBC_cone(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-inwards','.TRUE.')
      cos2 = COS(AuxBC_cone(AuxBCMap(iAuxBC))%halfangle)**2
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,1) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(1)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/cos2,0.,0./)
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,2) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(2)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/0.,cos2,0./)
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,3) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(3)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/0.,0.,cos2/)
    CASE ('parabol')
      AuxBC_parabol(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                              = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-axis',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_parabol(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
      AuxBC_parabol(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmin','0.')
      IF (AuxBC_parabol(AuxBCMap(iAuxBC))%lmin.LT.0.) CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-lmin is .lt. zero for AuxBC',iAuxBC)
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_parabol(AuxBCMap(iAuxBC))%lmin)
      AuxBC_parabol(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmax',TRIM(hilf2))
      AuxBC_parabol(AuxBCMap(iAuxBC))%zfac  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-zfac','1.')
      AuxBC_parabol(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-inwards','.TRUE.')

      norm(:,1)=AuxBC_parabol(AuxBCMap(iAuxBC))%axis
      IF (.NOT.ALMOSTZERO(SQRT(norm(1,1)**2+norm(3,1)**2))) THEN !collinear with y?
        alpha1=ATAN2(norm(1,1),norm(3,1))
        CALL roty(rot1,alpha1)
        norm1=MATMUL(rot1,norm)
      ELSE
        alpha1=0.
        CALL ident(rot1)
        norm1=norm
      END IF
      IF (.NOT.ALMOSTZERO(SQRT(norm1(2,1)**2+norm1(3,1)**2))) THEN !collinear with x?
        alpha2=-ATAN2(norm1(2,1),norm1(3,1))
        CALL rotx(rot2,alpha2)
        norm2=MATMUL(rot2,norm1)
      ELSE
        CALL abort(&
          __STAMP__&
          ,'vector is collinear with x-axis. this should not be possible... AuxBC:',iAuxBC)
      END IF
      AuxBC_parabol(AuxBCMap(iAuxBC))%rotmatrix(:,:)=MATMUL(rot2,rot1)
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(:,:)=0.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(1,1)=1.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(2,2)=1.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(3,3)=0.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(3,4)=-0.5*AuxBC_parabol(AuxBCMap(iAuxBC))%zfac
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(4,3)=-0.5*AuxBC_parabol(AuxBCMap(iAuxBC))%zfac
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist for AuxBC',iAuxBC)
    END SELECT
  END DO
  CALL MarkAuxBCElems()
ELSE
  UseAuxBCs=.FALSE.
END IF

END SUBROUTINE InitializeVariablesAuxBC

SUBROUTINE rotx(mat,a)
IMPLICIT NONE
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
REAL, INTENT(IN) :: a
mat(:,1)=(/1.0 , 0.     , 0.  /)
mat(:,2)=(/0.0 , COS(a) ,-SIN(a)/)
mat(:,3)=(/0.0 , SIN(a) , COS(a)/)
END SUBROUTINE


SUBROUTINE roty(mat,a)
IMPLICIT NONE
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
REAL, INTENT(IN) :: a
mat(:,1)=(/ COS(a) , 0., SIN(a)/)
mat(:,2)=(/ 0.     , 1., 0.  /)
mat(:,3)=(/-SIN(a) , 0., COS(a)/)
END SUBROUTINE


SUBROUTINE ident(mat)
IMPLICIT NONE
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
INTEGER :: j
mat = 0.
FORALL(j = 1:3) mat(j,j) = 1.
END SUBROUTINE


SUBROUTINE FinalizeParticleBoundary()
!----------------------------------------------------------------------------------------------------------------------------------!
!> Finalize particle boundary variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars
#if USE_MPI
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
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
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%WallTemp2)
SDEALLOCATE(PartBound%WallTempDelta)
SDEALLOCATE(PartBound%TempGradStart)
SDEALLOCATE(PartBound%TempGradEnd)
SDEALLOCATE(PartBound%TempGradVec)
SDEALLOCATE(PartBound%TransACC)
SDEALLOCATE(PartBound%VibACC)
SDEALLOCATE(PartBound%RotACC)
SDEALLOCATE(PartBound%ElecACC)
SDEALLOCATE(PartBound%Resample)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%RotVelo)
SDEALLOCATE(PartBound%RotFreq)
SDEALLOCATE(PartBound%RotOrg)
SDEALLOCATE(PartBound%RotAxi)
SDEALLOCATE(PartBound%RotPeriodicDir)
SDEALLOCATE(PartBound%Voltage)
SDEALLOCATE(PartBound%NbrOfSpeciesSwaps)
SDEALLOCATE(PartBound%ProbOfSpeciesSwaps)
SDEALLOCATE(PartBound%SpeciesSwaps)
SDEALLOCATE(PartBound%MapToPartBC)
SDEALLOCATE(PartBound%SurfaceModel)
SDEALLOCATE(PartBound%Reactive)
SDEALLOCATE(PartBound%Dielectric)
SDEALLOCATE(PartBound%BoundaryParticleOutputHDF5)
SDEALLOCATE(PartBound%RadiativeEmissivity)
SDEALLOCATE(PartStateBoundary)
! Rotational periodic boundary condition
SDEALLOCATE(RotPeriodicSide2GlobalSide)
SDEALLOCATE(NumRotPeriodicNeigh)
SDEALLOCATE(RotPeriodicSideMapping)
SDEALLOCATE(SurfSide2RotPeriodicSide)
! Adaptive wall temperature (e.g. calculate from sampled heat flux)
IF (ANY(PartBound%UseAdaptedWallTemp)) THEN
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
END SUBROUTINE FinalizeParticleBoundary

END MODULE MOD_Particle_Boundary_Init
