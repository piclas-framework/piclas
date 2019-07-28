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

MODULE MOD_Particle_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitParticleAnalyze
  MODULE PROCEDURE InitParticleAnalyze
END INTERFACE

INTERFACE FinalizeParticleAnalyze
  MODULE PROCEDURE FinalizeParticleAnalyze
END INTERFACE

INTERFACE AnalyzeParticles
  MODULE PROCEDURE AnalyzeParticles
END INTERFACE

INTERFACE CalcShapeEfficiencyR
  MODULE PROCEDURE CalcShapeEfficiencyR
END INTERFACE

INTERFACE CalcEkinPart
  MODULE PROCEDURE CalcEkinPart
END INTERFACE

INTERFACE CalcPowerDensity
  MODULE PROCEDURE CalcPowerDensity
END INTERFACE

INTERFACE PartIsElectron
  MODULE PROCEDURE PartIsElectron
END INTERFACE

INTERFACE CalculatePartElemData
  MODULE PROCEDURE CalculatePartElemData
END INTERFACE

INTERFACE WriteParticleTrackingData
  MODULE PROCEDURE WriteParticleTrackingData
END INTERFACE

#ifdef CODE_ANALYZE
INTERFACE AnalyticParticleMovement
  MODULE PROCEDURE AnalyticParticleMovement
END INTERFACE
#endif /*CODE_ANALYZE*/

PUBLIC:: InitParticleAnalyze, FinalizeParticleAnalyze!, CalcPotentialEnergy
PUBLIC:: CalcEkinPart, AnalyzeParticles, PartIsElectron
PUBLIC:: CalcPowerDensity
PUBLIC:: CalculatePartElemData
PUBLIC:: WriteParticleTrackingData
#ifdef CODE_ANALYZE
PUBLIC:: AnalyticParticleMovement
#if (PP_TimeDiscMethod==42)
PUBLIC :: ElectronicTransition, WriteEletronicTransition
#endif
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
PUBLIC::DefineParametersParticleAnalyze

CONTAINS

!==================================================================================================================================
!> Define parameters for analyze if particles (.csv output)
!==================================================================================================================================
SUBROUTINE DefineParametersParticleAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
!USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Analyze")

CALL prms%CreateIntOption(      'Part-AnalyzeStep'        , 'Analyze is performed each Nth time step','1')
CALL prms%CreateLogicalOption(  'CalcTotalEnergy'         , 'Calculate Total Energy. Output file is Database.csv','.FALSE.')
CALL prms%CreateLogicalOption(  'PIC-VerifyCharge'        , 'Validate the charge after each deposition'//&
                                                            'and write an output in std.out','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcIonizationDegree'    , 'Compute the ionization degree in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPointsPerShapeFunction','Compute the points per shape function in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPlasmaParameter'     ,'Compute the plasma parameter N_D in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPointsPerDebyeLength', 'Compute the points per Debye length in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcDebyeLength'         , 'Compute the Debye length in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPICTimeStep'         , 'Compute the HDG time step in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcElectronTemperature' , 'Compute the electron temperature in each cell','.FALSE.')
!CALL prms%CreateLogicalOption(  'ElectronTemperatureIsMaxwell', 'Flag if  electron temperature is assumed to be Maxwellian in each cell','.TRUE.')
CALL prms%CreateLogicalOption(  'CalcElectronIonDensity'     , 'Compute the electron density in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPlasmaFrequency'     , 'Compute the electron frequency in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcCharge'              , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Compute the whole deposited charge,'//&
                                                            ' absolute and relative charge error','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcKineticEnergy'       , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate Kinetic Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcInternalEnergy'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate Internal Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcTemp'                , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate Translational temperature.'&
                                                          ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPartBalance'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the Particle Power Balance'//&
                                                            '- input and outflow energy of all particles','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcVelos'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate thermal and flow velocities.'//&
                                                            'if CalcVelos = T VelocityDirections = (/[int],[int],[int],[int]/)  '//&
                                                            'Switching dimensions for CalcVelos on (1) or off (0)\n'//&
                                                            '(/v_x,v_y,v_z,|v|/) ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcLaserInteraction'     , 'Compute laser-plasma interaction properties such as maximum '//&
                                                             'particle energy per species.','.FALSE.')
CALL prms%CreateRealOption(     'LaserInteractionEkinMaxRadius','maximum radius (x- and y-dir) of particle to be considered for '//&
                                                                'Ekin maximum calculation (default is HUGE) '//&
                                                                'OR if LaserInteractionEkinMaxZPosMin condition is true')
CALL prms%CreateRealOption(     'LaserInteractionEkinMaxZPosMin','minimum z-position of particle to be considered for Ekin '//&
                                                                 'maximum calculation (default is -1.*HUGE) '//&
                                                                 'OR if LaserInteractionEkinMaxRadius condition is true')
CALL prms%CreateIntArrayOption( 'VelocityDirections'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'x,y,z,abs -> 0/1 = T/F. (please note: CalcVelos)'&
                                                          ,'1 , 1 , 1 , 1')
CALL prms%CreateLogicalOption(  'Part-TrackPosition'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Track particle position','.FALSE.')
CALL prms%CreateLogicalOption(  'printDiff'               , 'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateRealOption(     'printDiffTime'           , 'TODO-DEFINE-PARAMETER','12.')
CALL prms%CreateRealArrayOption('printDiffVec'            , 'TODO-DEFINE-PARAMETER','0. , 0. , 0. , 0. , 0. , 0.')
CALL prms%CreateLogicalOption(  'CalcNumSpec'             , 'Calculate the number of simulation particles per species for the '//&
                                                            'complete domain','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcNumDens'             , 'Calculate the number density per species for the complete domain' &
                                                          ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcCollRates'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the collision rates per '//&
                                                            'collision pair','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcReacRates'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the reaction rate per reaction'&
                                                          ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcShapeEfficiency'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Use efficiency methods for shape functions.'&
                                                          , '.FALSE.')
CALL prms%CreateStringOption(   'CalcShapeEfficiencyMethod', 'TODO-DEFINE-PARAMETER\n'//&
                                                             'Choose between "AllParts" and '//&
                                                             '"SomeParts", to either use all particles or a certain percentage'//&
                                                             ' (ShapeEfficiencyNumber) of the currently used particles','AllParts')
CALL prms%CreateIntOption(      'ShapeEfficiencyNumber'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                             'Percentage of currently used particles is used.'&
                                                           , '100')
CALL prms%CreateLogicalOption(  'IsRestart'                , 'TODO-DEFINE-PARAMETER\n'//&
                                                             'Flag, if the current calculation is a restart. '&
                                                           , '.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPorousBCInfo'         , 'Calculate output of porous BCs such pumping speed, removal '//&
                                                             'probability and pressure (normalized with the given pressure). '//&
                                                             'Values are averaged over the whole porous BC.' , '.FALSE.')

CALL prms%CreateLogicalOption(  'CalcCoupledPower'         , ' Calculate output of Power that is coupled into plasma' , '.FALSE.')

END SUBROUTINE DefineParametersParticleAnalyze

SUBROUTINE InitParticleAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: PI
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars
USE MOD_ReadInTools           ,ONLY: GETLOGICAL, GETINT, GETSTR, GETINTARRAY, GETREALARRAY, GETREAL
USE MOD_Particle_Vars         ,ONLY: nSpecies, VarTimeStep, PDM
USE MOD_PICDepo_Vars          ,ONLY: DoDeposition
USE MOD_IO_HDF5               ,ONLY: AddToElemData,ElementOut
USE MOD_PICDepo_Vars          ,ONLY: r_sf
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_ReadInTools           ,ONLY: PrintOption
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#if (PP_TimeDiscMethod == 42)
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
USE MOD_Particle_Vars         ,ONLY: ManualTimeStep
USE MOD_Restart_Vars          ,ONLY: RestartTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: dir, VeloDirs_hilf(4),iElem
REAL          :: DOF,VolumeShapeFunction
CHARACTER(32) :: hilf
!===================================================================================================================================
IF (ParticleAnalyzeInitIsDone) THEN
CALL abort(__STAMP__,&
'InitParticleAnalyse already called.',999,999.)
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

! Average number of points per shape function: max. number allowed is (PP_N+1)^3
CalcPointsPerShapeFunction = GETLOGICAL('CalcPointsPerShapeFunction','.FALSE.')
IF(CalcPointsPerShapeFunction)THEN
  ! calculate cell local number excluding neighbor DOFs
  ALLOCATE( PPSCell(1:PP_nElems) )
  PPSCell=0.0
  CALL AddToElemData(ElementOut,'PPSCell',RealArray=PPSCell(1:PP_nElems))
  ! assume Cartesian grid and calculate to total number including neighbor DOFs
  ALLOCATE( PPSCellEqui(1:PP_nElems) )
  PPSCellEqui=0.0
  CALL AddToElemData(ElementOut,'PPSCellEqui',RealArray=PPSCellEqui(1:PP_nElems))
  !
  VolumeShapeFunction = 4./3.*PI*(r_sf**3)
  CALL PrintOption('VolumeShapeFunction','OUTPUT',RealOpt=VolumeShapeFunction)
  DOF                 = REAL((PP_N+1)**3)
  CALL PrintOption('Max DOFs in Shape-Function per cell','OUTPUT',RealOpt=DOF)
  DO iElem = 1, nElems
    PPSCell(iElem)     = MIN(1.,VolumeShapeFunction/GEO%Volume(iElem)) * DOF
    PPSCellEqui(iElem) =       (VolumeShapeFunction/GEO%Volume(iElem)) * DOF
  END DO ! iElem = 1, nElems
END IF

!--------------------------------------------------------------------------------------------------------------------
! get derived particle properties
! (Note that for IMD/TTM initialization these values are calculated from the TTM grid values)
!--------------------------------------------------------------------------------------------------------------------
! PointsPerDebyeLength: PPD = (p+1)*lambda_D/L_cell
! p:        Polynomial degree
! lambda_D: Debye length
! L_cell:   Characteristic ceill length -> V_cell^(1/3)
CalcPointsPerDebyeLength       = GETLOGICAL('CalcPointsPerDebyeLength','.FALSE.')
IF(CalcPointsPerDebyeLength)THEN
  ALLOCATE( PPDCell(1:PP_nElems) )
  PPDCell=0.0
  CALL AddToElemData(ElementOut,'PPDCell',RealArray=PPDCell(1:PP_nElems))
END IF

! Plasma parameter
CalcPlasmaParameter   = GETLOGICAL('CalcPlasmaParameter','.FALSE.')
IF(CalcPlasmaParameter)THEN
  ALLOCATE( PlasmaParameterCell(1:PP_nElems) )
  PlasmaParameterCell=0.0
  CALL AddToElemData(ElementOut,'PlasmaParameterCell',RealArray=PlasmaParameterCell(1:PP_nElems))
END IF

! Debye Length
CalcDebyeLength       = GETLOGICAL('CalcDebyeLength','.FALSE.')
IF(CalcPointsPerDebyeLength.OR.CalcPlasmaParameter) CalcDebyeLength=.TRUE.
IF(CalcDebyeLength)THEN
  ALLOCATE( DebyeLengthCell(1:PP_nElems) )
  DebyeLengthCell=0.0
  CALL AddToElemData(ElementOut,'DebyeLengthCell',RealArray=DebyeLengthCell(1:PP_nElems))
END IF

! Ionization degree and quasi-neutrality
CalcIonizationDegree = GETLOGICAL('CalcIonizationDegree','.FALSE.')
IF(CalcDebyeLength) CalcIonizationDegree=.TRUE.
IF(CalcIonizationDegree)THEN
  ! degree of ionization
  ALLOCATE( IonizationCell(1:PP_nElems) )
  IonizationCell=0.0
  CALL AddToElemData(ElementOut,'IonizationCell',RealArray=IonizationCell(1:PP_nElems))
  ! quasi neutrality
  ALLOCATE( QuasiNeutralityCell(1:PP_nElems) )
  QuasiNeutralityCell=0.0
  CALL AddToElemData(ElementOut,'QuasiNeutralityCell',RealArray=QuasiNeutralityCell(1:PP_nElems))
END IF

! PIC Time Step Approximation
CalcPICTimeStep       = GETLOGICAL('CalcPICTimeStep','.FALSE.')
IF(CalcPICTimeStep)THEN
  ALLOCATE( PICTimeStepCell(1:PP_nElems) )
  PICTimeStepCell=0.0
  CALL AddToElemData(ElementOut,'PICTimeStepCell',RealArray=PICTimeStepCell(1:PP_nElems))
END IF

! Plasma Frequency
CalcPlasmaFrequency   = GETLOGICAL('CalcPlasmaFrequency','.FALSE.')
IF(CalcPICTimeStep) CalcPlasmaFrequency=.TRUE.
IF(CalcPlasmaFrequency)THEN
  ALLOCATE( PlasmaFrequencyCell(1:PP_nElems) )
  PlasmaFrequencyCell=0.0
  CALL AddToElemData(ElementOut,'PlasmaFrequencyCell',RealArray=PlasmaFrequencyCell(1:PP_nElems))
END IF

! Electron Density
CalcElectronIonDensity   = GETLOGICAL('CalcElectronIonDensity','.FALSE.')
IF(CalcDebyeLength.OR.CalcPlasmaFrequency.OR.CalcIonizationDegree) CalcElectronIonDensity=.TRUE.
IF(CalcElectronIonDensity) THEN
  ! electrons
  ALLOCATE( ElectronDensityCell(1:PP_nElems) )
  ElectronDensityCell=0.0
  CALL AddToElemData(ElementOut,'ElectronDensityCell',RealArray=ElectronDensityCell(1:PP_nElems))
  ! ions
  ALLOCATE( IonDensityCell(1:PP_nElems) )
  IonDensityCell=0.0
  CALL AddToElemData(ElementOut,'IonDensityCell',RealArray=IonDensityCell(1:PP_nElems))
  ! neutrals
  ALLOCATE( NeutralDensityCell(1:PP_nElems) )
  NeutralDensityCell=0.0
  CALL AddToElemData(ElementOut,'NeutralDensityCell',RealArray=NeutralDensityCell(1:PP_nElems))
  ! charge number
  ALLOCATE( ChargeNumberCell(1:PP_nElems) )
  ChargeNumberCell=0.0
  CALL AddToElemData(ElementOut,'ChargeNumberCell',RealArray=ChargeNumberCell(1:PP_nElems))
END IF

! Electron Temperature
CalcElectronTemperature   = GETLOGICAL('CalcElectronTemperature','.FALSE.')
IF(CalcDebyeLength.OR.CalcPlasmaFrequency) CalcElectronTemperature=.TRUE.
IF(CalcElectronTemperature)THEN
  !ElectronTemperatureIsMaxwell=GETLOGICAL('ElectronTemperatureIsMaxwell','.TRUE.')
  ALLOCATE( ElectronTemperatureCell(1:PP_nElems) )
  ElectronTemperatureCell=0.0
  CALL AddToElemData(ElementOut,'ElectronTemperatureCell',RealArray=ElectronTemperatureCell(1:PP_nElems))
END IF
!--------------------------------------------------------------------------------------------------------------------
! PartAnalyzeStep: The interval for the particle analyze output routines (write-out into PartAnalyze.csv)
!             = 1: Analyze and output every time step
!             = 0: Single output at the end, averaged over number of iterations (HUGE: MOD function can still be used to determine
!                  whether an output has to be performed)
!             = N: Analyze and output every Nth time step, average over N number of iterations
PartAnalyzeStep = GETINT('Part-AnalyzeStep','1')
IF (PartAnalyzeStep.EQ.0) PartAnalyzeStep = HUGE(PartAnalyzeStep)

#if (PP_TimeDiscMethod == 42)
  IF(PartAnalyzeStep.NE.HUGE(PartAnalyzeStep)) THEN
    IF(MOD(NINT((TEnd-RestartTime)/ManualTimeStep,8),PartAnalyzeStep).NE.0) THEN
      SWRITE(UNIT_stdOut,'(A,I0)') 'NINT((TEnd-RestartTime)/ManualTimeStep) = ',NINT((TEnd-RestartTime)/ManualTimeStep,8)
      SWRITE(UNIT_stdOut,'(A,I0)') '                        PartAnalyzeStep = ',PartAnalyzeStep
      CALL abort(&
        __STAMP__&
        ,'Please specify a PartAnalyzeStep, which is a factor of the total number of iterations!')
    END IF
  END IF
#endif

DoPartAnalyze = .FALSE.
! only verifycharge and CalcCharge if particles are deposited onto the grid
DoVerifyCharge= .FALSE.
CalcCharge = .FALSE.
IF(DoDeposition) THEN
  DoVerifyCharge = GETLOGICAL('PIC-VerifyCharge','.FALSE.')
  CalcCharge = GETLOGICAL('CalcCharge','.FALSE.')
  IF(CalcCharge) DoPartAnalyze = .TRUE.
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' Deposition is switched of. VerifyCharge and CalcCharge are deactivated!'
END IF

CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
! Laser-plasma interaction analysis
CalcLaserInteraction = GETLOGICAL('CalcLaserInteraction')
IF(CalcLaserInteraction)THEN
  ! set boundaries in order to exclude particles near the boundary (nonphysical velocities)
  WRITE(UNIT=hilf,FMT=WRITEFORMAT) 1.0E200!HUGE(1.0) -> HUGE produces IEEE overflow
  LaserInteractionEkinMaxRadius  = GETREAL('LaserInteractionEkinMaxRadius',TRIM(hilf))
  WRITE(UNIT=hilf,FMT=WRITEFORMAT) -1.0E200!-1.*HUGE(1.0) -> HUGE produces IEEE overflow
  LaserInteractionEkinMaxZPosMin = GETREAL('LaserInteractionEkinMaxZPosMin',TRIM(hilf))
  CalcEkin=.TRUE.
END IF

CalcEint = GETLOGICAL('CalcInternalEnergy','.FALSE.')
CalcTemp = GETLOGICAL('CalcTemp','.FALSE.')
IF(CalcTemp.OR.CalcEint) DoPartAnalyze = .TRUE.
IF(CalcEkin) DoPartAnalyze = .TRUE.
IF(nSpecies.GT.1) THEN
  nSpecAnalyze = nSpecies + 1
ELSE
  nSpecAnalyze = 1
END IF

CalcCoupledPower = GETLOGICAL('CalcCoupledPower','.FALSE.')

IF(CalcCoupledPower) THEN
  DoPartAnalyze = .TRUE.
#if !((PP_TimeDiscMethod==500) || (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506) || (PP_TimeDiscMethod==509))
  CALL abort(__STAMP__,&
            'ERROR: CalcCoupledPower is not implemented yet with the chosen time discretization method!')
#endif
END IF

! compute number of entering and leaving particles and their energy
CalcPartBalance = GETLOGICAL('CalcPartBalance','.FALSE.')
IF (CalcPartBalance) THEN
  DoPartAnalyze = .TRUE.
  SDEALLOCATE(nPartIn)
  SDEALLOCATE(nPartOut)
  SDEALLOCATE(PartEkinIn)
  SDEALLOCATE(PartEkinOut)
  ALLOCATE( nPartIn(1:nSpecAnalyze)     &
          , nPartOut(1:nSpecAnalyze)    &
          , PartEkinOut(1:nSpecAnalyze) &
          , PartEkinIn(1:nSpecAnalyze)  )
  nPartIn=0
  nPartOut=0
  PartEkinOut=0.
  PartEkinIn=0.
#if defined(LSERK) || defined(ROS) || defined(IMPA)
  SDEALLOCATE( nPartInTmp)
  SDEALLOCATE( PartEkinInTmp)
  ALLOCATE( nPartInTmp(1:nSpecAnalyze)     &
          , PartEkinInTmp(1:nSpecAnalyze)  )
  PartEkinInTmp=0.
  nPartInTmp=0
#endif
END IF
TrackParticlePosition = GETLOGICAL('Part-TrackPosition','.FALSE.')
IF(TrackParticlePosition)THEN
  IF(nProcessors.GT.1)THEN
    CALL abort(&
        __STAMP__&
        ,'Part-TrackPosition=T is currently not supported in combination with more than 1 proc!')
  ELSE
    IF(PDM%ParticleVecLength.GT.1)THEN
    CALL abort(&
        __STAMP__&
        ,'Part-TrackPosition=T is currently not supported in combination with more than 1 particle!')
    END IF
  END IF
  printDiff=GETLOGICAL('printDiff','.FALSE.')
  IF(printDiff)THEN
    printDiffTime=GETREAL('printDiffTime','12.')
    printDiffVec=GETREALARRAY('printDiffVec',6,'0.,0.,0.,0.,0.,0.')
  END IF
END IF
CalcNumSpec   = GETLOGICAL('CalcNumSpec','.FALSE.')
CalcNumDens   = GETLOGICAL('CalcNumDens','.FALSE.')
CalcCollRates = GETLOGICAL('CalcCollRates','.FALSE.')
CalcReacRates = GETLOGICAL('CalcReacRates','.FALSE.')

IF(CalcReacRates) THEN
  IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    CALL abort(&
      __STAMP__&
      ,'ERROR: CalcReacRates is not supported with radial weighting or variable time step yet!')
  END IF
END IF

IF(CalcNumSpec.OR.CalcNumDens.OR.CalcCollRates.OR.CalcReacRates) DoPartAnalyze = .TRUE.
! compute transversal or thermal velocity of whole computational domain
CalcVelos = GETLOGICAL('CalcVelos','.FALSE')
IF (CalcVelos) THEN
  IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    CALL abort(&
      __STAMP__&
      ,'ERROR: CalcVelos is not supported with radial weighting or variable time step yet!')
  END IF
  DoPartAnalyze=.TRUE.
  VeloDirs_hilf = GetIntArray('VelocityDirections',4,'1,1,1,1') ! x,y,z,abs -> 0/1 = T/F
  VeloDirs(:) = .FALSE.
  IF(.NOT.CalcNumSpec)THEN
    SWRITE(UNIT_stdOut,'(A)') ' Velocity computation requires NumSpec and SimNumSpec. Setting CalcNumSpec=.TRUE.'
    CalcNumSpec = .TRUE.
  END IF
  DO dir = 1,4
    IF (VeloDirs_hilf(dir) .EQ. 1) THEN
      VeloDirs(dir) = .TRUE.
    END IF
  END DO
  IF ((.NOT. VeloDirs(1)) .AND. (.NOT. VeloDirs(2)) .AND. &
      (.NOT. VeloDirs(3)) .AND. (.NOT. VeloDirs(4))) THEN
    CALL abort(&
      __STAMP__&
      ,'No VelocityDirections set in CalcVelos!')
  END IF
END IF
! Shape function efficiency
CalcShapeEfficiency = GETLOGICAL('CalcShapeEfficiency','.FALSE.')
IF (CalcShapeEfficiency) THEN
  DoPartAnalyze = .TRUE.
  CalcShapeEfficiencyMethod = GETSTR('CalcShapeEfficiencyMethod','AllParts')
  SELECT CASE(CalcShapeEfficiencyMethod)
  CASE('AllParts')  ! All currently available Particles are used
  CASE('SomeParts') ! A certain percentage of currently available Particles is used
    ShapeEfficiencyNumber = GETINT('ShapeEfficiencyNumber','100')  ! in percent
  CASE DEFAULT
    CALL abort(&
        __STAMP__&
        , ' CalcShapeEfficiencyMethod not implemented: ')

  END SELECT
END IF
! check if total energy should be computed
IF(DoPartAnalyze)THEN
  CalcEtot = GETLOGICAL('CalcTotalEnergy','.FALSE.')
END IF
IsRestart = GETLOGICAL('IsRestart','.FALSE.')

! Output for porous BC: Pump averaged values
CalcPorousBCInfo = GETLOGICAL('CalcPorousBCInfo','.FALSE.')
IF(CalcPorousBCInfo) DoPartAnalyze = .TRUE.

ParticleAnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTCILE ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleAnalyze


SUBROUTINE AnalyzeParticles(Time)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
! MPI-INFO:
! important change: all MPI-communication is done directly in the corresponding subroutines
! a) its easier and cleaner
! b) reduces the probability of errors (routines are again fixed for MPI...)
! The decision if an analysis is performed is done in PerformAnalysis.
! Furthermore:
! Each separate outputfile handler is called from within PerformAnalysis
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: CollInf, useDSMC, CollisMode, ChemReac
USE MOD_Restart_Vars           ,ONLY: DoRestart
USE MOD_Analyze_Vars           ,ONLY: CalcEpot,Wel,Wmag
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_TimeDisc_Vars          ,ONLY: iter, dt
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
USE MOD_DSMC_Analyze           ,ONLY: CalcMeanFreePath
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, BGGas
USE MOD_Particle_Vars          ,ONLY: Species
#endif
USE MOD_PIC_Analyze            ,ONLY: CalcDepositedCharge
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#if ( PP_TimeDiscMethod ==42)
#endif
USE MOD_Particle_Analyze_Vars  ,ONLY: ChemEnergySum
#ifdef CODE_ANALYZE
USE MOD_Analyze_Vars           ,ONLY: OutputErrorNorms
#endif /* CODE_ANALYZE */
USE MOD_Particle_Boundary_Vars, ONLY: nPorousBC, PorousBC
USE MOD_FPFlow_Vars            ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_FPFlow_Vars            ,ONLY: FP_PrandtlNumber, FPInitDone
USE MOD_BGK_Vars               ,ONLY: BGK_MaxRelaxFactor, BGK_MaxRotRelaxFactor, BGK_MeanRelaxFactor, BGK_MeanRelaxFactorCounter
USE MOD_BGK_Vars               ,ONLY: BGKInitDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: isOpen
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index, iSpec, OutputCounter, iPBC
INTEGER(KIND=8)     :: SimNumSpec(nSpecAnalyze)
REAL                :: NumSpec(nSpecAnalyze), NumDens(nSpecAnalyze)
REAL                :: Ekin(nSpecAnalyze), Temp(nSpecAnalyze)
REAL                :: EkinMax(nSpecies)
REAL                :: IntEn(nSpecAnalyze,3),IntTemp(nSpecies,3),TempTotal(nSpecAnalyze), Xi_Vib(nSpecies), Xi_Elec(nSpecies)
REAL                :: ETotal, totalChemEnergySum
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
REAL                :: MaxCollProb, MeanCollProb, MeanFreePath
REAL                :: NumSpecTmp(nSpecAnalyze)
#endif
#ifdef MPI
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
REAL                :: sumMeanCollProb
#endif
REAL                :: RECBR(nSpecies),RECBR1
INTEGER             :: RECBIM(nSpecies)
#endif /*MPI*/
REAL, ALLOCATABLE   :: CRate(:), RRate(:)
#if (PP_TimeDiscMethod ==42)
INTEGER             :: iCase, iTvib,jSpec
#ifdef CODE_ANALYZE
CHARACTER(LEN=64)   :: DebugElectronicStateFilename
INTEGER             :: ii, iunit
#endif
CHARACTER(LEN=350)  :: hilf
#endif
REAL                :: PartVtrans(nSpecies,4) ! macroscopic velocity (drift velocity) A. Frohn: kinetische Gastheorie
REAL                :: PartVtherm(nSpecies,4) ! microscopic velocity (eigen velocity) PartVtrans + PartVtherm = PartVtotal
INTEGER             :: dir
!===================================================================================================================================
  IF ( DoRestart ) THEN
    isRestart = .true.
  END IF
  IF (.NOT.DoPartAnalyze) RETURN
  IF (useDSMC) THEN
    IF (CollisMode.NE.0) THEN
      SDEALLOCATE(CRate)
      ALLOCATE(CRate(CollInf%NumCase + 1))
      IF (CollisMode.EQ.3) THEN
        SDEALLOCATE(RRate)
        ALLOCATE(RRate(ChemReac%NumOfReact))
        RRate = 0.0
      END IF
    END IF
  END IF
  OutputCounter = 2
  unit_index = 535
  IF (PartMPI%MPIRoot) THEN
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
#if (PP_TimeDiscMethod==42)
    ! if only the reaction rate is desired (resevoir) the initial temperature
    ! of the second species is added to the filename
      IF (DSMC%ReservoirSimuRate) THEN
        IF ( SpecDSMC(1)%InterID .EQ. 2 .OR. SpecDSMC(1)%InterID .EQ. 20 ) THEN
          iTvib = INT(SpecDSMC(1)%Init(0)%Tvib)
          WRITE( hilf, '(I5.5)') iTvib
          outfile = 'Database_Tvib_'//TRIM(hilf)//'.csv'
        ELSE
          !iTvib = INT(SpecDSMC(1)%Telec )
          iTvib = INT(Species(1)%Init(0)%MWTemperatureIC) !wrong name, if MWTemp is defined in %Init!!!
          WRITE( hilf, '(I5.5)') iTvib
          outfile = 'Database_Ttrans_'//TRIM(hilf)//'.csv'
        END IF
      ELSE
        outfile = 'PartAnalyze.csv'
      END IF
#else
      outfile = 'PartAnalyze.csv'
#endif

      IF (isRestart .and. FILEEXISTS(outfile)) THEN
        OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
        !CALL FLUSH (unit_index)
      ELSE
        OPEN(unit_index,file=TRIM(outfile))
        !CALL FLUSH (unit_index)
        !--- insert header
        WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
        IF (CalcNumSpec) THEN
          DO iSpec = 1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A12,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPart-Spec-', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcNumDens) THEN
          DO iSpec = 1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A12,I3.3,A5)',ADVANCE='NO') OutputCounter,'-NumDens-Spec-', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcCharge) THEN
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A12,A5)',ADVANCE='NO') OutputCounter,'-Charge',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A16,A5)',ADVANCE='NO') OutputCounter,'-Charge-absError',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A16,A5)',ADVANCE='NO') OutputCounter,'-Charge-relError',' '
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcPartBalance) THEN
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A14,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartIn-Spec-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A15,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartOut-Spec-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcEkin) THEN
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Ekin-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcCoupledPower) THEN
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-PCoupled',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-PCoupledMoAv',' '
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcLaserInteraction) THEN
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinMax-eV-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF(CalcEkin .AND. CalcEpot .AND. CalcEtot) THEN
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-E-kin+pot',' '
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcTemp) THEN
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempTra-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcVelos) THEN
          DO iSpec=1, nSpecies
            IF (VeloDirs(1)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Velo_Xtrans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Velo_Xtherm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
            IF (VeloDirs(2)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Velo_Ytrans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Velo_Ytherm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
            IF (VeloDirs(3)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Velo_Ztrans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Velo_Ztherm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
            IF (VeloDirs(4)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-AbsVelo_trans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-AbsVelo_therm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
          END DO
        END IF
        IF (CalcPartBalance) THEN
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A8,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinIn-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A9,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinOut-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
#if (PP_TimeDiscMethod==1000)
        IF (CollisMode.GT.1) THEN
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EVib',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-ERot',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          IF ( DSMC%ElectronicModel ) THEN
            DO iSpec = 1, nSpecAnalyze
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EElec',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempVib',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempRot',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
#endif
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
        IF (CollisMode.GT.1) THEN
          IF(CalcEint) THEN
            DO iSpec=1, nSpecAnalyze
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-E-Vib',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecAnalyze
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-E-Rot',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            IF (DSMC%ElectronicModel) THEN
              DO iSpec = 1, nSpecAnalyze
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-E-Elec',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
            END IF
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-E-TotalPart',' '
            OutputCounter = OutputCounter + 1
          END IF
          IF(CalcEpot .AND. CalcEtot .AND. CalcEint)THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-E-Tot',' '
            OutputCounter = OutputCounter + 1
          END IF
          IF(CalcTemp) THEN
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempVib',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-XiVibMean',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempRot',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            IF ( DSMC%ElectronicModel ) THEN
              DO iSpec=1, nSpecies
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempElec',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iSpec=1, nSpecies
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-XiElecMean',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
            END IF
            DO iSpec=1, nSpecAnalyze
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempTotal',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
        END IF
        IF(CalcPorousBCInfo) THEN
          DO iPBC = 1, nPorousBC
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I2.2,A,A5)',ADVANCE='NO') OutputCounter,'-PorousBC-',iPBC,'-PumpSpeed-Measure',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I2.2,A,A5)',ADVANCE='NO') OutputCounter,'-PorousBC-',iPBC,'-PumpSpeed-Control',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I2.2,A,A5)',ADVANCE='NO') OutputCounter,'-PorousBC-',iPBC,'-RemovalProbability',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I2.2,A,A5)',ADVANCE='NO') OutputCounter,'-PorousBC-',iPBC,'-PressureNorm',' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF(DSMC%CalcQualityFactors) THEN
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-Pmean',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-Pmax',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-MeanFreePath',' '
          OutputCounter = OutputCounter + 1
        END IF
#endif
        IF(FPInitDone) THEN
          IF(DSMC%CalcQualityFactors) THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-FP-MeanRelaxFactor',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-FP-MaxRelaxFactor',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-FP-MaxRotRelaxFactor',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-FP-MeanPrandtlNumber',' '
            OutputCounter = OutputCounter + 1
          END IF
        END IF
        IF(BGKInitDone) THEN
          IF(DSMC%CalcQualityFactors) THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-BGK-MeanRelaxFactor',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-BGK-MaxRelaxFactor',' '
            OutputCounter = OutputCounter + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-BGK-MaxRotRelaxFactor',' '
            OutputCounter = OutputCounter + 1
          END IF
        END IF
#if (PP_TimeDiscMethod==42)
        IF(CalcCollRates) THEN
          DO iSpec = 1, nSpecies
            DO jSpec = iSpec, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-CollRate', iSpec, '+', jSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END DO
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-TotalCollRate',' '
          OutputCounter = OutputCounter + 1
        END IF
        IF(CalcReacRates) THEN
          IF(CollisMode.EQ.3) THEN
            DO iCase=1, ChemReac%NumOfReact
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Reaction', iCase,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
        END IF
#endif
        WRITE(unit_index,'(A1)') ' '
      END IF
    END IF
  END IF

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================
  ! computes the real and simulated number of particles
  CALL CalcNumPartsOfSpec(NumSpec,SimNumSpec)
  IF(CalcNumDens) CALL CalcNumberDensity(NumSpec,NumDens)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate total temperature of each molecular species (Laux, p. 109)
  IF(CalcEkin.OR.CalcEint)THEN
    IF(CalcLaserInteraction)THEN
      CALL CalcKineticEnergyAndMaximum(Ekin,EkinMax)
    ELSE
      CALL CalcKineticEnergy(Ekin)
    END IF
  END IF
  IF(CalcTemp.OR.CalcEint.OR.DSMC%CalcQualityFactors) THEN
    CALL CalcTemperature(NumSpec,Temp,IntTemp,IntEn,TempTotal,Xi_Vib,Xi_Elec) ! contains MPI Communication
    IF(CalcEint.AND.(CollisMode.GT.1)) THEN
      ETotal = Ekin(nSpecAnalyze) + IntEn(nSpecAnalyze,1) + IntEn(nSpecAnalyze,2) + IntEn(nSpecAnalyze,3)
      IF(CollisMode.EQ.3) THEN
        totalChemEnergySum = 0.
#ifdef MPI
        IF(PartMPI%MPIRoot) THEN
          CALL MPI_REDUCE(ChemEnergySum,totalChemEnergySum,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        ELSE
          CALL MPI_REDUCE(ChemEnergySum,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        END IF
#else
        totalChemEnergySum = ChemEnergySum
#endif
        ETotal = ETotal - totalChemEnergySum
      END IF
    END IF
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Determine the maximal collision probability for whole reservoir and mean collision probability (only for one cell)
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
  MaxCollProb = 0.0
  MeanCollProb = 0.0
  MeanFreePath = 0.0
  IF(DSMC%CalcQualityFactors.OR.CalcReacRates) THEN
    NumSpecTmp = NumSpec
    IF(BGGas%BGGasSpecies.NE.0) THEN
      ! Calculation of mean free path and reactions rates requires the number of particles the background species would have if
      ! actually inserted at the chosen weighting factor, determined here and used later also for the ReacRates subroutine
      NumSpecTmp(BGGas%BGGasSpecies) = (BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor)
      IF(nSpecAnalyze.GT.1)THEN
        NumSpecTmp(nSpecAnalyze) = NumSpecTmp(nSpecAnalyze)+NumSpecTmp(BGGas%BGGasSpecies)
      END IF
    END IF
  END IF
  IF(DSMC%CalcQualityFactors) THEN
    IF(iter.GT.0) THEN
      MaxCollProb = DSMC%CollProbMax
      IF(DSMC%CollProbMeanCount.GT.0) MeanCollProb = DSMC%CollProbMean / DSMC%CollProbMeanCount
      IF (PartMPI%MPIRoot) THEN
        IF(TempTotal(nSpecAnalyze).GT.0.0) MeanFreePath = CalcMeanFreePath(NumSpecTmp(1:nSpecies), NumSpecTmp(nSpecAnalyze), &
                                                              GEO%MeshVolume, SpecDSMC(1)%omegaVHS, TempTotal(nSpecAnalyze))
      END IF
    END IF
  END IF
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! Other Analyze Routines
  IF(CalcCharge) CALL CalcDepositedCharge() ! mpi communication done in calcdepositedcharge
  ! get velocities
  IF(CalcVelos) CALL CalcVelocities(PartVtrans, PartVtherm,NumSpec,SimNumSpec)
!===================================================================================================================================
! MPI Communication for values which are not YET communicated
! All routines ABOVE contain the required MPI-Communication
!===================================================================================================================================
#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    IF (CalcPartBalance)THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,nPartIn(1:nSpecAnalyze)    ,nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,nPartOut(1:nSpecAnalyze)   ,nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinIn(1:nSpecAnalyze) ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinOut(1:nSpecAnalyze),nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
    IF(CalcCoupledPower) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,PCoupl,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PCouplAverage,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
    IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
      ! Determining the maximal (MPI_MAX) and mean (MPI_SUM) collision probabilities
      CALL MPI_REDUCE(MPI_IN_PLACE,MaxCollProb,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      CALL MPI_REDUCE(MeanCollProb,sumMeanCollProb,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
      MeanCollProb = sumMeanCollProb / REAL(PartMPI%nProcs)
    END IF
#endif
    IF(CalcPorousBCInfo) THEN
      DO iPBC = 1, nPorousBC
        CALL MPI_REDUCE(MPI_IN_PLACE,PorousBC(iPBC)%Output,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
      END DO
    END IF
    IF(FPInitDone) THEN
      IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,FP_MeanRelaxFactor,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(MPI_IN_PLACE,FP_MeanRelaxFactorCounter,1, MPI_INTEGER, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(MPI_IN_PLACE,FP_PrandtlNumber,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        ! Determining the maximal (MPI_MAX) relaxation factors
        CALL MPI_REDUCE(MPI_IN_PLACE,FP_MaxRelaxFactor,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(MPI_IN_PLACE,FP_MaxRotRelaxFactor,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      END IF
    END IF
    IF(BGKInitDone) THEN
      IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,BGK_MeanRelaxFactor,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(MPI_IN_PLACE,BGK_MeanRelaxFactorCounter,1, MPI_INTEGER, MPI_SUM,0, PartMPI%COMM, IERROR)
        ! Determining the maximal (MPI_MAX) relaxation factors
        CALL MPI_REDUCE(MPI_IN_PLACE,BGK_MaxRelaxFactor,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(MPI_IN_PLACE,BGK_MaxRotRelaxFactor,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      END IF
    END IF
  ELSE ! no Root
    IF (CalcPartBalance)THEN
      CALL MPI_REDUCE(nPartIn,RECBIM   ,nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(nPartOut,RECBIM  ,nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PartEkinIn,RECBR ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PartEkinOut,RECBR,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
    IF(CalcCoupledPower) THEN
      CALL MPI_REDUCE(PCoupl,PCoupl,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PCouplAverage,PCouplAverage,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
    IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
      CALL MPI_REDUCE(MaxCollProb,RECBR1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, PartMPI%COMM, IERROR)
      CALL MPI_REDUCE(MeanCollProb,sumMeanCollProb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, PartMPI%COMM, IERROR)
    END IF
#endif
    IF(CalcPorousBCInfo) THEN
      DO iPBC = 1, nPorousBC
        CALL MPI_REDUCE(PorousBC(iPBC)%Output,PorousBC(iPBC)%Output,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
      END DO
    END IF
    IF(FPInitDone) THEN
      IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
        CALL MPI_REDUCE(FP_MeanRelaxFactor,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(FP_MeanRelaxFactorCounter,RECBR1,1, MPI_INTEGER, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(FP_PrandtlNumber,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(FP_MaxRelaxFactor,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(FP_MaxRotRelaxFactor,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      END IF
    END IF
    IF(BGKInitDone) THEN
      IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
        CALL MPI_REDUCE(BGK_MeanRelaxFactor,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(BGK_MeanRelaxFactorCounter,RECBR1,1, MPI_INTEGER, MPI_SUM,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(BGK_MaxRelaxFactor,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
        CALL MPI_REDUCE(BGK_MaxRotRelaxFactor,RECBR1,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      END IF
    END IF
  END IF
#endif /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! Perform averaging/summation of the MPI communicated variables on the root only (and for the non-MPI case, MPIRoot is set to true)
IF(PartMPI%MPIRoot) THEN
  IF (CalcPartBalance)THEN
    IF(nSpecies.GT.1) THEN
      nPartIn(nSpecies+1)     = SUM(nPartIn(1:nSpecies))
      nPartOut(nSpecies+1)    = SUM(nPartOut(1:nSpecies))
      PartEkinIn(nSpecies+1)  = SUM(PartEkinIn(1:nSpecies))
      PartEkinOut(nSpecies+1) = SUM(PartEkinOut(1:nSpecies))
    END IF
  END IF
  IF(CalcPorousBCInfo) THEN
    DO iPBC = 1, nPorousBC
      IF(PorousBC(iPBC)%Output(1).GT.0.0) THEN
        ! Pumping Speed (Output(2)) is the sum of all elements (counter over particles exiting through pump)
        ! Other variales are averaged over the elements
        PorousBC(iPBC)%Output(3:5) = PorousBC(iPBC)%Output(3:5) / PorousBC(iPBC)%Output(1)
      END IF
    END DO
  END IF
  IF(FPInitDone) THEN
    IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
      IF(FP_MeanRelaxFactorCounter.GT.0) THEN
        FP_MeanRelaxFactor = FP_MeanRelaxFactor / REAL(FP_MeanRelaxFactorCounter)
        FP_PrandtlNumber = FP_PrandtlNumber / REAL(FP_MeanRelaxFactorCounter)
      END IF
    END IF
  END IF
  IF(BGKInitDone) THEN
    IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
      IF(BGK_MeanRelaxFactorCounter.GT.0) BGK_MeanRelaxFactor = BGK_MeanRelaxFactor / REAL(BGK_MeanRelaxFactorCounter)
    END IF
  END IF
  IF(CalcCoupledPower) THEN
  ! Moving Average of PCoupl:
    IF(iter.EQ.0) THEN
      PCouplAverage = 0.0
    ELSE
      PCouplAverage = PCouplAverage / Time
    END IF
    ! current PCoupl (Delta_E / Timestep)
    PCoupl = PCoupl / dt
  END IF
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==1000)
  IF (CollisMode.GT.1) CALL CalcIntTempsAndEn(NumSpec,IntTemp,IntEn)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate the collision rates and reaction rate coefficients (Arrhenius-type chemistry)
#if (PP_TimeDiscMethod==42)
  IF(CalcCollRates) CALL CollRates(CRate)
  IF(CalcReacRates) THEN
    IF ((CollisMode.EQ.3).AND.(iter.GT.0)) THEN
      CALL ReacRates(RRate, NumSpecTmp, iter)
    END IF
  END IF
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (CalcShapeEfficiency) CALL CalcShapeEfficiencyR()   ! This will NOT be placed in the file but directly in "out"
!===================================================================================================================================
! Output Routines
!===================================================================================================================================
#ifdef MPI
IF (PartMPI%MPIROOT) THEN
#endif    /* MPI */
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Time
    IF (CalcNumSpec) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(SimNumSpec(iSpec))
      END DO
    END IF
    IF (CalcNumDens) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(NumDens(iSpec))
      END DO
    END IF
    IF (CalcCharge) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartCharge(1)
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartCharge(2)
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartCharge(3)
    END IF
    IF (CalcPartBalance) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartIn(iSpec))
      END DO
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartOut(iSpec))
      END DO
    END IF
    IF (CalcEkin) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Ekin(iSpec)
      END DO
    END IF
    IF (CalcCoupledPower) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PCoupl
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PCouplAverage
    END IF
    IF (CalcLaserInteraction) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') EkinMax(iSpec)
      END DO
    END IF
    IF (CalcEpot .AND. CalcEkin .AND. CalcEtot) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Ekin(nSpecAnalyze) + WEl + WMag
    END IF
    IF (CalcTemp) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Temp(iSpec)
      END DO
    END IF
    IF (CalcVelos) THEN
      DO iSpec=1, nSpecies
        DO dir = 1,4
          IF (VeloDirs(dir)) THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartVtrans(iSpec,dir)
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartVtherm(iSpec,dir)
          END IF
        END DO
      END DO
    END IF
    IF (CalcPartBalance) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinIn(iSpec)
      END DO
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinOut(iSpec)
      END DO
    END IF

#if (PP_TimeDiscMethod==1000)
    IF (CollisMode.GT.1) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntEn(iSpec,1)
      END DO
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntEn(iSpec,2)
      END DO
      IF ( DSMC%ElectronicModel ) THEN
        DO iSpec=1, nSpecAnalyze
        ! currently set to one
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntEn(iSpec,3)
        END DO
      END IF
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntTemp(iSpec,1)
      END DO
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntTemp(iSpec,2)
      END DO
    END IF
#endif
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509))
    IF (CollisMode.GT.1) THEN
      IF(CalcEint) THEN
        DO iSpec=1, nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntEn(iSpec,1)
        END DO
        DO iSpec=1, nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntEn(iSpec,2)
        END DO
        IF (DSMC%ElectronicModel) THEN
          DO iSpec=1, nSpecAnalyze
          ! currently set to one
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntEn(iSpec,3)
          END DO
        END IF
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') ETotal
      END IF
      IF(CalcEpot .AND. CalcEtot .AND. CalcEint)THEN
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') ETotal+WEl+WMag
      END IF
      IF(CalcTemp) THEN
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntTemp(iSpec,1)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Xi_Vib(iSpec)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntTemp(iSpec,2)
        END DO
        IF ( DSMC%ElectronicModel ) THEN
          DO iSpec=1, nSpecies
          ! currently set to one
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') IntTemp(iSpec,3)
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Xi_Elec(iSpec)
          END DO
        END IF
        DO iSpec=1, nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') TempTotal(iSpec)
        END DO
      END IF
    END IF
    IF(CalcPorousBCInfo) THEN
      DO iPBC = 1, nPorousBC
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(2)
        OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(3)
        OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(4)
        OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(5)
        OutputCounter = OutputCounter + 1
      END DO
    END IF
    IF(DSMC%CalcQualityFactors) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') MeanCollProb
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') MaxCollProb
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') MeanFreePath
    END IF
#endif
    IF(FPInitDone) THEN
      IF(DSMC%CalcQualityFactors) THEN
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') FP_MeanRelaxFactor
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') FP_MaxRelaxFactor
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') FP_MaxRotRelaxFactor
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') FP_PrandtlNumber
      END IF
    END IF
    IF(BGKInitDone) THEN
      IF(DSMC%CalcQualityFactors) THEN
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') BGK_MeanRelaxFactor
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') BGK_MaxRelaxFactor
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') BGK_MaxRotRelaxFactor
      END IF
    END IF
#if (PP_TimeDiscMethod==42)
    IF(CalcCollRates) THEN
      DO iCase=1, CollInf%NumCase +1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') CRate(iCase)
      END DO
    END IF
    IF(CalcReacRates) THEN
      DO iCase=1, ChemReac%NumOfReact
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') RRate(iCase)
      END DO
    END IF
#endif /*(PP_TimeDiscMethod==42)*/
    WRITE(unit_index,'(A1)') ' '
#ifdef MPI
  END IF
#endif    /* MPI */

! Reset output variables
IF(CalcPorousBCInfo) THEN
  DO iPBC = 1,nPorousBC
    PorousBC(iPBC)%Output(1:5) = 0.
  END DO
END IF
IF (CalcCoupledPower) THEN                         ! if output of coupled power is active
  PCouplAverage = PCouplAverage * Time           ! PCouplAverage is reseted
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
  IF( CalcPartBalance) CALL CalcParticleBalance()
!-----------------------------------------------------------------------------------------------------------------------------------
#if ( PP_TimeDiscMethod ==42 )
#ifdef CODE_ANALYZE
IF (DSMC%ElectronicModel.AND.DSMC%ReservoirSimuRate) THEN
  ! Debug output for initialized electronic state
  IF(Time.GT.0.) CALL ElectronicTransition( Time, NumSpec )
  DO iSpec = 1, nSpecies
    IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      IF (  SpecDSMC(iSpec)%levelcounter(0) .ne. 0) THEN
        WRITE(DebugElectronicStateFilename,'(I2.2)') iSpec
        iunit = 485
        DebugElectronicStateFilename = 'End_Electronic_State_Species_'//trim(DebugElectronicStateFilename)//'.dat'
        OPEN(unit=iunit,file=DebugElectronicStateFilename,form='formatted',status='unknown')
        DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant - 1                         !has to be changed when using %Init definitions!!!
          WRITE(iunit,'(I3.1,3x,F12.7)') ii, REAL( SpecDSMC(iSpec)%levelcounter(ii) ) / &
                                              REAL( Species(iSpec)%Init(0)%initialParticleNumber)
        END DO
        CLOSE(iunit)
      END IF
    END IF
  END DO
END IF
#endif
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE AnalyzeParticles

! all other analysis with particles
SUBROUTINE CalcShapeEfficiencyR()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars ,ONLY: CalcShapeEfficiencyMethod, ShapeEfficiencyNumber
USE MOD_Mesh_Vars             ,ONLY: nElems, Elem_xGP
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_PreProc
#ifdef MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: NbrOfComps, NbrWithinRadius, NbrOfElems, NbrOfElemsWithinRadius
REAL                     :: RandVal1
LOGICAL                  :: chargedone(1:nElems), WITHIN
INTEGER                  :: kmin, kmax, lmin, lmax, mmin, mmax                           !
INTEGER                  :: kk, ll, mm, ppp,m,l,k, i                                             !
INTEGER                  :: ElemID
REAL                     :: radius, deltax, deltay, deltaz
!===================================================================================================================================

NbrOfComps = 0.
NbrOfElems = 0.
NbrWithinRadius = 0.
NbrOfElemsWithinRadius = 0.
SELECT CASE(CalcShapeEfficiencyMethod)
CASE('AllParts')
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      chargedone(:) = .FALSE.
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      kmax = INT((PartState(i,1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = INT((PartState(i,1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      lmax = INT((PartState(i,2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = INT((PartState(i,2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      mmax = INT((PartState(i,3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = INT((PartState(i,3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      !-- go through all these cells
      DO kk = kmin,kmax
        DO ll = lmin, lmax
          DO mm = mmin, mmax
            !--- go through all mapped elements not done yet
            DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              WITHIN=.FALSE.
              ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
              IF (.NOT.chargedone(ElemID)) THEN
                NbrOfElems = NbrOfElems + 1.
                !--- go through all gauss points
                DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                  NbrOfComps = NbrOfComps + 1.
                  !-- calculate distance between gauss and particle
                  deltax = PartState(i,1) - Elem_xGP(1,k,l,m,ElemID)
                  deltay = PartState(i,2) - Elem_xGP(2,k,l,m,ElemID)
                  deltaz = PartState(i,3) - Elem_xGP(3,k,l,m,ElemID)
                  radius = deltax * deltax + deltay * deltay + deltaz * deltaz
                  IF (radius .LT. r2_sf) THEN
                    WITHIN=.TRUE.
                    NbrWithinRadius = NbrWithinRadius + 1.
                  END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
              IF(WITHIN) NbrOfElemsWithinRadius = NbrOfElemsWithinRadius + 1.
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END IF ! inside
  END DO ! i
IF(NbrOfComps.GT.0.0)THEN
#ifdef MPI
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',PartMPI%MyRank,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',PartMPI%MyRank,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps, 100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
CASE('SomeParts')
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      CALL RANDOM_NUMBER(RandVal1)
      IF(RandVal1.LT.REAL(ShapeEfficiencyNumber)/100)THEN
        chargedone(:) = .FALSE.
        !-- determine which background mesh cells (and interpolation points within) need to be considered
        kmax = INT((PartState(i,1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = INT((PartState(i,1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = INT((PartState(i,2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = INT((PartState(i,2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = INT((PartState(i,3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = INT((PartState(i,3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              WITHIN=.FALSE.
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF (.NOT.chargedone(ElemID)) THEN
                  NbrOfElems = NbrOfElems + 1
                  !--- go through all gauss points
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    NbrOfComps = NbrOfComps + 1
                    !-- calculate distance between gauss and particle
                    deltax = PartState(i,1) - Elem_xGP(1,k,l,m,ElemID)
                    deltay = PartState(i,2) - Elem_xGP(2,k,l,m,ElemID)
                    deltaz = PartState(i,3) - Elem_xGP(3,k,l,m,ElemID)
                    radius = deltax * deltax + deltay * deltay + deltaz * deltaz
                    IF (radius .LT. r2_sf) THEN
                      NbrWithinRadius = NbrWithinRadius + 1
                      WITHIN=.TRUE.
                    END IF
                    END DO; END DO; END DO
                    chargedone(ElemID) = .TRUE.
                  END IF
                  IF(WITHIN) NbrOfElemsWithinRadius = NbrOfElemsWithinRadius + 1
                END DO ! ppp
              END DO ! mm
            END DO ! ll
          END DO ! kk
        END IF  ! RandVal
      END IF ! inside
  END DO ! i
IF(NbrOfComps.GT.0)THEN
#ifdef MPI
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',PartMPI%MyRank,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',PartMPI%MyRank,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
END SELECT
END SUBROUTINE CalcShapeEfficiencyR


SUBROUTINE CalcParticleBalance()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars,      ONLY : nPartIn,nPartOut,PartEkinIn,PartEkinOut
#if defined(LSERK) || defined(ROS) || defined(IMPA)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
USE MOD_Particle_Analyze_Vars ,ONLY: nPartInTmp,PartEkinInTmp
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if defined(LSERK) || defined(ROS) || defined(IMPA)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
nPartIn=nPartInTmp
nPartOut=0
PartEkinIn=PartEkinInTmp
PartEkinOut=0.
nPartInTmp=0
PartEkinInTmp=0.
#else
nPartIn=0
nPartOut=0
PartEkinIn=0.
PartEkinOut=0.
#endif

END SUBROUTINE CalcParticleBalance


SUBROUTINE CalcKineticEnergy(Ekin)
!===================================================================================================================================
! compute the kinetic energy of particles
! for velocity <1e3 non-relativistic formula is used, for larger velocities the relativistic kinetic energy is computed
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars         ,ONLY: c2, c2_inv
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, Species, PDM
USE MOD_PARTICLE_Vars         ,ONLY: usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#ifndef PP_HDG
USE MOD_PML_Vars              ,ONLY: DoPML,xyzPhysicalMinMax
#endif /*PP_HDG*/
#ifdef MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i
REAL(KIND=8)                    :: partV2, GammaFac
REAL                            :: Ekin_loc
!===================================================================================================================================
Ekin    = 0.!d0
IF (nSpecAnalyze.GT.1) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
#ifndef PP_HDG
      IF(DoPML)THEN
        IF (PartState(i,1) .GE. xyzPhysicalMinMax(1) .AND. PartState(i,1) .LE. xyzPhysicalMinMax(2) .AND. &
            PartState(i,2) .GE. xyzPhysicalMinMax(3) .AND. PartState(i,2) .LE. xyzPhysicalMinMax(4) .AND. &
            PartState(i,3) .GE. xyzPhysicalMinMax(5) .AND. PartState(i,3) .LE. xyzPhysicalMinMax(6)) THEN
          CYCLE
        END IF
      ENDIF
#endif /*PP_HDG*/
      partV2 = PartState(i,4) * PartState(i,4) &
              + PartState(i,5) * PartState(i,5) &
              + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000000
        Ekin_loc = 0.5 * Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! %MacroParticleFactor is included in the case of RadialWeighting (also in combination with variable time step)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          ! Case for variable time step without radial weighting (no regular weighting factor applied in GetParticleWeight)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF != usevMPF
      ELSE ! partV2 > 1e6
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! %MacroParticleFactor is included in the case of RadialWeighting (also in combination with variable time step)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          ! Case for variable time step without radial weighting (no regular weighting factor applied in GetParticleWeight)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF != usevMPF
      END IF ! partV2
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
ELSE ! nSpecAnalyze = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
#ifndef PP_HDG
      IF(DoPML)THEN
        IF (PartState(i,1) .GE. xyzPhysicalMinMax(1) .AND. PartState(i,1) .LE. xyzPhysicalMinMax(2) .AND. &
            PartState(i,2) .GE. xyzPhysicalMinMax(3) .AND. PartState(i,2) .LE. xyzPhysicalMinMax(4) .AND. &
            PartState(i,3) .GE. xyzPhysicalMinMax(5) .AND. PartState(i,3) .LE. xyzPhysicalMinMax(6)) THEN
          CYCLE
        END IF
      ENDIF
#endif /*PP_HDG*/
      partV2 = PartState(i,4) * PartState(i,4) &
             + PartState(i,5) * PartState(i,5) &
             + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000000
        Ekin_loc = 0.5 *  Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! usevMPF
      ELSE ! partV2 > 1e6
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting)THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! useuvMPF
      END IF ! par2
    END IF ! particle inside
  END DO ! particleveclength
END IF

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE , Ekin    , nSpecAnalyze , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , PartMPI%COMM , IERROR)
ELSE
  CALL MPI_REDUCE(Ekin         , 0.      , nSpecAnalyze , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , PartMPI%COMM , IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergy


SUBROUTINE CalcKineticEnergyAndMaximum(Ekin,EkinMax)
!===================================================================================================================================
! compute the kinetic energy of particles
! for velocity <1e3 non-relativistic formula is used, for larger velocities the relativistic kinetic energy is computed
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars         ,ONLY: c2, c2_inv
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, Species, PDM, nSpecies
USE MOD_PARTICLE_Vars         ,ONLY: usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze,LaserInteractionEkinMaxRadius,LaserInteractionEkinMaxZPosMin
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#ifndef PP_HDG
USE MOD_PML_Vars              ,ONLY: DoPML,xyzPhysicalMinMax
#endif /*PP_HDG*/
#ifdef MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(nSpecAnalyze)
REAL,INTENT(OUT)                :: EkinMax(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i
REAL(KIND=8)                    :: partV2, GammaFac
REAL                            :: Ekin_loc
!===================================================================================================================================
! default values
Ekin    =  0.
EkinMax = -1.
! set boundaries in order to exclude particles near the boundary (nonphysical velocities)
IF (nSpecAnalyze.GT.1) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
#ifndef PP_HDG
      IF(DoPML)THEN
        IF (PartState(i,1) .GE. xyzPhysicalMinMax(1) .AND. PartState(i,1) .LE. xyzPhysicalMinMax(2) .AND. &
            PartState(i,2) .GE. xyzPhysicalMinMax(3) .AND. PartState(i,2) .LE. xyzPhysicalMinMax(4) .AND. &
            PartState(i,3) .GE. xyzPhysicalMinMax(5) .AND. PartState(i,3) .LE. xyzPhysicalMinMax(6)) THEN
          CYCLE
        END IF
      ENDIF
#endif /*PP_HDG*/
      partV2 = PartState(i,4) * PartState(i,4) &
              + PartState(i,5) * PartState(i,5) &
              + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000000
        Ekin_loc = 0.5 * Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! %MacroParticleFactor is included in the case of RadialWeighting (also in combination with variable time step)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          ! Case for variable time step without radial weighting (no regular weighting factor applied in GetParticleWeight)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF != usevMPF
      ELSE ! partV2 > 1e6
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF !=usevMPF
      END IF ! partV2
      ! Determine energy of the most energetic particle in [eV]
      IF((SQRT(PartState(i,1)**2 + PartState(i,2)**2).LE.LaserInteractionEkinMaxRadius).OR.&
                                      (PartState(i,3).GE.LaserInteractionEkinMaxZPosMin))THEN
        EkinMax(PartSpecies(i)) = MAX(EkinMax(PartSpecies(i)),Ekin_loc*6.241509e18) ! 6.241509e18 is [J] -> [eV]
      END IF
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
ELSE ! nSpecAnalyze = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
#ifndef PP_HDG
      IF(DoPML)THEN
        IF (PartState(i,1) .GE. xyzPhysicalMinMax(1) .AND. PartState(i,1) .LE. xyzPhysicalMinMax(2) .AND. &
            PartState(i,2) .GE. xyzPhysicalMinMax(3) .AND. PartState(i,2) .LE. xyzPhysicalMinMax(4) .AND. &
            PartState(i,3) .GE. xyzPhysicalMinMax(5) .AND. PartState(i,3) .LE. xyzPhysicalMinMax(6)) THEN
          CYCLE
        END IF
      ENDIF
#endif /*PP_HDG*/
      partV2 = PartState(i,4) * PartState(i,4) &
             + PartState(i,5) * PartState(i,5) &
             + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000000
        Ekin_loc = 0.5 *  Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! usevMPF
      ELSE ! partV2 > 1e6
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting)THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! useuvMPF

      END IF ! par2
      ! Determine energy of the most energetic particle in [eV]
      IF((SQRT(PartState(i,1)**2 + PartState(i,2)**2).LE.LaserInteractionEkinMaxRadius).OR.&
                                      (PartState(i,3).GE.LaserInteractionEkinMaxZPosMin))THEN
        EkinMax(PartSpecies(i)) = MAX(EkinMax(PartSpecies(i)),Ekin_loc*6.241509e18) ! 6.241509e18 is [J] -> [eV]
      END IF
    END IF ! particle inside
  END DO ! particleveclength
END IF

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE , Ekin    , nSpecAnalyze , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , PartMPI%COMM , IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE , EkinMax , nSpecies     , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , PartMPI%COMM , iError)
ELSE
  CALL MPI_REDUCE(Ekin         , 0.      , nSpecAnalyze , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , PartMPI%COMM , IERROR)
  CALL MPI_REDUCE(EkinMax      , 0.      , nSpecies     , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , PartMPI%COMM , iError)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergyAndMaximum



SUBROUTINE CalcNumPartsOfSpec(NumSpec,SimNumSpec)
!===================================================================================================================================
! Computes the number of simulated particles AND number of real particles within the domain
! Last section of the routine contains the MPI-communication
! CAUTION: SimNumSpec equals NumSpec only for constant weighting factor
! NOTE: Background gas particles are not considered
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PDM,PartSpecies
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_DSMC_Vars             ,ONLY: BGGas
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#ifdef MPI
USE MOD_Particle_Analyze_Vars ,ONLY: CalcNumSpec
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: NumSpec(nSpecAnalyze)
INTEGER(KIND=8),INTENT(OUT)        :: SimNumSpec(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iPart
!===================================================================================================================================

NumSpec    = 0.
SimNumSpec = 0
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! NumSpec = real particle number
    NumSpec(PartSpecies(iPart))    = NumSpec(PartSpecies(iPart)) + GetParticleWeight(iPart)
    ! SimNumSpec = simulated particle number
    SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1
  END IF
END DO
IF(BGGas%BGGasSpecies.NE.0) THEN
  NumSpec(BGGas%BGGasSpecies) = 0.
  SimNumSpec(BGGas%BGGasSpecies) = 0
END IF
IF(nSpecAnalyze.GT.1)THEN
  NumSpec(nSpecAnalyze)    = SUM(NumSpec(1:nSpecies))
  SimNumSpec(nSpecAnalyze) = SUM(SimNumSpec(1:nSpecies))
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcNumSpec) &
  CALL MPI_REDUCE(MPI_IN_PLACE,SimNumSpec ,nSpecAnalyze,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(NumSpec     ,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcNumSpec) &
  CALL MPI_REDUCE(SimNumSpec  ,SimNumSpec ,nSpecAnalyze,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcNumPartsOfSpec


SUBROUTINE CalcNumberDensity(NumSpec,NumDens)
!===================================================================================================================================
!> Computes the number density per species using the total mesh volume and if neccessary particle weights
!> Background gas density is saved as given in the input
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: Species,usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_DSMC_Vars             ,ONLY: BGGas, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)                   :: NumSpec(nSpecAnalyze)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: NumDens(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (PartMPI%MPIRoot) THEN
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
    NumDens(1:nSpecies) = NumSpec(1:nSpecies) / GEO%MeshVolume
  ELSE
    NumDens(1:nSpecies) = NumSpec(1:nSpecies) * Species(1:nSpecies)%MacroParticleFactor / GEO%MeshVolume
  END IF

  IF(BGGas%BGGasSpecies.NE.0) NumDens(BGGas%BGGasSpecies) = BGGas%BGGasDensity

  IF(nSpecAnalyze.GT.1) NumDens(nSpecAnalyze) = SUM(NumDens(1:nSpecies))
END IF

END SUBROUTINE CalcNumberDensity


SUBROUTINE CalcTemperature(NumSpec,Temp,IntTemp,IntEn,TempTotal,Xi_Vib,Xi_Elec)
!===================================================================================================================================
! computes the temperature, subroutine performs all the MPI communication, do to it at ONE place and CORRECT
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PARTICLE_Vars         ,ONLY: nSpecies
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, PolyatomMolDSMC,CollisMode
USE MOD_DSMC_ElectronicModel  ,ONLY: CalcXiElec
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL, INTENT(IN)                   :: NumSpec(nSpecAnalyze)    ! number of real particles (already GLOBAL number)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)                  :: Temp(nSpecAnalyze)
REAL, INTENT(OUT)                  :: IntEn(nSpecAnalyze,3)
REAL, INTENT(OUT)                  :: IntTemp(nSpecies,3)
REAL, INTENT(OUT)                  :: TempTotal(nSpecAnalyze)
REAL, INTENT(OUT)                  :: Xi_Vib(nSpecies), Xi_Elec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iPolyatMole,iDOF,iSpec
!===================================================================================================================================


! next, calctranstemp

#if (PP_TimeDiscMethod!=1000)
CALL CalcTransTemp(NumSpec, Temp)
#else
CALL CalcTransTemp(Temp)
#endif

IF (CollisMode.GT.1) THEN
  CALL CalcIntTempsAndEn(NumSpec,IntTemp,IntEn)
  IF(PartMPI%MPIRoot)THEN
    TempTotal = 0.0
    Xi_Vib = 0.0
    Xi_Elec = 0.0
    DO iSpec = 1, nSpecies
      IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)).AND.(NumSpec(iSpec).GT.0)) THEN
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          IF(IntTemp(iSpec,1).GT.0) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
              Xi_Vib(iSpec) = Xi_Vib(iSpec) + 2*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/IntTemp(iSpec,1) &
                            /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/IntTemp(iSpec,1)) - 1)
            END DO
          ELSE
            Xi_Vib(iSpec) = 0.0
          END IF
        ELSE
          IF(IntTemp(iSpec,1).GT.0) THEN
            Xi_Vib(iSpec) = 2*SpecDSMC(iSpec)%CharaTVib/IntTemp(iSpec,1)/(EXP(SpecDSMC(iSpec)%CharaTVib/IntTemp(iSpec,1)) - 1)
          ELSE
            Xi_Vib(iSpec) = 0.0
          END IF
        END IF
        IF(IntTemp(iSpec,3).GT.0.0) Xi_Elec(iSpec) = CalcXiElec(IntTemp(iSpec,3), iSpec)
        TempTotal(iSpec) = (3*Temp(iSpec)+SpecDSMC(iSpec)%Xi_Rot*IntTemp(iSpec,2) &
                            + Xi_Vib(iSpec)*IntTemp(iSpec,1) + Xi_Elec(iSpec)*IntTemp(iSpec,3)) &
                            / (3+SpecDSMC(iSpec)%Xi_Rot+Xi_Vib(iSpec)+Xi_Elec(iSpec))
      ELSE
        IF(IntTemp(iSpec,3).GT.0.0) Xi_Elec(iSpec) = CalcXiElec(IntTemp(iSpec,3), iSpec)
        TempTotal(iSpec) = (3*Temp(iSpec) + Xi_Elec(iSpec)*IntTemp(iSpec,3)) / (3+Xi_Elec(iSpec))
      END IF
      IF(nSpecAnalyze.GT.1)THEN
        TempTotal(nSpecAnalyze) = TempTotal(nSpecAnalyze) + TempTotal(iSpec)*NumSpec(iSpec)
      END IF
    END DO
    IF(nSpecAnalyze.GT.1)THEN
      IF(NumSpec(iSpec).NE.0) THEN
        TempTotal(nSpecAnalyze) = TempTotal(nSpecAnalyze) / NumSpec(nSpecAnalyze)
      ELSE
        TempTotal(nSpecAnalyze)= 0.
      END IF
    END IF
  END IF
ELSE
  IF(PartMPI%MPIRoot)THEN
    TempTotal = Temp
    IntTemp   = 0.0
    IntEn     = 0.0
    Xi_Vib    = 0.0
    Xi_Elec   = 0.0
  END IF
END IF

END SUBROUTINE CalcTemperature



#if (PP_TimeDiscMethod!=1000)
SUBROUTINE CalcTransTemp(NumSpec, Temp)
#else
SUBROUTINE CalcTransTemp(Temp)
#endif
!===================================================================================================================================
! calculate the translational temperature of each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY : BoltzmannConst
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: nSpecies
#if (PP_TimeDiscMethod==1000)
USE MOD_LD_Vars               ,ONLY: BulkValues
#endif
#if (PP_TimeDiscMethod!=1000)
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PartState, Species, PDM
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if (PP_TimeDiscMethod!=1000)
REAL, INTENT(IN)                :: NumSpec(:)    !< global number of REAL particles in domain
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Temp(:)       !< output value is already the GLOBAL temperature
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec
REAL              :: TempDirec(nSpecies,3)
#if (PP_TimeDiscMethod!=1000)
REAL              :: PartVandV2(nSpecies, 6), Mean_PartV2(nSpecies, 3), MeanPartV_2(nSpecies,3)
INTEGER           :: i
#endif
!===================================================================================================================================

! Compute velocity averages
Temp = 0.0
! Sum up velocity
#if (PP_TimeDiscMethod==1000)
DO iSpec=1, nSpecies
  TempDirec(iSpec,1:3) = BulkValues(1)%BulkTemperature
  Temp(iSpec) = BulkValues(1)%BulkTemperature
END DO
#else
PartVandV2 = 0.
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    PartVandV2(PartSpecies(i),1:3) = PartVandV2(PartSpecies(i),1:3) + PartState(i,4:6) * GetParticleWeight(i)
    PartVandV2(PartSpecies(i),4:6) = PartVandV2(PartSpecies(i),4:6) + PartState(i,4:6)**2 * GetParticleWeight(i)
  END IF
END DO

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PartVandV2,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(PartVandV2  ,PartVandV2,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

IF(PartMPI%MPIRoot)THEN
  DO iSpec=1, nSpecies
    IF(NumSpec(iSpec).NE.0) THEN
      ! Compute velocity averages
      MeanPartV_2(iSpec,1:3)  = (PartVandV2(iSpec,1:3) / NumSpec(iSpec))**2       ! < |v| >**2
      Mean_PartV2(iSpec,1:3)  =  PartVandV2(iSpec,4:6) / NumSpec(iSpec)           ! < |v|**2 >
    ELSE
      MeanPartV_2(iSpec,1:3) = 0.
      Mean_PartV2(iSpec,1:3) = 0.
    END IF
    ! Compute temperatures
    TempDirec(iSpec,1:3) = Species(iSpec)%MassIC * (Mean_PartV2(iSpec,1:3) - MeanPartV_2(iSpec,1:3)) &
         / BoltzmannConst ! Trans Temp calculation is limited to one species
    Temp(iSpec) = (TempDirec(iSpec,1) + TempDirec(iSpec,2) + TempDirec(iSpec,3))/3
    IF(nSpecAnalyze.GT.1)THEN
      Temp(nSpecAnalyze) = Temp(nSpecAnalyze) + Temp(iSpec)*NumSpec(iSpec)
    END IF
  END DO
  IF(nSpecAnalyze.GT.1)THEN
    IF(NumSpec(iSpec).NE.0) THEN
      Temp(nSpecAnalyze)= Temp(nSpecAnalyze) / NumSpec(nSpecAnalyze)
    ELSE
      Temp(nSpecAnalyze)= 0.
    END IF
  END IF
END IF
#endif
END SUBROUTINE CalcTransTemp


SUBROUTINE CalcVelocities(PartVtrans, PartVtherm,NumSpec,SimNumSpec)
!===================================================================================================================================
! Calculates the drift and eigen velocity of all particles: PartVtotal = PartVtrans + PartVtherm
! PartVtrans(nSpecies,4) ! macroscopic velocity (drift velocity) A. Frohn: kinetische Gastheorie
! PartVtherm(nSpecies,4) ! microscopic velocity (eigen velocity)
!
! Note that the thermal velocity corresponds to the root mean square of the total velocity (in three dimensions), which is given by
!
!      v_th = SQRT(3 * kB * T / m)
!
! with kB : Boltzmann's constant
!      T  : temperature
!      m  : mass of the particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Analyze_Vars ,ONLY: VeloDirs
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, PDM, nSpecies, PartMPF, usevMPF
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: NumSpec(nSpecAnalyze)
INTEGER(KIND=8),INTENT(IN)     :: SimNumSpec(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: PartVtrans(nSpecies,4), PartVtherm(nSpecies,4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSpec
INTEGER                        :: i
INTEGER                        :: dir
#ifdef MPI
REAL                           :: RD(nSpecies*4)
#endif /*MPI*/
!===================================================================================================================================
! Compute velocity averages
  PartVtrans = 0.
  PartVtherm = 0.

  ! compute trans. velocity
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          IF (usevMPF) THEN
            PartVtrans(PartSpecies(i),dir) = PartVtrans(PartSpecies(i),dir) + PartState(i,dir+3) * PartMPF(i)
          ELSE
            PartVtrans(PartSpecies(i),dir) = PartVtrans(PartSpecies(i),dir) + PartState(i,dir+3)
          END IF
        END IF
      END DO
    END IF
  END DO

#ifdef MPI
  IF(PartMPI%MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,PartVtrans ,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  ELSE
    CALL MPI_REDUCE(PartVtrans  ,RD         ,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  END IF
#endif /*MPI*/

  IF(PartMPI%MPIRoot)THEN
    IF (usevMPF) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          DO iSpec = 1,nSpecies
            IF(NumSpec(iSpec).EQ.0)THEN
              PartVtrans(iSpec,dir) = 0.
            ELSE
              PartVtrans(iSpec,dir) = PartVtrans(iSpec,dir)/NumSpec(iSpec)
            END IF
          END DO ! iSpec = 1,nSpecies
        END IF
      END DO
    ELSE !no vMPF
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          DO iSpec = 1,nSpecies
            IF(SimNumSpec(iSpec).EQ.0)THEN
              PartVtrans(iSpec,dir) = 0.
            ELSE
              PartVtrans(iSpec,dir) = PartVtrans(iSpec,dir)/REAL(SimNumSpec(iSpec),8)
            END IF
          END DO ! iSpec = 1,nSpecies
        END IF
      END DO
    END IF !usevMPF
  END IF

#ifdef MPI
  CALL MPI_BCAST(PartVtrans,4*nSpecies, MPI_DOUBLE_PRECISION,0,PartMPI%COMM,iERROR)
#endif /*MPI*/

  ! calculate thermal velocity
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          IF (usevMPF) THEN
            PartVtherm(PartSpecies(i),dir) = PartVtherm(PartSpecies(i),dir) + PartMPF(i) * &
                (PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))*(PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))
          ELSE
            PartVtherm(PartSpecies(i),dir) = PartVtherm(PartSpecies(i),dir) + &
                (PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))*(PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))
          END IF
        END IF
      END DO
    END IF
  END DO

#ifdef MPI
  IF(PartMPI%MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,PartVtherm,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  ELSE
    CALL MPI_REDUCE(PartVtherm  ,RD        ,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  END IF
#endif /*MPI*/

  IF(PartMPI%MPIRoot)THEN
    DO dir = 1,3
      IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
        DO iSpec = 1,nSpecies
          IF (usevMPF) THEN
            IF(NumSpec(iSpec).EQ.0)THEN
              PartVtherm(iSpec,dir)=0.
            ELSE
              PartVtherm(iSpec,dir)=PartVtherm(iSpec,dir)/NumSpec(iSpec)
            END IF
          ELSE
            IF(SimNumSpec(iSpec).EQ.0)THEN
              PartVtherm(iSpec,dir)=0.
            ELSE
              PartVtherm(iSpec,dir)=PartVtherm(iSpec,dir)/REAL(SimNumSpec(iSpec),8)
            END IF
          END IF
        END DO ! iSpec = 1,nSpecies
      END IF
    END DO
 !   calc absolute value
    IF (VeloDirs(4)) THEN
      PartVtrans(:,4) = SQRT(PartVtrans(:,1)*PartVtrans(:,1) + PartVtrans(:,2)*PartVtrans(:,2) + PartVtrans(:,3)*PartVtrans(:,3))
      PartVtherm(:,4) = PartVtherm(:,1) + PartVtherm(:,2) + PartVtherm(:,3)
    END IF
    PartVtherm(:,:) = SQRT(PartVtherm(:,:))
  END IF
END SUBROUTINE CalcVelocities


SUBROUTINE CalcIntTempsAndEn(NumSpec,IntTemp,IntEn)
!===================================================================================================================================
! Calculation of internal Temps (TVib, TRot, Telec) and gives back the global values
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, PDM, nSpecies, usevMPF
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC
USE MOD_DSMC_Analyze          ,ONLY: CalcTVib, CalcTelec, CalcTVibPoly
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: NumSpec(nSpecAnalyze)    ! number of real particles (already GLOBAL number)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: IntTemp(nSpecies,3) , IntEn(nSpecAnalyze,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPart, iSpec
REAL                           :: EVib(nSpecies), ERot(nSpecies), Eelec(nSpecies), tempVib, NumSpecTemp
#ifdef MPI
REAL                           :: RD(nSpecies)
#endif /*MPI*/
!===================================================================================================================================
EVib    = 0.
ERot    = 0.
Eelec   = 0.
! set electronic state to zero
IntEn(:,:) = 0.
IntTemp(:,:) = 0.

! Sum up internal energies
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    iSpec = PartSpecies(iPart)
    EVib(iSpec) = EVib(iSpec) + PartStateIntEn(iPart,1) * GetParticleWeight(iPart)
    ERot(iSpec) = ERot(iSpec) + PartStateIntEn(iPart,2) * GetParticleWeight(iPart)
    IF (DSMC%ElectronicModel) THEN
      IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
        Eelec(iSpec) = Eelec(iSpec) + PartStateIntEn(iPart,3) * GetParticleWeight(iPart)
      END IF
    END IF
  END IF
END DO

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,EVib ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,ERot ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  IF(DSMC%ElectronicModel) CALL MPI_REDUCE(MPI_IN_PLACE,Eelec,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(EVib        ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  CALL MPI_REDUCE(ERot        ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  IF(DSMC%ElectronicModel) CALL MPI_REDUCE(Eelec       ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

! final computation is only done for the root
IF(PartMPI%MPIRoot)THEN
  ! Calc TVib, TRot
  DO iSpec = 1, nSpecies
    NumSpecTemp = NumSpec(iSpec)
    IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)).AND.(NumSpecTemp.GT.0.0)) THEN
      IF (SpecDSMC(iSpec)%PolyatomicMol.AND.(SpecDSMC(iSpec)%Xi_Rot.EQ.3)) THEN
        IntTemp(iSpec,2) = 2.0*ERot(iSpec)/(3.0*BoltzmannConst*NumSpecTemp)  !Calc TRot
      ELSE
        IntTemp(iSpec,2) = ERot(iSpec)/(BoltzmannConst*NumSpecTemp)  !Calc TRot
      END IF
      IF (EVib(iSpec)/NumSpecTemp.GT.SpecDSMC(iSpec)%EZeroPoint) THEN
        IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
          IntTemp(iSpec,1) = CalcTVibPoly(EVib(iSpec)/NumSpecTemp, iSpec)
        ELSE
          IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
            tempVib = (EVib(iSpec)/(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant)
            IF ((tempVib.GT.0.0) &
              .OR.(EVib(iSpec)/(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib).GT.DSMC%GammaQuant)) THEN
              IntTemp(iSpec,1) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(EVib(iSpec) &
                                /(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant))
            END IF
          ELSE                                            ! TSHO-model
            IntTemp(iSpec,1) = CalcTVib(SpecDSMC(iSpec)%CharaTVib, EVib(iSpec)/NumSpecTemp, SpecDSMC(iSpec)%MaxVibQuant)
          END IF
        END IF
      ELSE
        IntTemp(iSpec,1) = 0
      END IF
    ELSE
      IntTemp(iSpec,1) = 0
      IntTemp(iSpec,2) = 0
    END IF
    IF(DSMC%ElectronicModel) THEN
      IF(NumSpecTemp.GT.0) THEN
        IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
          IntTemp(iSpec,3) = CalcTelec(Eelec(iSpec)/NumSpecTemp,iSpec)
        END IF
      ELSE
        IntEn(iSpec,3) = 0.0
      END IF
    END IF
    IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
      ! MacroParticleFactor is included in the case of RadialWeighting (also in combination with variable time step)
      IntEn(iSpec,1) = EVib(iSpec)
      IntEn(iSpec,2) = ERot(iSpec)
      IF(DSMC%ElectronicModel) IntEn(iSpec,3) = Eelec(iSpec)
    ELSE
      IntEn(iSpec,1) = EVib(iSpec) * Species(iSpec)%MacroParticleFactor
      IntEn(iSpec,2) = ERot(iSpec) * Species(iSpec)%MacroParticleFactor
      IF(DSMC%ElectronicModel) IntEn(iSpec,3) = Eelec(iSpec) * Species(iSpec)%MacroParticleFactor
    END IF
  END DO
  ! Sums of the energy values
  IF(nSpecAnalyze.GT.1) THEN
    IntEn(nSpecAnalyze,1) = SUM(IntEn(:,1))
    IntEn(nSpecAnalyze,2) = SUM(IntEn(:,2))
    IF(DSMC%ElectronicModel) IntEn(nSpecAnalyze,3) = SUM(IntEn(:,3))
  END IF
END IF

END SUBROUTINE CalcIntTempsAndEn

#if (PP_TimeDiscMethod==42)
SUBROUTINE CollRates(CRate)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars     ,ONLY: CollInf, DSMC
USE MOD_TimeDisc_Vars ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: CRate(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iCase
!===================================================================================================================================

  DO iCase=1, CollInf%NumCase + 1
    CRate(iCase) =  DSMC%NumColl(iCase) / dt
  END DO
  DSMC%NumColl = 0
END SUBROUTINE CollRates

SUBROUTINE ReacRates(RRate, NumSpec,iter)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: ChemReac, DSMC
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_Particle_Vars         ,ONLY: Species, nSpecies
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
USE MOD_Particle_Analyze_Vars ,ONLY: PartAnalyzeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: RRate(:)
REAL,INTENT(IN)                 :: NumSpec(:) ! is the global number of real particles
INTEGER(KIND=8),INTENT(IN)      :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iReac
#ifdef MPI
REAL                            :: RD(1:ChemReac%NumOfReact)
#endif /*MPI*/
!===================================================================================================================================

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE    ,ChemReac%NumReac,ChemReac%NumOfReact,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(ChemReac%NumReac,RD              ,ChemReac%NumOfReact,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(PartMPI%MPIRoot)THEN
  DO iReac=1, ChemReac%NumOfReact
    IF ((NumSpec(ChemReac%DefinedReact(iReac,1,1)).GT.0).AND.(NumSpec(ChemReac%DefinedReact(iReac,1,2)).GT.0)) THEN
      SELECT CASE(TRIM(ChemReac%ReactType(iReac)))
      CASE('R','r','rQK')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%MeshVolume**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(nSpecies+1))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) * GEO%MeshVolume**2 &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor / (dt * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))     &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2))    &
               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(nSpecies+1))
        END IF
      CASE('D','E','i','iQK','x')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor* GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      END SELECT
    END IF
  END DO
END IF
ChemReac%NumReac = 0.
ChemReac%ReacCount = 0
ChemReac%ReacCollMean = 0.0
ChemReac%ReacCollMeanCount = 0
! Consider Part-AnalyzeStep
IF(PartAnalyzeStep.GT.1)THEN
  IF(PartAnalyzeStep.EQ.HUGE(PartAnalyzeStep))THEN
    DO iReac=1, ChemReac%NumOfReact
      RRate(iReac) = RRate(iReac) / iter
    END DO ! iReac=1, ChemReac%NumOfReact
  ELSE
    DO iReac=1, ChemReac%NumOfReact
      RRate(iReac) = RRate(iReac) / MIN(PartAnalyzeStep,iter)
    END DO ! iReac=1, ChemReac%NumOfReact
  END IF
END IF


END SUBROUTINE ReacRates
#endif

#if ( PP_TimeDiscMethod == 42)
#ifdef CODE_ANALYZE
SUBROUTINE ElectronicTransition (  Time, NumSpec )
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars          ,ONLY: DSMC, SpecDSMC
USE MOD_TimeDisc_Vars      ,ONLY: dt
USE MOD_Particle_Vars      ,ONLY: nSpecies, Species
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)                :: Time
REAL,INTENT(IN)               :: NumSpec(:)
INTEGER                        :: iSpec, iSpec2, iQua1, iQua2, MaxElecQua
! accary of kf
!===================================================================================================================================

IF ( DSMC%ElectronicModel ) THEN
! kf = d n_of_N^i / dt / ( n_of_N^i n_of_M) )
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      DO iSpec2 = 1, nSpecies
   ! calculaction of kf for each reaction
!      MaxElecQua = SpecDSMC(iSpec)%MaxElecQuant
      MaxElecQua = 2
      ! for first tests only consider the first 10 transition levels
        DO iQua1 = 0, MaxElecQua
          DO iQua2 = 0, MaxElecQua
        ! calculate kf
        ! kf = ( d n_of_N^i / d t )  / ( n_of_N^i n_of_M )
            IF ( (NumSpec(iSpec2) .ne. 0) .and. (SpecDSMC(iSpec)%levelcounter(iQua1) .ne. 0 ) ) THEN
              SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2) = &
                                    SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2) * GEO%MeshVolume &
                                  / ( dt * SpecDSMC(iSpec)%levelcounter(iQua1) * NumSpec(iSpec2) *           &
                                      Species(iSpec2)%MacroParticleFactor )
            END IF
!             print*,SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2)
          END DO
        END DO
      END DO
    END IF
  END DO
  CALL WriteEletronicTransition( Time )
  ! nullyfy
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      SpecDSMC(iSpec)%ElectronicTransition = 0
    END IF
  END DO
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ElectronicTransition
#endif
#endif

#if ( PP_TimeDiscMethod == 42)
#ifdef CODE_ANALYZE
SUBROUTINE WriteEletronicTransition ( Time )
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: Getfreeunit
USE MOD_DSMC_Vars     ,ONLY: DSMC, SpecDSMC
USE MOD_Particle_Vars ,ONLY: nSpecies
!   ! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)                :: Time
INTEGER                        :: iSpec, iSpec2, iunit, iQua1, iQua2, MaxElecQua, ii
CHARACTER(LEN=128)             :: FileNameTransition
LOGICAL                        :: bExist
! accary of kf
!-----------------------------------------------------------------------------------------------------------------------------------
bExist = .false.
IF ( DSMC%ElectronicModel ) THEN
! kf = d n_of_N^i / dt / ( n_of_N^i n_of_M) )
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
!        MaxElecQua = SpecDSMC(iSpec)%MaxElecQuant
      MaxElecQua = 2
      ! output to transition file
      FileNameTransition = trim(SpecDSMC(iSpec)%Name)//'_Transition.csv'
      INQUIRE( FILE = FileNameTransition, EXIST=bExist)
!-----------------------------------------------------------------------------------------------------------------------------------
      IF ( bExist .EQV. .false. ) THEN
        iunit=GETFREEUNIT()
        OPEN(UNIT=iunit,FILE=FileNameTransition,FORM='FORMATTED',STATUS='UNKNOWN')
!         ! writing header
        WRITE(iunit,'(A6,A5)',ADVANCE='NO') '001-TIME', ' '
        ii = 2
        DO iSpec2 = 1, nSpecies
          DO iQua1 = 0, MaxElecQua
            DO iQua2 = 0, MaxElecQua
              WRITE(iunit,'(I3.3,A,I2.2,A,I2.2,A,I2.2,A)', ADVANCE='NO') ii,'_Species',iSpec2,'_',iQua1,'_to_',iQua2,'  '
              WRITE(iunit,'(A1)',ADVANCE='NO') ','
              ii = ii + 1
            END DO
          END DO
        END DO
      ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
!         ! writing header
        iunit=GETFREEUNIT()
        OPEN(unit=iunit,FILE=FileNameTransition,FORM='Formatted',POSITION='APPEND',STATUS='old')
        WRITE(iunit,104,ADVANCE='NO') TIME
        DO iSpec2 = 1, nSpecies
!       ! calculaction of kf for each reaction
          DO iQua1 = 0, MaxElecQua
            DO iQua2 = 0, MaxElecQua
            ! write values to file
              WRITE(iunit,'(A1)',ADVANCE='NO') ','
              WRITE(iunit,104,ADVANCE='NO') SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2)
            END DO
          END DO
        END DO
        WRITE(iunit,'(A)') ' '
      END IF
      CLOSE(unit = iunit)
    END IF
  END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
104    FORMAT (e25.14)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE WriteEletronicTransition
#endif
#endif

!  SUBROUTINE TrackingParticlePosition(time)
!  !===================================================================================================================================
!  ! Initializes variables necessary for analyse subroutines
!  !===================================================================================================================================
!  ! MODULES
!  USE MOD_Globals
!  USE MOD_Preproc
!  USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PEM
!  USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
!  USE MOD_Particle_Analyze_Vars ,ONLY: printDiff,printDiffVec,printDiffTime
!  #if defined(LSERK) || defined(IMPA) || defined(ROS)
!  USE MOD_Equation_Vars         ,ONLY: c2_inv
!  #endif
!  ! IMPLICIT VARIABLE HANDLING
!  IMPLICIT NONE
!  !-----------------------------------------------------------------------------------------------------------------------------------
!  ! INPUT VARIABLES
!  REAL,INTENT(IN)                :: time
!  !-----------------------------------------------------------------------------------------------------------------------------------
!  ! OUTPUT VARIABLES
!  !-----------------------------------------------------------------------------------------------------------------------------------
!  ! LOCAL VARIABLES
!  INTEGER            :: i,iunit,iPartState
!  CHARACTER(LEN=60) :: TrackingFilename!,hilf
!  LOGICAL            :: fexist
!  REAL               :: diffPos,diffVelo
!  !===================================================================================================================================
!
!  !WRITE(UNIT=hilf,FMT='(I6.6)') MyRank
!  !TrackingFilename = ('MyRank'//TRIM(hilf)//'_ParticlePosition.csv')
!  TrackingFilename = ('ParticlePosition.csv')
!
!  INQUIRE(FILE = TrackingFilename, EXIST=fexist)
!  IF(.NOT.fexist) THEN
!   IF(PartMPI%MPIRoot)THEN
!     OPEN(NEWUNIT=iunit,FILE=TrackingFilename,FORM='FORMATTED',STATUS='UNKNOWN')
!     !CALL FLUSH (iunit)
!      ! writing header
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') '001-TIME', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartNum', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosX', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosY', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosZ', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelX', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelY', ' '
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelZ', ' '
!  #if defined(LSERK) || defined(IMPA) || defined(ROS)
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'gamma', ' '
!  #endif
!      WRITE(iunit,'(A1)',ADVANCE='NO') ','
!      WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'Element '
!      CLOSE(iunit)
!    END IF
!  ELSE
!    OPEN(NEWUNIT=iunit,FILE=TrackingFileName,FORM='Formatted',POSITION='APPEND',STATUS='old')
!    !CALL FLUSH (iunit)
!    DO i=1,PDM%ParticleVecLength
!      IF (PDM%ParticleInside(i)) THEN
!        WRITE(iunit,104,ADVANCE='NO') TIME
!        WRITE(iunit,'(A1)',ADVANCE='NO') ','
!        WRITE(iunit,'(I12)',ADVANCE='NO') i
!        DO iPartState=1,6
!          WRITE(iunit,'(A1)',ADVANCE='NO') ','
!          WRITE(iunit,104,ADVANCE='NO') PartState(i,iPartState)
!        END DO
!  #if defined(LSERK) || defined(IMPA) || defined(ROS)
!          WRITE(iunit,'(A1)',ADVANCE='NO') ','
!          WRITE(iunit,104,ADVANCE='NO') 1./SQRT(1-(DOT_PRODUCT(PartState(i,4:6),PartState(i,4:6))*c2_inv))
!  #endif
!        WRITE(iunit,'(A1)',ADVANCE='NO') ','
!        WRITE(iunit,'(I12)',ADVANCE='NO') PEM%Element(i)
!        WRITE(iunit,'(A)') ' '
!       END IF
!    END DO
!    CLOSE(iunit)
!  END IF
!  IF (printDiff) THEN
!    diffPos=0.
!    diffVelo=0.
!    IF (time.GE.printDiffTime) THEN
!      printDiff=.FALSE.
!      DO iPartState=1,3
!        diffPos=diffPos+(printDiffVec(iPartState)-PartState(1,iPartState))**2
!        diffVelo=diffVelo+(printDiffVec(iPartState+3)-PartState(1,iPartState+3))**2
!      END DO
!      WRITE(*,'(A,e24.14,x,e24.14)') 'L2-norm from printDiffVec: ',SQRT(diffPos),SQRT(diffVelo)
!    END IF
!  END IF
!  104    FORMAT (e25.14)
!
!  END SUBROUTINE TrackingParticlePosition

!----------------------------------------------------------------------------------------------------------------------------------!
!> Write particle info to ParticlePosition.csv file
!> time, pos, velocity, gamma, element
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteParticleTrackingData(time,iter)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals          ,ONLY: MPIRoot,FILEEXISTS,unit_stdout
USE MOD_Restart_Vars     ,ONLY: DoRestart
USE MOD_Globals          ,ONLY: abort

USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PEM
USE MOD_Particle_Analyze_Vars ,ONLY: printDiff,printDiffVec,printDiffTime
#if defined(LSERK) || defined(IMPA) || defined(ROS)
USE MOD_Equation_Vars         ,ONLY: c2_inv
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                  :: time
INTEGER(KIND=8),INTENT(IN)       :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=20),PARAMETER              :: outfile='ParticlePosition.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=150)                       :: formatStr
#if defined(LSERK) || defined(IMPA) || defined(ROS)
INTEGER,PARAMETER                        :: nOutputVar=10
#else
INTEGER,PARAMETER                        :: nOutputVar=9
#endif
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    '001-time',     &
    'PartNum',  &
    'PartPosX', &
    'PartPosY', &
    'PartPosZ', &
    'PartVelX', &
    'PartVelY', &
    'PartVelZ', &
#if defined(LSERK) || defined(IMPA) || defined(ROS)
    'gamma',    &
#endif
    'Element'/)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2
CHARACTER(LEN=1),PARAMETER               :: delimiter=","
LOGICAL                                  :: FileExist,CreateFile
REAL                                     :: diffPos,diffVelo
INTEGER                                  :: iPartState
!===================================================================================================================================
! only the root shall write this file
IF(.NOT.MPIRoot)RETURN

! check if file is to be created
CreateFile=.TRUE.
IF(iter.GT.0)CreateFile=.FALSE.                             ! don't create new file if this is not the first iteration
IF((DoRestart).AND.(FILEEXISTS(outfile)))CreateFile=.FALSE. ! don't create new file if this is a restart and the file already exists
!                                                           ! assume continued simulation and old load balance data is still needed

! check if new file with header is to be created
INQUIRE(FILE = outfile, EXIST=FileExist)
IF(.NOT.FileExist)CreateFile=.TRUE.                         ! if no file exists, create one

! create file with header
IF(CreateFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=""
  DO I=1,nOutputVar
    WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO I=1,nOutputVar
    IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF

! Print info to file
IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14)')'(',nOutputVar,CSVFORMAT
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      WRITE(tmpStr2,formatStr)&
          " ",time, &                                                                     ! time
          delimiter,REAL(i), &                                                            ! PartNum
          delimiter,PartState(i,1), &                                                     ! PartPosX
          delimiter,PartState(i,2), &                                                     ! PartPosY
          delimiter,PartState(i,3), &                                                     ! PartPosZ
          delimiter,PartState(i,4), &                                                     ! PartVelX
          delimiter,PartState(i,5), &                                                     ! PartVelY
          delimiter,PartState(i,6), &                                                     ! PartVelZ
#if defined(LSERK) || defined(IMPA) || defined(ROS)
          delimiter,1./SQRT(1-(DOT_PRODUCT(PartState(i,4:6),PartState(i,4:6))*c2_inv)), & ! gamma
#endif
          delimiter,REAL(PEM%Element(i))                                                  ! Element
      WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
    END IF
  END DO
  CLOSE(ioUnit)
ELSE
  SWRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write particle tracking info!"
END IF

! printDiff
IF (printDiff) THEN
  diffPos=0.
  diffVelo=0.
  IF (time.GE.printDiffTime) THEN
    printDiff=.FALSE.
    DO iPartState=1,3
      diffPos=diffPos+(printDiffVec(iPartState)-PartState(1,iPartState))**2
      diffVelo=diffVelo+(printDiffVec(iPartState+3)-PartState(1,iPartState+3))**2
    END DO
    WRITE(*,'(A,e24.14,x,e24.14)') 'L2-norm from printDiffVec: ',SQRT(diffPos),SQRT(diffVelo)
  END IF
END IF
END SUBROUTINE WriteParticleTrackingData


#ifdef CODE_ANALYZE
!----------------------------------------------------------------------------------------------------------------------------------!
!> Write analytic particle info to ParticlePositionAnalytic.csv file
!> time, pos, velocity
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteParticleTrackingDataAnalytic(time,iter,PartStateAnalytic)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals               ,ONLY: MPIRoot,FILEEXISTS,unit_stdout
USE MOD_Restart_Vars          ,ONLY: DoRestart
USE MOD_Globals               ,ONLY: abort
USE MOD_PICInterpolation_Vars ,ONLY: L_2_Error_Part
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                  :: time
INTEGER(KIND=8),INTENT(IN)       :: iter
REAL(KIND=8),INTENT(IN)          :: PartStateAnalytic(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=28),PARAMETER              :: outfile='ParticlePositionAnalytic.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=150)                       :: formatStr
INTEGER,PARAMETER                        :: nOutputVar=13
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    '001-time',     &
    'PartPosX_Analytic', &
    'PartPosY_Analytic', &
    'PartPosZ_Analytic', &
    'PartVelX_Analytic', &
    'PartVelY_Analytic', &
    'PartVelZ_Analytic', &
    'L2_PartPosX'      , &
    'L2_PartPosY'      , &
    'L2_PartPosZ'      , &
    'L2_PartVelX'      , &
    'L2_PartVelY'      , &
    'L2_PartVelZ'        &
    /)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2
CHARACTER(LEN=1),PARAMETER               :: delimiter=","
LOGICAL                                  :: FileExist,CreateFile
!===================================================================================================================================
! only the root shall write this file
IF(.NOT.MPIRoot)RETURN

! check if file is to be created
CreateFile=.TRUE.
IF(iter.GT.0)CreateFile=.FALSE.                             ! don't create new file if this is not the first iteration
IF((DoRestart).AND.(FILEEXISTS(outfile)))CreateFile=.FALSE. ! don't create new file if this is a restart and the file already exists
!                                                           ! assume continued simulation and old load balance data is still needed

! check if new file with header is to be created
INQUIRE(FILE = outfile, EXIST=FileExist)
IF(.NOT.FileExist)CreateFile=.TRUE.                         ! if no file exists, create one

! create file with header
IF(CreateFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=""
  DO I=1,nOutputVar
    WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO I=1,nOutputVar
    IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible relimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF

! Print info to file
IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14)')'(',nOutputVar,CSVFORMAT
  WRITE(tmpStr2,formatStr)&
      " ",time, &                           ! time
      delimiter,PartStateAnalytic(1), &     ! PartPosX analytic solution
      delimiter,PartStateAnalytic(2), &     ! PartPosY analytic solution
      delimiter,PartStateAnalytic(3), &     ! PartPosZ analytic solution
      delimiter,PartStateAnalytic(4), &     ! PartVelX analytic solution
      delimiter,PartStateAnalytic(5), &     ! PartVelY analytic solution
      delimiter,PartStateAnalytic(6), &     ! PartVelZ analytic solution
      delimiter,L_2_Error_Part(1), &     ! L2 error for PartPosX solution
      delimiter,L_2_Error_Part(2), &     ! L2 error for PartPosY solution
      delimiter,L_2_Error_Part(3), &     ! L2 error for PartPosZ solution
      delimiter,L_2_Error_Part(4), &     ! L2 error for PartVelX solution
      delimiter,L_2_Error_Part(5), &     ! L2 error for PartVelY solution
      delimiter,L_2_Error_Part(6)        ! L2 error for PartVelZ solution
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  SWRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write particle tracking (analytic) info!"
END IF

END SUBROUTINE WriteParticleTrackingDataAnalytic
#endif /* CODE_ANALYZE */


PURE Function CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars ,ONLY: c2, c2_inv
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species
USE MOD_PARTICLE_Vars ,ONLY: usevMPF
USE MOD_Particle_Vars ,ONLY: PartLorentzType
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: iPart
REAL                               :: CalcEkinPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: partV2, gamma1, WeightingFactor
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  WeightingFactor = GetParticleWeight(iPart)
ELSE
  WeightingFactor = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
END IF

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv)
  CalcEkinPart=(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2 * WeightingFactor
ELSE
  partV2 = PartState(iPart,4) * PartState(iPart,4) &
         + PartState(iPart,5) * PartState(iPart,5) &
         + PartState(iPart,6) * PartState(iPart,6)
  IF (partV2.LT.1e6)THEN
    CalcEkinPart= 0.5 * Species(PartSpecies(iPart))%MassIC * partV2 * WeightingFactor
  ELSE
    gamma1=partV2*c2_inv
    gamma1=1.0/SQRT(1.-gamma1)
    CalcEkinPart=(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2 * WeightingFactor
  END IF ! ipartV2
END IF
END FUNCTION CalcEkinPart


SUBROUTINE CalcPowerDensity()
!===================================================================================================================================
! Used to average the source terms per species
!   * compute the power density of the considered species
!     scalar product of < j, E >, with j- current density and E - electric field
!   * compute the charge density of the considered species
!   * compute the current density of the considered species
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Timeaverage_Vars  ,ONLY: DoPowerDensity,PowerDensity
USE MOD_Particle_Vars     ,ONLY: nSpecies,PartSpecies,PDM
USE MOD_PICDepo_Vars      ,ONLY: PartSource
USE MOD_Part_RHS          ,ONLY: PartVeloToImp
USE MOD_Preproc
USE MOD_PICDepo           ,ONLY: Deposition
#ifndef PP_HDG
USE MOD_DG_Vars           ,ONLY: U
#else
#if PP_nVar==1
USE MOD_Equation_Vars     ,ONLY: E
#else
#endif
#endif
#ifdef MPI
USE MOD_Particle_MPI      ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars ,ONLY: PartMPIExchange
USE MOD_Particle_MPI_Vars ,ONLY: DoExternalParts
USE MOD_Particle_MPI_Vars ,ONLY: ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iSpec,iSpec2
INTEGER              :: iElem,i,j,k,iPart
LOGICAL              :: doParticle(1:PDM%MaxParticleNumber)
!===================================================================================================================================

iSpec2=0
PowerDensity=0.
DO iSpec=1,nSpecies
  IF(.NOT.DoPowerDensity(iSpec)) CYCLE
  iSpec2=iSpec2+1
  ! mark particle
  DoParticle(:)=.FALSE.
  DO iPart=1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(iPart))THEN
      IF(PartSpecies(iPart).EQ.iSpec)THEN
        DoParticle(iPart)=.TRUE.
      END IF
    END IF ! ParticleInside
  END DO ! iPart


  ! map particle from gamma v to v
  CALL PartVeloToImp(VeloToImp=.FALSE.,doParticle_In=DoParticle(1:PDM%ParticleVecLength))

  ! communicate shape function particles
#ifdef MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles(doParticle_In=DoParticle)
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
    ! set exchanged number of particles to zero
    PartMPIExchange%nMPIParticles=0
  END IF
#endif /*MPI*/

  ! compute source terms
  ! compute particle source terms on field solver of considered species
  CALL Deposition(doInnerParts=.TRUE.,doParticle_In=DoParticle(1:PDM%ParticleVecLength))
  CALL Deposition(doInnerParts=.FALSE.,doParticle_In=DoParticle(1:PDM%ParticleVecLength))
  ! map particle from v to gamma v
  CALL PartVeloToImp(VeloToImp=.TRUE.,doParticle_In=DoParticle(1:PDM%ParticleVecLength))

  ! compute power density
  DO iElem=1,PP_nElems
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! 1:3 PowerDensity, 4 charge density
#ifndef PP_HDG
          PowerDensity(1,i,j,k,iElem,iSpec2)=PartSource(1,i,j,k,iElem)*U(1,i,j,k,iElem)
          PowerDensity(2,i,j,k,iElem,iSpec2)=PartSource(2,i,j,k,iElem)*U(2,i,j,k,iElem)
          PowerDensity(3,i,j,k,iElem,iSpec2)=PartSource(3,i,j,k,iElem)*U(3,i,j,k,iElem)
          PowerDensity(4,i,j,k,iElem,iSpec2)=PartSource(4,i,j,k,iElem)
#else
#if PP_nVar==1
          PowerDensity(1,i,j,k,iElem,iSpec2)=PartSource(1,i,j,k,iElem)*E(1,i,j,k,iElem)
          PowerDensity(2,i,j,k,iElem,iSpec2)=PartSource(2,i,j,k,iElem)*E(2,i,j,k,iElem)
          PowerDensity(3,i,j,k,iElem,iSpec2)=PartSource(3,i,j,k,iElem)*E(3,i,j,k,iElem)
#else
          PowerDensity(1:3,i,j,k,iElem,iSpec2)=0.
#endif
          PowerDensity(4,i,j,k,iElem,iSpec2)=PartSource(4,i,j,k,iElem)
#endif
          ! 5:7 current density
          PowerDensity(5,i,j,k,iElem,iSpec2)=PartSource(1,i,j,k,iElem)
          PowerDensity(6,i,j,k,iElem,iSpec2)=PartSource(2,i,j,k,iElem)
          PowerDensity(7,i,j,k,iElem,iSpec2)=PartSource(3,i,j,k,iElem)
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
  END DO ! iElem=1,PP_nElems
END DO

END SUBROUTINE CalcPowerDensity


PURE FUNCTION PARTISELECTRON(PartID)
!===================================================================================================================================
! check if particle is an electron (species-charge = -1.609)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: ElementaryCharge
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL            :: PartIsElectron  !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SpeciesID
!===================================================================================================================================

PartIsElectron=.FALSE.
SpeciesID = PartSpecies(PartID)
IF(Species(SpeciesID)%ChargeIC.GT.0.0) RETURN
IF(NINT(Species(SpeciesID)%ChargeIC/(-ElementaryCharge)).EQ.1) PartIsElectron=.TRUE.

END FUNCTION PARTISELECTRON


SUBROUTINE CalculateElectronIonDensityCell()
!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY:ElementaryCharge
USE MOD_Particle_Mesh_Vars     ,ONLY:GEO,NbrOfRegions
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell,IonDensityCell,NeutralDensityCell,ChargeNumberCell
USE MOD_Particle_Vars          ,ONLY:Species,PartSpecies,PDM,PEM,usevMPF
USE MOD_Preproc                ,ONLY:PP_nElems
USE MOD_PIC_Analyze            ,ONLY:CalculateBRElectronsPerCell
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_part_tools             ,ONLY: GetParticleWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPart,iElem,RegionID
REAL                 :: charge, MPF
!===================================================================================================================================

! nullify
ElectronDensityCell=0.
     IonDensityCell=0.
 NeutralDensityCell=0.
   ChargeNumberCell=0.

! loop over all particles and count the number of electrons per cell
! CAUTION: we need the number of all real particle instead of simulated particles
DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
    MPF = GetParticleWeight(iPart)
  ELSE
    MPF = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF
  ASSOCIATE ( &
    ElemID  => PEM%Element(iPart)                              )  ! Element ID
    ASSOCIATE ( &
      n_e    => ElectronDensityCell(ElemID),& ! Electron density (cell average)
      n_i    => IonDensityCell(ElemID)     ,& ! Ion density (cell average)
      n_n    => NeutralDensityCell(ElemID) ,& ! Neutral density (cell average)
      Z      => ChargeNumberCell(ElemID)   )  ! Charge number (cell average)
      charge = Species(PartSpecies(iPart))%ChargeIC/ElementaryCharge
      IF(PARTISELECTRON(iPart))THEN ! electrons
        n_e = n_e + MPF
      ELSEIF(ABS(charge).GT.0.0)THEN ! ions (positive or negative)
        n_i = n_i + MPF
        Z   = Z   + charge*MPF
      ELSE ! neutrals
        n_n  = n_n + MPF
      END IF
    END ASSOCIATE
  END ASSOCIATE
END DO ! iPart
IF (NbrOfRegions .GT. 0) THEN !check for BR electrons
  DO iElem=1,PP_nElems
    RegionID=GEO%ElemToRegion(iElem)
    IF (RegionID.GT.0) THEN
      IF (ElectronDensityCell(iElem).NE.0.) CALL abort(&
__STAMP__&
,'Mixed BR and kinetic electrons are not implemented in CalculateElectronIonDensityCell yet!')
      CALL CalculateBRElectronsPerCell(iElem,RegionID,ElectronDensityCell(iElem))
    END IF
  END DO ! iElem=1,PP_nElems
END IF

! loop over all elements and divide by volume
DO iElem=1,PP_nElems
  ElectronDensityCell(iElem)=ElectronDensityCell(iElem)/GEO%Volume(iElem)
       IonDensityCell(iElem)=IonDensityCell(iElem)     /GEO%Volume(iElem)
   NeutralDensityCell(iElem)=NeutralDensityCell(iElem) /GEO%Volume(iElem)
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateElectronIonDensityCell


SUBROUTINE CalculateElectronTemperatureCell()
!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst,ElectronMass,ElementaryCharge
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO,NbrOfRegions
USE MOD_Preproc               ,ONLY: PP_nElems
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronTemperatureCell
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,usevMPF,Species,PartSpecies,PartState,RegionElectronRef
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
USE MOD_part_tools            ,ONLY: GetParticleWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart,iElem,ElemID,Method,RegionID
REAL    :: nElectronsPerCell(1:PP_nElems)
REAL    ::  PartVandV2(1:PP_nElems,1:6)
REAL    :: Mean_PartV2(1:3)
REAL    :: MeanPartV_2(1:3)
REAL    ::   TempDirec(1:3)
REAL    :: WeightingFactor
!===================================================================================================================================
IF (NbrOfRegions .GT. 0) THEN ! check for BR electrons
  DO iElem=1,PP_nElems
    RegionID=GEO%ElemToRegion(iElem)
    IF (RegionID.GT.0) THEN
      ElectronTemperatureCell(iElem) = RegionElectronRef(3,RegionID)*ElementaryCharge/BoltzmannConst ! convert eV to K
    END IF
  END DO ! iElem=1,PP_nElems
  RETURN ! Mixed BR and kinetic electrons are not implemented yet!
END IF

! nullify
ElectronTemperatureCell=0.
nElectronsPerCell      =0.

! hard-coded
Method=1 ! 0: <E> = (3/2)*<k_B*T> (keeps drift velocity, hence, over-estimates the temperature when drift becomes important)
!        ! 1: remove drift from temperature calculation

PartVandV2 = 0.
! 1.   loop over all particles and sum-up the electron energy per cell and count the number of electrons per cell
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    IF(.NOT.PARTISELECTRON(iPart)) CYCLE  ! ignore anything that is not an electron
    ElemID                      = PEM%Element(iPart)
    IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
      WeightingFactor = GetParticleWeight(iPart)
    ELSE
      WeightingFactor = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
    END IF
    nElectronsPerCell(ElemID) = nElectronsPerCell(ElemID) + WeightingFactor
    ! Determine velocity or kinetic energy
    SELECT CASE(Method)
    CASE(0) ! 1.0   for distributions where the drift is negligible
      ElectronTemperatureCell(ElemID) = ElectronTemperatureCell(ElemID)+CalcEkinPart(iPart)
    CASE(1) ! 1.1   remove drift from distribution
      PartVandV2(ElemID,1:3) = PartVandV2(ElemID,1:3) + PartState(iPart,4:6)    * WeightingFactor
      PartVandV2(ElemID,4:6) = PartVandV2(ElemID,4:6) + PartState(iPart,4:6)**2 * WeightingFactor
    END SELECT
  END IF ! ParticleInside
END DO ! iPart

! 2.   loop over all elements and divide by electrons per cell to get average kinetic energy
SELECT CASE(Method)
CASE(0) ! 2.0   for distributions where the drift is negligible
  DO iElem=1,PP_nElems
    IF(nElectronsPerCell(iElem).GT.0.) THEN
      ! <E> = (3/2)*<k_B*T>
      ElectronTemperatureCell(iElem)  = 2.*ElectronTemperatureCell(iElem)/(3.*nElectronsPerCell(iElem)*BoltzmannConst)
    END IF
  END DO ! iElem=1,PP_nElems
CASE(1) ! 2.1   remove drift from distribution
  DO iElem=1,PP_nElems
    IF(nElectronsPerCell(iElem).LT.2.) THEN ! only calculate the temperature when more than one electron are present
      ElectronTemperatureCell(iElem) = 0.0
    ELSE
      ! Compute velocity averages
      MeanPartV_2(1:3)  = (PartVandV2(iElem,1:3) / nElectronsPerCell(iElem))**2 ! < |v| >**2
      Mean_PartV2(1:3)  =  PartVandV2(iElem,4:6) / nElectronsPerCell(iElem)     ! < |v|**2 >
      ! Compute temperatures
      TempDirec(1:3) = ElectronMass * (Mean_PartV2(1:3) - MeanPartV_2(1:3)) / BoltzmannConst
      ElectronTemperatureCell(iElem) = (TempDirec(1) + TempDirec(2) + TempDirec(3))/3.0
      IF(ElectronTemperatureCell(iElem).LT.0.0)THEN
        ElectronTemperatureCell(iElem)=0.0
      END IF
    END IF
  END DO
END SELECT

END SUBROUTINE CalculateElectronTemperatureCell


SUBROUTINE CalculatePlasmaFrequencyCell()
!===================================================================================================================================
! use the number of electron density to compute the plasma frequency per cell using fixed and global values for the
! electron charge, electronmass and eps0
! CAUTION: if c!=3e8 m/s the computed frequency may be wrong
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc                ,ONLY:PP_nElems
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell,PlasmaFrequencyCell
USE MOD_Globals_Vars           ,ONLY:ElementaryCharge,ElectronMass
USE MOD_Equation_Vars          ,ONLY:Eps0
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
PlasmaFrequencyCell=0.

! loop over all elements and compute the plasma frequency with the use of the electron density
DO iElem=1,PP_nElems
  PlasmaFrequencyCell(iElem) = SQRT((ElectronDensityCell(iElem)*ElementaryCharge*ElementaryCharge)/(ElectronMass*Eps0))
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePlasmaFrequencyCell

SUBROUTINE CalculatePICTimeStepCell()
!===================================================================================================================================
! use the plasma frequency per cell to estimate the pic time step
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc                ,ONLY:PP_nElems
USE MOD_Particle_Analyze_Vars  ,ONLY:PlasmaFrequencyCell,PICTimeStepCell
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
PICTimeStepCell=0

! loop over all elements and compute the PIC-timestep with the plasma frequency
DO iElem=1,PP_nElems
  IF(PlasmaFrequencyCell(iElem).LE.0) CYCLE
  PICTimeStepCell(iElem) = 0.2 / PlasmaFrequencyCell(iElem)
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePICTimeStepCell


SUBROUTINE CalculateDebyeLengthCell()
!===================================================================================================================================
! use the number of electron density and electron temperature to compute the cold Debye-length per cell
! CAUTION: use SI-units
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc                ,ONLY:PP_nElems
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell,ElectronTemperatureCell,DebyeLengthCell,QuasiNeutralityCell
USE MOD_Globals_Vars           ,ONLY:ElementaryCharge, BoltzmannConst
USE MOD_Equation_Vars          ,ONLY:Eps0
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
DebyeLengthCell=0.

! loop over all elements and compute the plasma frequency with the use of the electron density
DO iElem=1,PP_nElems
  IF(ElectronDensityCell(iElem).LE.0.0) CYCLE ! ignore cells in which no electrons are present
  IF(QuasiNeutralityCell(iElem).LE.0.0) CYCLE ! ignore cells in which quasi neutrality is not possible
  DebyeLengthCell(iElem) = SQRT( (eps0*BoltzmannConst*ElectronTemperatureCell(iElem))/&
                                 (ElectronDensityCell(iElem)*(ElementaryCharge**2))       )
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateDebyeLengthCell


SUBROUTINE CalculatePPDCell()
!===================================================================================================================================
! Calculate the points per Debye length for each cell
! PointsPerDebyeLength: PPD = (p+1)*lambda_D/L_cell
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc                ,ONLY:PP_nElems,PP_N
USE MOD_Particle_Analyze_Vars  ,ONLY:DebyeLengthCell,PPDCell
USE MOD_Particle_Mesh_Vars     ,ONLY:GEO
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! loop over all elements
DO iElem=1,PP_nElems
  PPDCell(iElem) = (REAL(PP_N)+1.0)*DebyeLengthCell(iElem)/GEO%CharLength(iElem)
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePPDCell


SUBROUTINE CalculateIonizationCell()
!===================================================================================================================================
! 1.) Count the number of ions per DG cell and divide it by element-volume -> ion density n_i
! 2.) Count the number of neutrals per DG cell and divide it by element-volume -> neutral density n_n
! 3.) Calculate the ionization degree: alpha = n_i/(n_i + n_n)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Analyze_Vars  ,ONLY:IonizationCell,QuasiNeutralityCell,NeutralDensityCell,ElectronDensityCell,IonDensityCell
USE MOD_Particle_Analyze_Vars  ,ONLY:ChargeNumberCell
USE MOD_Preproc                ,ONLY:PP_nElems
USE MOD_Particle_Mesh_Vars     ,ONLY:GEO
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! Nullify
IonizationCell      = 0.
QuasiNeutralityCell = 0.

! Loop over all elements
DO iElem=1,PP_nElems
  ASSOCIATE(&
    Q   => QuasiNeutralityCell(iElem) ,& ! Quasi neutral condition approximation
    n_e => ElectronDensityCell(iElem) ,& ! Electron number density (cell average)
    X   => IonizationCell(iElem)      ,& ! Ionization degree (cell average)
    n_i => IonDensityCell(iElem)      ,& ! Ion number density (cell average)
    n_n => NeutralDensityCell(iElem)   ) ! Neutral number density (cell average)

    IF(ABS(n_i + n_n).LE.0.0)THEN ! no particles in cell
      X = 0.0
      Q = 0.0
    ELSE
      ! 0.  Set degree of ionization: X = n_i / (n_i + n_n)
      X  = n_i / (n_i + n_n)

      ! Set quasi neutrality between zero and unity depending on which density is larger
      ! Quasi neutrality holds, when n_e ~ Z_i*n_i (electron density approximately equal to ion density multiplied with charge number)
      ! 1.  Calculate Z_i*n_i (Charge density cell average)
      Q = ChargeNumberCell(iElem) / GEO%Volume(iElem)

      ! 2.  Calculate the quasi neutrality parameter: should be near to 1 for quasi-neutrality
      IF(Q.GT.n_e)THEN
        ! 2.1  if Z_i*n_i > n_e -> calculate n_e/Z_i*n_i
        Q = n_e / Q
      ELSE
        ! 2.2  if Z_i*n_i < n_e -> calculate Z_i*n_i/n_e
        IF(ABS(n_e).GT.0.0)THEN
          Q = Q / n_e
        ELSE
          Q = 0.0
        END IF
      END IF
    END IF
  END ASSOCIATE
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateIonizationCell


SUBROUTINE CalculatePlasmaParameter()
!===================================================================================================================================
! Calculate the points per Debye length for each cell
! PointsPerDebyeLength: PPD = (p+1)*lambda_D/L_cell
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars          ,ONLY: PI
USE MOD_Preproc               ,ONLY: PP_nElems
USE MOD_Particle_Analyze_Vars ,ONLY: DebyeLengthCell,ElectronDensityCell,PlasmaParameterCell
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! loop over all elements
DO iElem=1,PP_nElems
  IF((DebyeLengthCell(iElem).GT.0.0).AND.(ElectronDensityCell(iElem).GT.0.0))THEN
    PlasmaParameterCell(iElem) = (4.0/3.0) * PI * ElectronDensityCell(iElem) * (DebyeLengthCell(iElem)**3)
  ELSE
    PlasmaParameterCell(iElem) = 0.0
  END IF
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePlasmaParameter


SUBROUTINE CalculatePartElemData()
!===================================================================================================================================
! use the plasma frequency per cell to estimate the pic time step
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Analyze_Vars  ,ONLY:CalcPlasmaFrequency,CalcPICTimeStep,CalcElectronIonDensity
USE MOD_Particle_Analyze_Vars  ,ONLY:CalcElectronTemperature,CalcDebyeLength,CalcIonizationDegree,CalcPointsPerDebyeLength
USE MOD_Particle_Analyze_Vars  ,ONLY:CalcPlasmaParameter
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! electron density
IF(CalcElectronIonDensity) CALL CalculateElectronIonDensityCell()

! Ionization degree: n_i / (n_i + n_n)
! ion density versus sum of ion and neutral density
IF(CalcIonizationDegree) CALL CalculateIonizationCell()

! electron temperature
IF(CalcElectronTemperature) CALL CalculateElectronTemperatureCell()

! plasma frequency
IF(CalcPlasmaFrequency) CALL CalculatePlasmaFrequencyCell()

! Debye length
IF(CalcDebyeLength) CALL CalculateDebyeLengthCell()

! Plasma parameter: 4/3 * pi * n_e * lambda_D^3
IF(CalcPlasmaParameter) CALL CalculatePlasmaParameter()

! PIC time step
IF(CalcPICTimeStep) CALL CalculatePICTimeStepCell()

! PointsPerDebyeLength: PPD = (p+1)*lambda_D/L_cell
IF(CalcPointsPerDebyeLength) CALL CalculatePPDCell()

END SUBROUTINE CalculatePartElemData

#ifdef CODE_ANALYZE
!===================================================================================================================================
!> Calculate the analytical position and velocity depending on the pre-defined function
!===================================================================================================================================
SUBROUTINE CalcAnalyticalParticleState(t,PartStateAnalytic)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: PI
USE MOD_PreProc
USE MOD_PICInterpolation_Vars ,ONLY: AnalyticInterpolationType,AnalyticInterpolationSubType,AnalyticInterpolationP
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t                        !< simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: PartStateAnalytic(1:6)   !< analytic position and velocity
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL    :: p
REAL    :: gamma_0
REAL    :: phi_0
REAL    :: Theta
REAL    :: beta
!===================================================================================================================================
PartStateAnalytic=0. ! default

! select the analytical solution
SELECT CASE(AnalyticInterpolationType)
CASE(1)
  SELECT CASE(AnalyticInterpolationSubType)
  CASE(1,2)
    ASSOCIATE( p       => AnalyticInterpolationP , &
               Theta_0 => -PI/2.0                     , &
               t       => t - TEnd/2. )
               !t       => t )
      ! gamma
      gamma_0 = SQRT(ABS(p*p - 1.))

      ! angle
      Theta   = -2.*ATAN( SQRT((1.+p)/(1.-p)) * TANH(0.5*gamma_0*t) ) + Theta_0

      ! x-pos
      PartStateAnalytic(1) = LOG(-SIN(Theta) + p )

      ! y-pos
      PartStateAnalytic(2) = p*t + Theta - Theta_0
    END ASSOCIATE
  CASE(3)
    ASSOCIATE( p       => AnalyticInterpolationP , &
               Theta_0 => -PI/2.0                      &
                )
      ! gamma
      gamma_0 = SQRT(ABS(p*p - 1.))

      ! angle
      Theta   = -2.*ATAN( SQRT((p+1.)/(p-1.)) * TAN(0.5*gamma_0*t) ) -2.*PI*REAL(NINT((gamma_0*t)/(2.*PI))) + Theta_0

      ! x-pos
      PartStateAnalytic(1) = LOG(-SIN(Theta) + p )

      ! y-pos
      PartStateAnalytic(2) = p*t + Theta - Theta_0
    END ASSOCIATE
  CASE(11,21) ! old CASE(1,2)
    ASSOCIATE( p       => AnalyticInterpolationP , &
          Theta_0 => 0.d0 ) !0.785398163397448d0    )
      beta = ACOS(p)
      !beta = ASIN(-p)
      ! phase shift
      phi_0   = ATANH( (1./TAN(beta/2.)) * TAN(Theta_0/2.) )
      ! angle
      Theta   = -2.*ATANH( TAN(beta/2.) * TANH(0.5*t*SIN(beta)-phi_0) )
      Theta   = -2.*ATANH( TAN(beta/2.) * TANH(0.5*SIN(beta*t)-phi_0) )
      ! x-pos
      PartStateAnalytic(1) = LOG((COS(Theta)-p)/(COS(Theta_0)-p))
      ! y-pos
      PartStateAnalytic(2) = p*t - (Theta-Theta_0)
    END ASSOCIATE
  CASE(31) ! old CASE(3)
    ASSOCIATE( p       => AnalyticInterpolationP , &
          Theta_0 => 0.d0                   )
      gamma_0 = SQRT(p*p-1.)
      ! phase shift
      phi_0   = ATAN( (gamma_0/(p-1.)) * TAN(Theta_0/2.) )
      ! angle
      Theta   = 2.*ATAN( SQRT((p-1)/(p+1)) * TAN(0.5*gamma_0*t - phi_0) ) + 2*Pi*REAL(NINT((t*gamma_0)/(2*Pi) - phi_0/Pi))
      ! x-pos
      PartStateAnalytic(1) = LOG((COS(Theta)-p)/(COS(Theta_0)-p))
      ! y-pos
      PartStateAnalytic(2) = p*t - (Theta-Theta_0)
    END ASSOCIATE
    !WRITE (*,*) "PartStateAnalytic =", PartStateAnalytic
    !read*
  END SELECT
END SELECT

END SUBROUTINE CalcAnalyticalParticleState


!===================================================================================================================================
!> Calculates "running" L_2 norms
!> running means: use the old L_2 error from the previous iteration in order to determine the L_2 error over time (simulation time)
!>
!> -------------------------------------------------------------------------
!> OLD METHOD: assuming constant timestep (ignoring the total time tEnd -> Delta t = tEnd / Niter)
!> L_2(t) = SQRT( ( L_2(t-1)^2 * (iter-1) + delta(t)^2 ) / iter )
!>
!> -------------------------------------------------------------------------
!> NEW METHOD: assuming variable timestep
!> L_2(t) = SQRT(  L_2(t-1)^2   +   (t - t_old) * delta(t)^2  )
!>
!> t     : simulation time
!> t_old : simulation time of the last iteration
!> L_2   : error norm
!> delta : difference numerical to analytical solution
!> iter  : simulation iteration counter
!===================================================================================================================================
SUBROUTINE CalcErrorParticle(t,iter,PartStateAnalytic)
! MODULES
USE MOD_PICInterpolation_Vars ,ONLY: L_2_Error_Part,L_2_Error_Part_time
USE MOD_Particle_Vars         ,ONLY: PartState, PDM
! OLD METHOD: considering TEnd:
! USE MOD_TimeDisc_Vars         ,ONLY: TEnd
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)    :: iter                     !< simulation iteration counter
REAL,INTENT(IN)               :: t                        !< simulation time
REAL,INTENT(INOUT)            :: PartStateAnalytic(1:6)   !< analytic position and velocity
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: i,j
!===================================================================================================================================
! Get analytic particle position
CALL CalcAnalyticalParticleState(t,PartStateAnalytic)

! Depending on the iteration counter, set the L_2 error (re-use the value in the next loop)
IF(iter.LT.1)THEN ! first iteration
  L_2_Error_Part(1:6) = 0.
  L_2_Error_Part_time = 0.
ELSE
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      DO j = 1, 6
        ! OLD METHOD: original
        ! L_2_Error_Part(j) = SQRT( ( (L_2_Error_Part(j))**2*REAL(iter-1) + (PartStateAnalytic(j)-PartState(i,j))**2 )/ REAL(iter))

        ! OLD METHOD: considering TEnd
        ! L_2_Error_Part(j) = SQRT( Tend * ( (L_2_Error_Part(j))**2*REAL(iter-1) + (PartStateAnalytic(j)-PartState(i,j))**2 ) &
        !                      / REAL(iter))

        ! NEW METHOD: considering variable time step
        L_2_Error_Part(j) = SQRT(  (L_2_Error_Part(j))**2 + (t-L_2_Error_Part_time)*(PartStateAnalytic(j)-PartState(i,j))**2 )
      END DO ! j = 1, 6
      L_2_Error_Part_time = t
    ELSE
      L_2_Error_Part(1:6) = -1.0
    END IF
  END DO
END IF

END SUBROUTINE CalcErrorParticle


!===================================================================================================================================
!> Calculate the analytical position and velocity depending on the pre-defined function
!===================================================================================================================================
SUBROUTINE AnalyticParticleMovement(time,iter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars           ,ONLY: OutputErrorNorms
USE MOD_Particle_Analyze_Vars  ,ONLY: TrackParticlePosition
USE MOD_PICInterpolation_Vars  ,ONLY: L_2_Error_Part
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: time                        !< simulation time
INTEGER(KIND=8),INTENT(IN)    :: iter                        !< iteration
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: PartStateAnalytic(1:6)   !< analytic position and velocity
CHARACTER(LEN=40)             :: formatStr
!===================================================================================================================================

CALL CalcErrorParticle(time,iter,PartStateAnalytic)
IF(PartMPI%MPIRoot.AND.OutputErrorNorms) THEN
  WRITE(UNIT_StdOut,'(A13,ES16.7)')' Sim time  : ',time
  WRITE(formatStr,'(A5,I1,A7)')'(A13,',6,'ES16.7)'
  WRITE(UNIT_StdOut,formatStr)' L2_Part   : ',L_2_Error_Part
  OutputErrorNorms=.FALSE.
END IF
IF(TrackParticlePosition) CALL WriteParticleTrackingDataAnalytic(time,iter,PartStateAnalytic) ! new function

END SUBROUTINE AnalyticParticleMovement
#endif /*CODE_ANALYZE*/


SUBROUTINE FinalizeParticleAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars ,ONLY: ParticleAnalyzeInitIsDone,DebyeLengthCell,PICTimeStepCell &
                                    ,ElectronTemperatureCell,ElectronDensityCell,PlasmaFrequencyCell,PPSCell,PPSCellEqui
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ParticleAnalyzeInitIsDone = .FALSE.
SDEALLOCATE(DebyeLengthCell)
SDEALLOCATE(PICTimeStepCell)
SDEALLOCATE(ElectronDensityCell)
SDEALLOCATE(ElectronTemperatureCell)
SDEALLOCATE(PlasmaFrequencyCell)
SDEALLOCATE(PPSCell)
SDEALLOCATE(PPSCellEqui)
END SUBROUTINE FinalizeParticleAnalyze
#endif /*PARTICLES*/


END MODULE MOD_Particle_Analyze
