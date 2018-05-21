#include "boltzplatz.h"

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

INTERFACE CalcKineticEnergy
  MODULE PROCEDURE CalcKineticEnergy
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

PUBLIC:: InitParticleAnalyze, FinalizeParticleAnalyze!, CalcPotentialEnergy
PUBLIC:: CalcKineticEnergy, CalcEkinPart, AnalyzeParticles, PartIsElectron
PUBLIC:: CalcPowerDensity
PUBLIC:: CalculatePartElemData
#if (PP_TimeDiscMethod==42)
PUBLIC :: ElectronicTransition, WriteEletronicTransition
#endif
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

CALL prms%CreateIntOption(      'Part-AnalyzeStep'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Analyze is performed each Nth time step','1') 
CALL prms%CreateLogicalOption(  'CalcPotentialEnergy', 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Flag to calculate Potential Energy.','.FALSE.')
CALL prms%CreateLogicalOption(  'PIC-VerifyCharge'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Validate the charge after each deposition'//&
                                                       'and produces an output in std.out','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcDebyeLength'   ,  'Flag to compute the Debye length (min and max) in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPICTimeStep'   ,  'Flag to compute the HDG time step (min and max) in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcElectronTemperature', 'Flag to compute the electron temperature in each cell','.FALSE.')
!CALL prms%CreateLogicalOption(  'ElectronTemperatureIsMaxwell', 'Flag if  electron temperature is assumed to be Maxwellian in each cell','.TRUE.')
CALL prms%CreateLogicalOption(  'CalcElectronDensity', 'Flag to compute the electron density in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPlasmaFrequency', 'Flag to compute the electron frequency in each cell','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcCharge'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Flag to compute the whole deposited charge,'//&
                                                       ' absolute and relative charge error','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcKineticEnergy'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate Kinetic Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcInternalEnergy' , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate Internal Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcTemp'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate Translational temperature.'&
                                                     ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPartBalance'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the Particle Power Balance'//&
                                                       '- input and outflow energy of all particles','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcVelos'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate thermal and flow velocities.'//&
                                                       'if CalcVelos = T VelocityDirections = (/[int],[int],[int],[int]/)  '//&
                                                       'Switching dimensions for CalcVelos on (1) or off (0)\n'//&
                                                       '(/v_x,v_y,v_z,|v|/) ','.FALSE.')
CALL prms%CreateIntArrayOption( 'VelocityDirections' , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'x,y,z,abs -> 0/1 = T/F. (please note: CalcVelos)'&
                                                     ,'1 , 1 , 1 , 1')
CALL prms%CreateLogicalOption(  'Part-TrackPosition' , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Track particle position','.FALSE.')
CALL prms%CreateLogicalOption(  'printDiff'          , 'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateRealOption(     'printDiffTime'      , 'TODO-DEFINE-PARAMETER','12.')
CALL prms%CreateRealArrayOption('printDiffVec'       , 'TODO-DEFINE-PARAMETER','0. , 0. , 0. , 0. , 0. , 0.')
CALL prms%CreateLogicalOption(  'CalcNumSpec'        , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate species count.','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcCollRates'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the collision rates per '//&
                                                       'collision pair','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcReacRates'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the reaction rate per reaction'&
                                                     ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcSurfNumSpec'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the number of simulated'//&
                                                       'particles per species on surfaces','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcSurfCoverage'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the surface coverages for'//&
                                                       'each species','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcAccomodation'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the surface accomodation coefficient'&
                                                     ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcEvaporation'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate rate of evaporation [kg/s]','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcAdsorbRates'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calcualte the adsorption probabilities of species'&
                                                     ,'.FALSE.')
CALL prms%CreateLogicalOption(  'CalcSurfRates'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Calculate the surface reaction rate per reaction'//&
                                                       ' (k_r)','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcShapeEfficiency', 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Use efficiency methods for shape functions.'&
                                                     , '.FALSE.')
CALL prms%CreateStringOption(   'CalcShapeEfficiencyMethod'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Choose between "AllParts" and '//&
                                                       '"SomeParts", to either use all particles or a certain percentage'//&
                                                       ' (ShapeEfficiencyNumber) of the currently used particles','AllParts')
CALL prms%CreateIntOption(      'ShapeEfficiencyNumber'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Percentage of currently used particles is used.'&
                                                     ,'100')
CALL prms%CreateLogicalOption(  'IsRestart'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Flag, if the current calculation is a restart. '&
                                                     ,'.FALSE.')

END SUBROUTINE DefineParametersParticleAnalyze

SUBROUTINE InitParticleAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars          ,ONLY: DoAnalyze,CalcEpot
USE MOD_Particle_Analyze_Vars 
USE MOD_ReadInTools           ,ONLY: GETLOGICAL, GETINT, GETSTR, GETINTARRAY, GETREALARRAY, GETREAL
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_PICDepo_Vars          ,ONLY: DoDeposition
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_Vars             ,ONLY: Adsorption
#endif
USE MOD_IO_HDF5               ,ONLY: AddToElemData,ElementOut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: dir, VeloDirs_hilf(4)
!===================================================================================================================================
IF (ParticleAnalyzeInitIsDone) THEN
CALL abort(__STAMP__,&
'InitParticleAnalyse already called.',999,999.)
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

PartAnalyzeStep = GETINT('Part-AnalyzeStep','1')
IF (PartAnalyzeStep.EQ.0) PartAnalyzeStep = 123456789

DoAnalyze = .FALSE.
CalcEpot = GETLOGICAL('CalcPotentialEnergy','.FALSE.')
IF(CalcEpot) DoAnalyze = .TRUE.
! only verifycharge and CalcCharge if particles are deposited onto the grid
DoVerifyCharge= .FALSE.
CalcCharge = .FALSE.
IF(DoDeposition) THEN
  DoVerifyCharge = GETLOGICAL('PIC-VerifyCharge','.FALSE.')
  CalcCharge = GETLOGICAL('CalcCharge','.FALSE.')
  IF(CalcCharge) DoAnalyze = .TRUE. 
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' Deposition is switched of. VerifyCharge and CalcCharge are deactivated!'
END IF

!--------------------------------------------------------------------------------------------------------------------
! get derived particle properties 
! (Note that for IMD/TTM initialization these values are calculated from the TTM grid values)
!--------------------------------------------------------------------------------------------------------------------
! Debye Length
CalcDebyeLength       = GETLOGICAL('CalcDebyeLength','.FALSE.')
IF(CalcDebyeLength)THEN
  ALLOCATE( DebyeLengthCell(1:PP_nElems) )
  DebyeLengthCell=0.0
  CALL AddToElemData(ElementOut,'DebyeLengthCell',RealArray=DebyeLengthCell(1:PP_nElems))
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
CalcElectronDensity   = GETLOGICAL('CalcElectronDensity','.FALSE.')
IF(CalcDebyeLength.OR.CalcPlasmaFrequency) CalcElectronDensity=.TRUE.
IF(CalcElectronDensity) THEN
  ALLOCATE( ElectronDensityCell(1:PP_nElems) )
  ElectronDensityCell=0.0
  CALL AddToElemData(ElementOut,'ElectronDensityCell',RealArray=ElectronDensityCell(1:PP_nElems))
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


CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
CalcEint = GETLOGICAL('CalcInternalEnergy','.FALSE.')
CalcTemp = GETLOGICAL('CalcTemp','.FALSE.')
IF(CalcTemp.OR.CalcEint) DoAnalyze = .TRUE.
IF(CalcEkin) DoAnalyze = .TRUE.
IF(nSpecies.GT.1) THEN
  nSpecAnalyze = nSpecies + 1
ELSE
  nSpecAnalyze = 1
END IF
CalcPartBalance = GETLOGICAL('CalcPartBalance','.FALSE.')
IF (CalcPartBalance) THEN
  DoAnalyze = .TRUE.
  SDEALLOCATE(nPartIn)
  SDEALLOCATE(nPartOut)
  SDEALLOCATE(PartEkinIn)
  SDEALLOCATE(PartEkinOut)
  ALLOCATE( nPartIn(nSpecies)     &
          , nPartOut(nSpecies)    &
          , PartEkinOut(nSpecies) &
          , PartEkinIn(nSpecies)  )
  nPartIn=0
  nPartOut=0
  PartEkinOut=0.
  PartEkinIn=0.
#if defined(LSERK) || defined(ROS) || defined(IMPA) 
  SDEALLOCATE( nPartInTmp)
  SDEALLOCATE( PartEkinInTmp)
  ALLOCATE( nPartInTmp(nSpecies)     &
          , PartEkinInTmp(nSpecies)  )
  PartEkinInTmp=0.
  nPartInTmp=0
#endif
END IF
TrackParticlePosition = GETLOGICAL('Part-TrackPosition','.FALSE.')
IF(TrackParticlePosition)THEN
  printDiff=GETLOGICAL('printDiff','.FALSE.')
  IF(printDiff)THEN
    printDiffTime=GETREAL('printDiffTime','12.')
    printDiffVec=GETREALARRAY('printDiffVec',6,'0.,0.,0.,0.,0.,0.')
  END IF
END IF
CalcNumSpec   = GETLOGICAL('CalcNumSpec','.FALSE.')
CalcCollRates = GETLOGICAL('CalcCollRates','.FALSE.')
CalcReacRates = GETLOGICAL('CalcReacRates','.FALSE.')
IF(CalcNumSpec.OR.CalcCollRates.OR.CalcReacRates) DoAnalyze = .TRUE.
CalcVelos = GETLOGICAL('CalcVelos','.FALSE')
IF (CalcVelos) THEN
  DoAnalyze=.TRUE.
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
#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
CalcSurfNumSpec = GETLOGICAL('CalcSurfNumSpec','.FALSE.')
CalcSurfCoverage = GETLOGICAL('CalcSurfCoverage','.FALSE.')
CalcAccomodation = GETLOGICAL('CalcAccomodation','.FALSE.')
#if (PP_TimeDiscMethod==42)
CalcEvaporation = GETLOGICAL('CalcEvaporation','.FALSE.')
IF (CalcEvaporation) DoAnalyze = .TRUE.
CalcAdsorbRates = GETLOGICAL('CalcAdsorbRates','.FALSE.')
CalcSurfRates = GETLOGICAL('CalcSurfRates','.FALSE.')
IF(CalcSurfNumSpec.OR.CalcSurfRates.OR.CalcSurfCoverage.OR.CalcAccomodation.OR.Adsorption%TPD.OR.CalcAdsorbRates) &
    DoAnalyze = .TRUE.
IF (Adsorption%TPD.AND.((.NOT.CalcSurfRates))) CalcSurfRates = .TRUE.
#else
IF(CalcSurfNumSpec.OR.CalcSurfCoverage.OR.CalcAccomodation) DoAnalyze = .TRUE.
#endif
#endif
CalcShapeEfficiency = GETLOGICAL('CalcShapeEfficiency','.FALSE.')
IF (CalcShapeEfficiency) THEN
  DoAnalyze = .TRUE.
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

IsRestart = GETLOGICAL('IsRestart','.FALSE.')

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
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars,          ONLY: DoAnalyze,CalcEpot
USE MOD_Particle_Analyze_Vars!,ONLY: ParticleAnalyzeInitIsDone,CalcCharge,CalcEkin,IsRestart
USE MOD_PARTICLE_Vars,         ONLY: nSpecies, BoltzmannConst
USE MOD_DSMC_Vars,             ONLY: CollInf, useDSMC, CollisMode, ChemReac
USE MOD_Restart_Vars,          ONLY: DoRestart
USE MOD_AnalyzeField,          ONLY: CalcPotentialEnergy,CalcPotentialEnergy_Dielectric
USE MOD_DSMC_Vars,             ONLY: DSMC
USE MOD_Dielectric_Vars,       ONLY: DoDielectric
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506))
USE MOD_TimeDisc_Vars          ,ONLY: iter
USE MOD_DSMC_Analyze           ,ONLY: CalcMeanFreePath
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
#endif
USE MOD_PIC_Analyze            ,ONLY: CalcDepositedCharge
#ifdef MPI
#if USE_LOADBALANCE
USE MOD_LoadBalance_tools      ,ONLY: LBStartTime, LBSplitTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
#if ( PP_TimeDiscMethod ==42)
USE MOD_DSMC_Vars              ,ONLY: Adsorption,BGGas
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
#endif
USE MOD_Particle_Analyze_Vars  ,ONLY: ChemEnergySum
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
INTEGER             :: unit_index, iSpec, OutputCounter
INTEGER(KIND=8)     :: SimNumSpec(nSpecAnalyze)
REAL                :: WEl, WMag, NumSpec(nSpecAnalyze)
REAL                :: Ekin(nSpecAnalyze), Temp(nSpecAnalyze)
REAL                :: IntEn(nSpecAnalyze,3),IntTemp(nSpecies,3),TempTotal(nSpecAnalyze), Xi_Vib(nSpecies), Xi_Elec(nSpecies)
REAL                :: MaxCollProb, MeanCollProb, ETotal, totalChemEnergySum, MeanFreePath
#ifdef MPI
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42)
REAL                :: sumMeanCollProb
#endif
REAL                :: RECBR(nSpecies),RECBR1
INTEGER             :: RECBIM(nSpecies)
#endif /*MPI*/
REAL, ALLOCATABLE   :: CRate(:), RRate(:)
#if (PP_TimeDiscMethod ==42)
INTEGER             :: ii, iunit, iCase, iTvib,jSpec
CHARACTER(LEN=64)   :: DebugElectronicStateFilename
CHARACTER(LEN=350)  :: hilf
REAL                :: NumSpecTmp(nSpecAnalyze)
REAL                :: Adsorptionrate(nSpecies), Desorptionrate(nSpecies)
REAL,ALLOCATABLE    :: SurfReactRate(:), AdsorptionReactRate(:), AdsorptionActE(:), SurfaceActE(:)
#endif
#if (PP_TimeDiscMethod ==42) || (PP_TimeDiscMethod ==4)
INTEGER(KIND=8)     :: WallNumSpec(nSpecies), WallNumSpec_SurfDist(nSpecies)
INTEGER             :: iCov
REAL                :: WallCoverage(nSpecies), Accomodation(nSpecies)
REAL                :: EvaporationRate(nSpecies)
INTEGER             :: SurfCollNum(nSpecies),AdsorptionNum(nSpecies),DesorptionNum(nSpecies)
#endif
REAL                :: PartVtrans(nSpecies,4) ! macroscopic velocity (drift velocity) A. Frohn: kinetische Gastheorie
REAL                :: PartVtherm(nSpecies,4) ! microscopic velocity (eigen velocity) PartVtrans + PartVtherm = PartVtotal
INTEGER             :: dir
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
  IF ( DoRestart ) THEN
    isRestart = .true.
  END IF
  IF (.NOT.DoAnalyze) RETURN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
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
#ifdef MPI
!#ifdef PARTICLES
  IF (PartMPI%MPIRoot) THEN
#endif    /* MPI */
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
      ELSE IF (Adsorption%TPD) THEN
          iCov = INT(Adsorption%Coverage(1,1,1,1)*1000)
          WRITE( hilf, '(I4.4)') iCov
        outfile = 'Database_Cov_'//TRIM(hilf)//'.csv'
      ELSE
        outfile = 'Database.csv'
      END IF
#else
      outfile = 'Database.csv'
#endif

      IF (isRestart .and. FILEEXISTS(outfile)) THEN
        OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
        !CALL FLUSH (unit_index)
      ELSE
        OPEN(unit_index,file=TRIM(outfile))
        !CALL FLUSH (unit_index)
        !--- insert header
        WRITE(unit_index,'(A6,A5)',ADVANCE='NO') 'TIME', ' '
        IF (CalcNumSpec) THEN
          DO iSpec = 1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A12,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPart-Spec-', iSpec,' '
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
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A14,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartIn-Spec-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A15,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartOut-Spec-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcEpot) THEN 
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-W-El',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-W-Mag',' '
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcEkin) THEN
          DO iSpec=1, nSpecAnalyze
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Ekin-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
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
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A8,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinIn-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
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
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506))
        IF (CollisMode.GT.1) THEN
          IF(CalcEint) THEN
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
            IF (DSMC%ElectronicModel) THEN
              DO iSpec = 1, nSpecAnalyze
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EElec',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
            END IF
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-ETotal',' '
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
#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
        IF (DSMC%WallModel.EQ.3) THEN
          IF (CalcSurfNumSpec) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nSimPart-Wall-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nSurfPart-Wall-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
          IF (CalcSurfCoverage) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Surf-Cov-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
#if (PP_TimeDiscMethod==42)
          IF (CalcAccomodation) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Alpha-Spec', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
          IF (CalcAdsorbRates) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nSurfColl-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-N_Ads-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Prob_adsorption-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-P_Molec-Adsorb-Spec-',iSpec,' '
                OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P_Dissoc-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1, Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P_ER-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            DO iSpec = 1, nSpecies
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-diss-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-ER-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
          END IF
          IF (CalcSurfRates) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-N_Des-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-P_Des-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-P-SurfDesorb-Molec-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
              DO iCase = 1, Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P-SurfDissoc-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1, Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P-SurfLH-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            DO iCase = 1, Adsorption%NumOfExchReact
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-P-Surfexch-Case-', iCase,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-E-Desorb-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-Diss-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-LH-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            DO iCase = 1,Adsorption%NumOfExchReact
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') &
                  OutputCounter,'-E-Exch-Reaction-', iCase,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
          IF (Adsorption%TPD) THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-WallTemp',' '
          END IF
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcEvaporation) THEN
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' -Evap-Mass-Spec', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
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
#endif
        END IF
#endif
        WRITE(unit_index,'(A1)') ' '
      END IF
    END IF
#ifdef MPI
  END IF
#endif    /* MPI */

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================
  ! computes the real and simulated number of particles
  CALL CalcNumPartsofSpec(NumSpec,SimNumSpec)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate total temperature of each molecular species (Laux, p. 109)
  IF(CalcEkin) CALL CalcKineticEnergy(Ekin)
  IF(CalcTemp.OR.CalcEint.OR.DSMC%CalcQualityFactors) THEN
    CALL CalcTemperature(NumSpec,Temp,IntTemp,IntEn,TempTotal,Xi_Vib,Xi_Elec) ! contains MPI Communication
    IF(CalcEint.AND.(CollisMode.GT.1)) THEN
      CALL CalcIntTempsAndEn(NumSpec,IntTemp,IntEn)
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
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506))
    IF(DSMC%CalcQualityFactors) THEN
      MeanFreePath = 0.0
#ifdef MPI
      IF((iter.GT.0).AND.PartMPI%MPIRoot) THEN
#else
      IF((iter.GT.0)) THEN
#endif
        IF(TempTotal(nSpecAnalyze).GT.0.0) MeanFreePath = CalcMeanFreePath(NumSpec(1:nSpecies), NumSpec(nSpecAnalyze), &
                                                              GEO%MeshVolume, SpecDSMC(1)%omegaVHS, TempTotal(nSpecAnalyze))
      END IF
    END IF
#endif
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Determine the maximal collision probability for whole reservoir and mean collision probability (only for one cell)
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506))
  IF(DSMC%CalcQualityFactors) THEN
    MaxCollProb = 0.0
    MeanCollProb = 0.0
    IF(iter.GT.0) THEN
      MaxCollProb = DSMC%CollProbMax
      IF(DSMC%CollProbMeanCount.GT.0) MeanCollProb = DSMC%CollProbMean / DSMC%CollProbMeanCount
    END IF
  END IF
#else
  MaxCollProb = 0.0
  MeanCollProb = 0.0
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! Other Analyze Routines
  IF(CalcCharge) CALL CalcDepositedCharge() ! mpi communication done in calcdepositedcharge
  IF(CalcEpot)THEN
    IF(DoDielectric)THEN
      CALL CalcPotentialEnergy_Dielectric(WEl,WMag)
    ELSE
      CALL CalcPotentialEnergy(WEl,WMag)
    END IF
  END IF
  IF(TrackParticlePosition) CALL TrackingParticlePosition(time)
  IF(CalcVelos) CALL CalcVelocities(PartVtrans, PartVtherm,NumSpec,SimNumSpec)
!===================================================================================================================================
! MPI Communication for values which are not YET communicated
! all routines ABOVE contains the required MPI-Communication
!===================================================================================================================================
#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    IF (CalcPartBalance)THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,nPartIn(:)    ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,nPartOUt(:)   ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinIn(:) ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinOut(:),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42)
    IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
      ! Determining the maximal (MPI_MAX) and mean (MPI_SUM) collision probabilities
      CALL MPI_REDUCE(MPI_IN_PLACE,MaxCollProb,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      CALL MPI_REDUCE(MeanCollProb,sumMeanCollProb,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
      MeanCollProb = sumMeanCollProb / REAL(PartMPI%nProcs)
    END IF
#endif
  ELSE ! no Root
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    IF (CalcPartBalance)THEN
      CALL MPI_REDUCE(nPartIn,RECBIM   ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(nPartOut,RECBIM  ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PartEkinIn,RECBR ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PartEkinOut,RECBR,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42)
    IF((iter.GT.0).AND.(DSMC%CalcQualityFactors)) THEN
      CALL MPI_REDUCE(MaxCollProb,RECBR1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, PartMPI%COMM, IERROR)
      CALL MPI_REDUCE(MeanCollProb,sumMeanCollProb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, PartMPI%COMM, IERROR)
    END IF
#endif
  END IF
#endif /*MPI*/
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTANALYZE,tLBStart)
#endif /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==1000)
  IF (CollisMode.GT.1) CALL CalcIntTempsAndEn(NumSpec,IntTemp,IntEn,Xi_Vib,Xi_Elec)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate the collision rates and reaction rate coefficients (Arrhenius-type chemistry)
#if (PP_TimeDiscMethod==42)
  IF(CalcCollRates) CALL CollRates(CRate)
  IF(CalcReacRates) THEN
    IF ((CollisMode.EQ.3).AND.(iter.GT.0)) THEN
      NumSpecTmp = NumSpec
      IF(BGGas%BGGasSpecies.NE.0) THEN
        NumSpecTmp(BGGas%BGGasSpecies) = (BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor)
        IF(nSpecAnalyze.GT.1)THEN
          NumSpecTmp(nSpecAnalyze) = NumSpecTmp(nSpecAnalyze)+NumSpecTmp(BGGas%BGGasSpecies)
        END IF
      END IF
      CALL ReacRates(RRate, NumSpecTmp)
    END IF
  END IF
#endif
#if (PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)
IF (DSMC%WallModel.EQ.3) THEN
  IF (CalcSurfNumSpec.OR.CalcSurfCoverage) CALL GetWallNumSpec(WallNumSpec,WallCoverage,WallNumSpec_SurfDist)
#if (PP_TimeDiscMethod==42)
  IF (CalcAccomodation) CALL GetAccCoeff(Accomodation)
  IF (CalcAdsorbRates) THEN
    SDEALLOCATE(AdsorptionReactRate)
    SDEALLOCATE(AdsorptionActE)
    ALLOCATE(AdsorptionReactRate(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(AdsorptionActE(1:nSpecies*(Adsorption%ReactNum)))
    CALL GetAdsRates(Adsorptionrate,SurfCollNum,AdsorptionNum,AdsorptionReactRate,AdsorptionActE)
  ELSE
    IF(SurfMesh%SurfOnProc)THEN
      DO iSpec = 1,nSpecies
        Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
      END DO
    END IF
  END IF
  IF (CalcSurfRates) THEN
    SDEALLOCATE(SurfReactRate)
    SDEALLOCATE(SurfaceActE)
    ALLOCATE(SurfReactRate(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(SurfaceActE(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    CALL GetSurfRates(Desorptionrate,DesorptionNum,SurfReactRate,SurfaceActE)
  END IF
#endif
END IF
#endif
#if (PP_TimeDiscMethod==42)
IF (CalcEvaporation) CALL GetEvaporationRate(EvaporationRate)
#endif /*PP_TimeDiscMethod==42*/
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
    IF (CalcCharge) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartCharge(1)
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartCharge(2)
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartCharge(3)
    END IF
    IF (CalcPartBalance) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartIn(iSpec))
      END DO
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartOut(iSpec))
      END DO
    END IF
    IF (CalcEpot) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') WEl
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') WMag
    END IF
    IF (CalcEkin) THEN
      DO iSpec=1, nSpecAnalyze
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Ekin(iSpec)
      END DO
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
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinIn(iSpec)
      END DO
      DO iSpec=1, nSpecies
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
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==300||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506))
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
    IF(DSMC%CalcQualityFactors) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') MeanCollProb
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') MaxCollProb
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') MeanFreePath
    END IF
#endif
#if ((PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4))
! output for adsorption
    IF (DSMC%WallModel.EQ.3) THEN
      IF (CalcSurfNumSpec) THEN
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I18.1)',ADVANCE='NO') WallNumSpec(iSpec)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I18.1)',ADVANCE='NO') WallNumSpec_SurfDist(iSpec)
        END DO
      END IF
      IF (CalcSurfCoverage) THEN
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') WallCoverage(iSpec)
        END DO
      END IF
#if (PP_TimeDiscMethod==42)
      IF (CalcAccomodation) THEN
        DO iSpec = 1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Accomodation(iSpec)
        END DO
      END IF
      IF (CalcAdsorbRates) THEN
        DO iSpec = 1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I18.1)',ADVANCE='NO') SurfCollNum(iSpec)
        END DO
        DO iSpec = 1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I18.1)',ADVANCE='NO') AdsorptionNum(iSpec)
        END DO
        DO iSpec = 1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Adsorptionrate(iSpec)
        END DO
        DO iCase = 1, nSpecies*(Adsorption%ReactNum+1)
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') AdsorptionReactRate(iCase)
        END DO
        DO iCase = 1, nSpecies*Adsorption%ReactNum
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') AdsorptionActE(iCase)
        END DO
      END IF
      IF (CalcSurfRates) THEN
        DO iSpec = 1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I18.1)',ADVANCE='NO') DesorptionNum(iSpec)
        END DO
        DO iSpec = 1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Desorptionrate(iSpec)
        END DO
        DO iCase = 1, nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') SurfReactRate(iCase)
        END DO
        DO iCase = 1, nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') SurfaceActE(iCase)
        END DO
      END IF
      IF (Adsorption%TPD) THEN
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Adsorption%TPD_Temp
      END IF
    END IF
    IF (CalcEvaporation) THEN
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') EvaporationRate(iSpec)
      END DO
    END IF
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
#endif /*(PP_TimeDiscMethod==42)*/
    END IF
#endif /*(PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)*/
    WRITE(unit_index,'(A1)') ' '
#ifdef MPI
  END IF
#endif    /* MPI */

  !104    FORMAT (e25.14)

!-----------------------------------------------------------------------------------------------------------------------------------
  IF( CalcPartBalance) CALL CalcParticleBalance()
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTANALYZE,tLBStart)
#endif /*USE_LOADBALANCE*/

#if ( PP_TimeDiscMethod ==42 )
! hard coded
! array not allocated
  IF ( DSMC%ElectronicModel ) THEN
  IF (DSMC%ReservoirSimuRate) THEN
    IF(Time.GT.0.) CALL ElectronicTransition( Time, NumSpec )
  END IF
END IF
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
#if ( PP_TimeDiscMethod ==42 )
  ! Debug Output for initialized electronic state
  IF ( DSMC%ElectronicModel .AND. DSMC%ReservoirSimuRate) THEN
    DO iSpec = 1, nSpecies
      IF ( SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
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
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
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


#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
SUBROUTINE GetWallNumSpec(WallNumSpec,WallCoverage,WallNumSpec_SurfDist)
!===================================================================================================================================
! Calculate number of wallparticles for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies, PDM, nSpecies, KeepWallParticles
USE MOD_Particle_Analyze_Vars
USE MOD_DSMC_Vars              ,ONLY: Adsorption, SurfDistInfo!, DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8), INTENT(OUT)    :: WallNumSpec(nSpecies),WallNumSpec_SurfDist(nSpecies)
REAL           , INTENT(OUT)    :: WallCoverage(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i, iSpec, iSurfSide, p, q, SideID, PartBoundID
REAL                            :: SurfPart
REAL                            :: Coverage(nSpecies)
#ifdef MPI
REAL                            :: RD(nSpecies)
INTEGER(KIND=8)                 :: IDR(nSpecies), ID1(nSpecies), ID2(nSpecies), ID3(nSpecies*2)
#endif /*MPI*/
INTEGER                         :: Coord, AdsorbID, Surfpos, SpecID
INTEGER                         :: adsorbates(nSpecies)
REAL                            :: SubWallNumSpec(nSpecies), WallNumSpec_tmp(2*nSpecies)
!===================================================================================================================================
WallNumSpec = 0
WallNumSpec_SurfDist = 0
SurfPart = 0.
Coverage(:) = 0.
WallCoverage(:) = 0.
WallNumSpec_tmp = 0.
SubWallNumSpec = 0.

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec=1,nSpecies
  DO iSurfSide=1,SurfMesh%nSides
    SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
    PartboundID = PartBound%MapToPartBC(BC(SideID))
    IF (PartBound%SolidCatalytic(PartboundID)) THEN
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        Coverage(iSpec) = Coverage(iSpec) + Adsorption%Coverage(p,q,iSurfSide,iSpec)
        IF ((.NOT.KeepWallParticles) .AND. CalcSurfNumSpec) THEN
          SurfPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
!          WallNumSpec(iSpec) = WallNumSpec(iSpec) + INT( Adsorption%Coverage(p,q,iSurfSide,iSpec) &
!              * SurfPart/Species(iSpec)%MacroParticleFactor)
          ! calculate number of adsorbates for each species
          adsorbates = 0
          DO Coord = 1,3
          DO AdsorbID = 1,SurfDistInfo(p,q,iSurfSide)%nSites(Coord)-SurfDistInfo(p,q,iSurfSide)%SitesRemain(Coord)
            Surfpos = SurfDistInfo(p,q,iSurfSide)%AdsMap(Coord)%UsedSiteMap(SurfDistInfo(p,q,iSurfSide)%SitesRemain(Coord)+AdsorbID)
            SpecID = SurfDistInfo(p,q,iSurfSide)%AdsMap(Coord)%Species(Surfpos)
            adsorbates(SpecID) = adsorbates(SpecID) + 1
          END DO
          END DO
          ! discret simulated particles on surface distribution
          WallNumSpec_SurfDist(iSpec) = WallNumSpec_SurfDist(iSpec) + adsorbates(iSpec)
          ! simulated (gas) particles from discret surface distribution
          SubWallNumSpec(iSpec) = SubWallNumSpec(iSpec) + REAL(adsorbates(iSpec)) / REAL(SurfDistInfo(p,q,iSurfSide)%nSites(3))&
              * SurfPart/Species(iSpec)%MacroParticleFactor
          ! simulated gas particles safed in temporary arrays            
          WallNumSpec_tmp(iSpec) = WallNumSpec_tmp(iSpec) + &
              ( SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) / SurfDistInfo(p,q,iSurfSide)%nSites(3) &
              * SurfPart / Species(iSpec)%MacroParticleFactor )
          WallNumSpec_tmp(iSpec+nSpecies) = WallNumSpec_tmp(iSpec+nSpecies) + SurfDistInfo(p,q,iSurfSide)%desorbnum_tmp(iSpec)&
              - SurfDistInfo(p,q,iSurfSide)%reactnum_tmp(iSpec)
        END IF
      END DO
    END DO
    END IF
  END DO
  END DO
  IF (CalcSurfCoverage) THEN
    WallCoverage(:) = Coverage(:) / (SurfMesh%nSides*nSurfSample*nSurfSample)
  END IF
END IF

#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    IF (CalcSurfNumSpec)  THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,SubWallNumSpec      ,nSpecies  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec_SurfDist,nSpecies  ,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec_tmp     ,nSpecies*2,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
    IF (CalcSurfCoverage) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,WallCoverage,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      WallCoverage = WallCoverage / REAL(SurfCOMM%nProcs)
    END IF
  ELSE
    IF (CalcSurfNumSpec) THEN
      CALL MPI_REDUCE(SubWallNumSpec      ,ID1,nSpecies  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(WallNumSpec_SurfDist,ID2,nSpecies  ,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(WallNumSpec_tmp     ,ID3,nSpecies*2,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
    IF (CalcSurfCoverage) CALL MPI_REDUCE(WallCoverage,RD,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
#endif /*MPI*/

  IF (KeepWallParticles.AND.CalcSurfNumSpec) THEN
    DO i=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i) .AND. PDM%ParticleAtWall(i)) THEN
        WallNumSpec(PartSpecies(i)) = WallNumSpec(PartSpecies(i)) + 1
      END IF
    END DO
#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    IF (CalcSurfNumSpec) CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec,nSpecies,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
  ELSE
    IF (CalcSurfNumSpec) CALL MPI_REDUCE(WallNumSpec ,IDR        ,nSpecies,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
#endif /*MPI*/
  ELSE
    WallNumSpec = INT(SubWallNumSpec)+INT(WallNumSpec_tmp(1:nSpecies))+INT(WallNumSpec_tmp(nSpecies+1:nSpecies*2))
  END IF

END SUBROUTINE GetWallNumSpec

#if (PP_TimeDiscMethod==42)
SUBROUTINE GetAccCoeff(Accomodation)
!===================================================================================================================================
! Calculate accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: Adsorption, DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: Accomodation(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec
#ifdef MPI
REAL                            :: AC(nSpecies)
#endif /*MPI*/
!===================================================================================================================================

Accomodation(:) = 0.
IF(SurfMesh%SurfOnProc)THEN
  IF (DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        Accomodation(iSpec) = Adsorption%AdsorpInfo(iSpec)%Accomodation / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        Accomodation(iSpec) = 0.
      END IF
    END DO
  ELSE IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        Accomodation(iSpec) = Adsorption%AdsorpInfo(iSpec)%Accomodation / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        Accomodation(iSpec) = 0.
      END IF
    END DO
  END IF
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Accomodation,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  Accomodation= Accomodation/ REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(Accomodation,AC          ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%Accomodation = 0.
  END DO
END IF

END SUBROUTINE GetAccCoeff


SUBROUTINE GetAdsRates(AdsorbRate,SurfCollNum,AdsorbNum,ReactRate,AdsorbActE)
!===================================================================================================================================
! Calculate adsorption, desorption and accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: Adsorption, DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: AdsorbRate(nSpecies)
REAL   , INTENT(OUT)            :: ReactRate(nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: AdsorbActE(nSpecies*Adsorption%ReactNum)
INTEGER, INTENT(OUT)            :: SurfCollNum(nSpecies), AdsorbNum(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iCase, iReact
#ifdef MPI
REAL                            :: AD(nSpecies),RR(nSpecies*Adsorption%ReactNum)
INTEGER                         :: ADN(nSpecies)
#endif /*MPI*/
!===================================================================================================================================

IF(SurfMesh%SurfOnProc)THEN
  IF (DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        AdsorbRate(iSpec) = REAL(Adsorption%AdsorpInfo(iSpec)%NumOfAds) / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        AdsorbRate(iSpec) = 0.
      END IF
    END DO
  ELSE IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        IF (DSMC%WallModel.EQ.1) THEN
          AdsorbRate(iSpec) = Adsorption%AdsorpInfo(iSpec)%MeanProbAds / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
        ELSE IF (DSMC%WallModel.EQ.3) THEN
          AdsorbRate(iSpec) = Adsorption%AdsorpInfo(iSpec)%MeanProbAds / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
        END IF
      ELSE
        AdsorbRate(iSpec)= 0.
      END IF
    END DO
  END IF

  iCase = 1
  DO iSpec = 1, nSpecies
    DO iReact = 1, Adsorption%ReactNum+1
      IF (Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact).GT.0) THEN
        ReactRate(iCase) = Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iReact) &
            / REAL(Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact)) !* REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        ReactRate(iCase) = 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO

  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum
      IF (Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iReact).GT.0) THEN
        IF (DSMC%WallModel.EQ.1) THEN
          AdsorbActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iReact) &
              / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
        ELSE IF (DSMC%WallModel.EQ.3) THEN
          IF (Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact).GT.0) THEN
            AdsorbActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iReact) &
                / REAL(Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact))
          END IF
        END IF
      ELSE
        AdsorbActE(iCase)= 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO
  DO iSpec = 1,nSpecies
    SurfCollNum(iSpec) = Adsorption%AdsorpInfo(iSpec)%WallCollCount
    AdsorbNum(iSpec) = Adsorption%AdsorpInfo(iSpec)%NumOfAds
  END DO
ELSE
  SurfCollNum(:) = 0
  AdsorbRate(:) = 0.
  AdsorbNum(:) = 0
  ReactRate(:) = 0.
  AdsorbActE(:) = 0.
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbRate  ,nSpecies                    ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,SurfCollNum ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbNum   ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,ReactRate   ,nSpecies*(Adsorption%ReactNum+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbActE  ,nSpecies*Adsorption%ReactNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  AdsorbRate = AdsorbRate  / REAL(SurfCOMM%nProcs)
  SurfCollNum= SurfCollNum / REAL(SurfCOMM%nProcs)
  AdsorbNum  = AdsorbNum   / REAL(SurfCOMM%nProcs)
  ReactRate  = ReactRate   / REAL(SurfCOMM%nProcs)
  AdsorbActE = AdsorbActE  / REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(AdsorbRate  ,AD          ,nSpecies                    ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(SurfCollNum ,ADN         ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(AdsorbNum   ,ADN         ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(ReactRate   ,RR          ,nSpecies*(Adsorption%ReactNum+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(AdsorbActE  ,RR          ,nSpecies*Adsorption%ReactNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
    Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
    Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
    Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(:) = 0.
  END DO
END IF

END SUBROUTINE GetAdsRates


SUBROUTINE GetSurfRates(DesorbRate,DesorbNum,ReactRate,SurfaceActE)
!===================================================================================================================================
! Calculate adsorption, desorption and accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: Adsorption, DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: DesorbRate(nSpecies)
INTEGER, INTENT(OUT)            :: DesorbNum(nSpecies)
REAL   , INTENT(OUT)            :: ReactRate(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: SurfaceActE(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iReact, iCase
#ifdef MPI
INTEGER                         :: commSize
REAL                            :: DE(nSpecies)
REAL                            :: RR(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
INTEGER                         :: DEN(nSpecies)
#endif /*MPI*/
!===================================================================================================================================

IF(SurfMesh%SurfOnProc)THEN
  IF (DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount.GT.0) THEN
        DesorbRate(iSpec) = REAL(Adsorption%AdsorpInfo(iSpec)%NumOfDes) / REAL(Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount)
      ELSE
        DesorbRate(iSpec) = 0.
      END IF
    END DO
  ELSE IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    iCase = 1
    DO iSpec = 1,nSpecies
      DesorbRate(iSpec)= 0.
      IF (DSMC%WallModel.EQ.1) THEN
        DO iReact = 1, Adsorption%ReactNum+1
          ReactRate(iCase) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iReact) &
              / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
          iCase = iCase + 1
        END DO
      ELSE IF (DSMC%WallModel.EQ.3) THEN
        DO iReact = 1, Adsorption%ReactNum+1
          IF (Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
            ReactRate(iCase) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iReact) &
                / REAL(Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact))
          ELSE
            ReactRate(iCase) = 0.
          END IF
          iCase = iCase + 1
        END DO
      END IF
    END DO
  END IF
  DO iSpec = 1,nSpecies
    DesorbNum(iSpec) = Adsorption%AdsorpInfo(iSpec)%NumOfDes
  END DO
ELSE
  DesorbNum(:)  = 0
  DesorbRate(:) = 0.
  ReactRate(:)  = 0.
  SurfaceActE(:)= 0.
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbRate  ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbNum   ,nSpecies,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  DesorbRate  = DesorbRate  / REAL(SurfCOMM%nProcs)
  DesorbNum   = DesorbNum   / REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(DesorbRate  ,DE          ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(DesorbNum   ,DEN         ,nSpecies,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
    Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
  END DO
END IF

IF(SurfMesh%SurfOnProc)THEN
  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum+1
      IF (Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
        IF (DSMC%WallModel.EQ.1) THEN
          SurfaceActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iReact) &
              / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
        ELSE IF (DSMC%WallModel.EQ.3) THEN
          SurfaceActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iReact) &
              / REAL(Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact))
        END IF
      ELSE
        SurfaceActE(iCase) = 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO
ELSE
  SurfaceActE(:)= 0.
END IF

#ifdef MPI
commSize = nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE ,ReactRate  ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE ,SurfaceActE,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ReactRate   = ReactRate   / REAL(SurfCOMM%nProcs)
  SurfaceActE = SurfaceActE / REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(ReactRate    ,RR         ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(SurfaceActE  ,RR         ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
    Adsorption%AdsorpInfo(iSPec)%NumOfDes = 0
    Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE = 0.
    Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact = 0.
    Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount = 0.
  END DO
END IF

END SUBROUTINE GetSurfRates


SUBROUTINE GetEvaporationRate(EvaporationRate)
!===================================================================================================================================
! Calculate evaporation rate from number of particles of a species evaporating from surface in the defined analyze time [kg/s]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: Species, nSpecies
USE MOD_Particle_Analyze_Vars
USE MOD_DSMC_Vars             ,ONLY: Liquid
#ifdef MPI
!USE MOD_Particle_Boundary_Vars, ONLY : SurfCOMM
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
USE MOD_TimeDisc_Vars         ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: EvaporationRate(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec
#ifdef MPI
REAL                            :: RD(nSpecies)
#endif /*MPI*/
!===================================================================================================================================
EvaporationRate = 0.

DO iSpec=1,nSpecies
  EvaporationRate(iSpec) = Species(iSpec)%MassIC * Species(iSpec)%MacroParticleFactor &
                        * REAL(Liquid%Info(iSpec)%NumOfDes - Liquid%Info(iSpec)%NumOfAds) / dt
END DO

Liquid%Info(:)%NumOfAds = 0
Liquid%Info(:)%NumOfDes = 0

#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,EvaporationRate,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ELSE
    CALL MPI_REDUCE(EvaporationRate,RD,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
#endif /*MPI*/

END SUBROUTINE GetEvaporationRate
#endif /*(PP_TimeDiscMethod==42)*/
#endif /*(PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)*/

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
USE MOD_PARTICLE_Vars         ,ONLY: PartMPF, usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
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
INTEGER           :: i
REAL(KIND=8)              :: partV2, GammaFac
#ifdef MPI
REAL                      :: RD(nSpecAnalyze)
#endif /*MPI*/
!===================================================================================================================================

Ekin = 0.!d0
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
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000
  !       Ekin = Ekin + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 * PartMPF(i)            
        IF(usevMPF) THEN
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + 0.5 * Species(PartSpecies(i))%MassIC * partV2 * PartMPF(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 * Species(PartSpecies(i))%MassIC * partV2 * PartMPF(i)
        ELSE
          Ekin(nSpecAnalyze) = Ekin(nSpecAnalyze) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
        END IF != usevMPF
      ELSE ! partV2 > 1e6
  !       Ekin = Ekin + (GammaFac - 1) * mass * MPF *c^2
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        IF(usevMPF) THEN
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + PartMPF(i) * (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + PartMPF(i) * (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        ELSE
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + (GammaFac-1.) * Species(PartSpecies(i))%MassIC &
                               * Species(PartSpecies(i))%MacroParticleFactor * c2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + (GammaFac-1.) * Species(PartSpecies(i))%MassIC &
                               * Species(PartSpecies(i))%MacroParticleFactor * c2
        END IF !=usevMPF
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
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000
        IF(usevMPF) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            * PartMPF(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
        END IF ! usevMPF
      ELSE ! partV2 > 1e6
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        IF(usevMPF)THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + PartMPF(i) * (GammaFac-1.) &
                      * Species(PartSpecies(i))%MassIC * c2
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + (GammaFac -1.) &
                      * Species(PartSpecies(i))%MassIC &
                      * Species(PartSpecies(i))%MacroParticleFactor * c2
        END IF ! useuvMPF

      END IF ! par2
    END IF ! particle inside
  END DO ! particleveclength
END IF

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Ekin,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(Ekin  ,RD        ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergy


SUBROUTINE CalcNumPartsOfSpec(NumSpec,SimNumSpec)
!===================================================================================================================================
! computes the number of simulated particles AND number of real particles within the domain
! CAUTION: SimNumSpec equals NumSpec only for constant MPF, not vMPF
! Last section of the routine contains the MPI-communication
! Please not, we do not take the paricle number of the BGGas into account
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PartMPF, usevMPF, PDM,Species,PartSpecies
USE MOD_Particle_Analyze_Vars ,ONLY: CalcNumSpec,nSpecAnalyze
USE MOD_DSMC_Vars             ,ONLY: BGGas
USE MOD_Particle_Mesh_Vars    ,ONLY: Geo
USE MOD_Particle_Vars         ,ONLY: nSpecies
#ifdef MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
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
#ifdef MPI
REAL                               :: RD(nSpecAnalyze)
INTEGER(KIND=8)                    :: ID(nSpecAnalyze)
#endif /*MPI*/
!===================================================================================================================================

IF(usevMPF)THEN ! for MPF differentiate between real particle number and simulated particles
  NumSpec    = 0.
  SimNumSpec = 0
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
        NumSpec(PartSpecies(iPart))    = NumSpec(PartSpecies(iPart)) + PartMPF(iPart)          ! NumSpec = real particle number
        SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1                    ! NumSpec =  particle number
    END IF
  END DO
  IF(BGGas%BGGasSpecies.NE.0) THEN
    !NumSpec(BGGas%BGGasSpecies) = BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor
    !SimNumSpec(BGGas%BGGasSpecies) = INT(NumSpec(BGGas%BGGasSpecies))
    NumSpec(BGGas%BGGasSpecies) = 0.
    SimNumSpec(BGGas%BGGasSpecies) = 0
  END IF
  IF(nSpecAnalyze.GT.1)THEN
    NumSpec(nSpecAnalyze)    = SUM(NumSpec(1:nSpecies))
    SimNumSpec(nSpecAnalyze) = SUM(SimNumSpec(1:nSpecies))
  END IF
ELSE ! no mpf-simulated particle number = real particle number
  SimNumSpec = 0
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1                      ! NumSpec =  particle number
    END IF
  END DO
  NumSpec = REAL(SimNumSpec)
  IF(BGGas%BGGasSpecies.NE.0) THEN
  !  NumSpec(BGGas%BGGasSpecies) = (BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor)
  !  SimNumSpec(BGGas%BGGasSpecies) = 1.
    NumSpec(BGGas%BGGasSpecies) = 0.
    SimNumSpec(BGGas%BGGasSpecies) = 0
  END IF
  IF(nSpecAnalyze.GT.1)THEN
    SimNumSpec(nSpecAnalyze) = SUM(SimNumSpec(1:nSpecies))
    NumSpec(nSpecAnalyze) = SUM(NumSpec(1:nSpecies))
  END IF
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcNumSpec) &
  CALL MPI_REDUCE(MPI_IN_PLACE,SimNumSpec ,nSpecAnalyze,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(NumSpec     ,RD         ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcNumSpec) &
  CALL MPI_REDUCE(SimNumSpec  ,ID         ,nSpecAnalyze,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcNumPartsOfSpec


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
CALL CalcTransTemp(NumSpec, Temp)

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
        TempTotal(iSpec) = Temp(iSpec)
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


SUBROUTINE CalcTransTemp(NumSpec, Temp)
!===================================================================================================================================
! calculate the translational temperature of each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, Species, PDM, nSpecies, BoltzmannConst, PartMPF, usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
#if (PP_TimeDiscMethod==1000)
USE MOD_LD_Vars               ,ONLY: BulkValues
#endif
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: NumSpec(:)    !< global number of REAL particles in domain
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Temp(:)       !< output value is already the GLOBAL temperature
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i, iSpec
REAL              ::  TempDirec(nSpecies,3)
#if (PP_TimeDiscMethod!=1000)
REAL              :: PartVandV2(nSpecies, 6), Mean_PartV2(nSpecies, 3), MeanPartV_2(nSpecies,3)
#endif
#ifdef MPI
REAL              :: RD(nSpecies*6)
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
    IF (usevMPF) THEN
      PartVandV2(PartSpecies(i),1:3) = PartVandV2(PartSpecies(i),1:3) + PartState(i,4:6) * PartMPF(i)
      PartVandV2(PartSpecies(i),4:6) = PartVandV2(PartSpecies(i),4:6) + PartState(i,4:6)**2 * PartMPF(i)
    ELSE
      PartVandV2(PartSpecies(i),1:3) = PartVandV2(PartSpecies(i),1:3) + PartState(i,4:6)
      PartVandV2(PartSpecies(i),4:6) = PartVandV2(PartSpecies(i),4:6) + PartState(i,4:6)**2
    END IF
  END IF
END DO

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PartVandV2,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(PartVandV2  ,RD        ,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
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
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, PDM, nSpecies, BoltzmannConst, PartMPF, usevMPF
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC
USE MOD_DSMC_Analyze          ,ONLY: CalcTVib, CalcTelec, CalcTVibPoly
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
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
    IF (usevMPF) THEN
      EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + PartStateIntEn(iPart,1) * PartMPF(iPart)
      ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(iPart,2) * PartMPF(iPart)
      IF ( DSMC%ElectronicModel .AND. SpecDSMC(PartSpecies(iPart))%InterID .NE. 4) THEN
        Eelec(PartSpecies(iPart)) = Eelec(PartSpecies(iPart)) + PartStateIntEn(iPart,3) * PartMPF(iPart)
      END IF
    ELSE
      EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + PartStateIntEn(iPart,1)
      ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(iPart,2)
      IF ( DSMC%ElectronicModel .AND. SpecDSMC(PartSpecies(iPart))%InterID .NE. 4) THEN
        Eelec(PartSpecies(iPart)) = Eelec(PartSpecies(iPart)) + PartStateIntEn(iPart,3)
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
#ifdef MPI
IF(PartMPI%MPIRoot)THEN
#endif
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
    IF ( DSMC%ElectronicModel ) THEN
      IF ((NumSpecTemp.GT.0).AND.(SpecDSMC(iSpec)%InterID.NE.4) ) THEN
        IntTemp(iSpec,3) = CalcTelec(Eelec(iSpec)/NumSpecTemp,iSpec)
        IntEn(iSpec,3) = Eelec(iSpec)
      ELSE
        IntEn(iSpec,3) = 0.0
      END IF
    END IF
    IntEn(iSpec,1) = EVib(iSpec)
    IntEn(iSpec,2) = ERot(iSpec)
    IF(.NOT.usevMPF) IntEn(iSpec,:) = IntEn(iSpec,:) * Species(iSpec)%MacroParticleFactor
  END DO
  ! Sums of the energy values
  IF(nSpecAnalyze.GT.1) THEN
    IntEn(nSpecAnalyze,1) = SUM(IntEn(:,1))
    IntEn(nSpecAnalyze,2) = SUM(IntEn(:,2))
    IF(DSMC%ElectronicModel) IntEn(nSpecAnalyze,3) = SUM(IntEn(:,3))
  END IF
#ifdef MPI
END IF
#endif

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

SUBROUTINE ReacRates(RRate, NumSpec)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars          ,ONLY: ChemReac, DSMC
USE MOD_TimeDisc_Vars      ,ONLY: dt
USE MOD_Particle_Vars      ,ONLY: Species, nSpecies
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_Particle_MPI_Vars  ,ONLY: PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: RRate(:)
REAL,INTENT(IN)                 :: NumSpec(:) ! is the global number of real particles
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
      CASE('R')
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
      CASE('D')
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
      CASE('E')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN 
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('i')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('r')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%MeshVolume**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
!                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(nSpecies+1))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) * GEO%MeshVolume**2 &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor / (dt * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))     &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2))    &
               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(nSpecies+1))
!               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
        END IF
      CASE('x')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
      RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
           * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('iQK')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
      RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
           * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('rQK')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%MeshVolume**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) * GEO%MeshVolume**2 &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor / (dt * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))     &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2))    &
               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
        END IF
      END SELECT
    END IF
  END DO
END IF
ChemReac%NumReac = 0.

END SUBROUTINE ReacRates
#endif 

#if ( PP_TimeDiscMethod == 42)
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
    IF ( SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
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
    IF(SpecDSMC(iSpec)%InterID.ne.4)THEN
      SpecDSMC(iSpec)%ElectronicTransition = 0
    END IF
  END DO
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ElectronicTransition
#endif

#if ( PP_TimeDiscMethod == 42)
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
    IF (SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
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
        WRITE(iunit,'(A6,A5)',ADVANCE='NO') 'TIME', ' '
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

SUBROUTINE TrackingParticlePosition(time)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PEM
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
USE MOD_Particle_Analyze_Vars ,ONLY: printDiff,printDiffVec,printDiffTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iunit,iPartState
CHARACTER(LEN=60) :: TrackingFilename!,hilf
LOGICAL            :: fexist
REAL               :: diffPos,diffVelo
!===================================================================================================================================

iunit=GETFREEUNIT()
!WRITE(UNIT=hilf,FMT='(I6.6)') MyRank
!TrackingFilename = ('MyRank'//TRIM(hilf)//'_ParticlePosition.csv')
TrackingFilename = ('ParticlePosition.csv')

INQUIRE(FILE = TrackingFilename, EXIST=fexist)
IF(.NOT.fexist) THEN 
 IF(PartMPI%MPIRoot)THEN
   iunit=GETFREEUNIT()
   OPEN(UNIT=iunit,FILE=TrackingFilename,FORM='FORMATTED',STATUS='UNKNOWN')
   !CALL FLUSH (iunit)
    ! writing header
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'TIME', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartNum', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosX', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosY', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosZ', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelX', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelY', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelZ', ' '
    CLOSE(iunit)
  END IF
ELSE
  iunit=GETFREEUNIT()
  OPEN(unit=iunit,FILE=TrackingFileName,FORM='Formatted',POSITION='APPEND',STATUS='old')
  !CALL FLUSH (iunit)
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      WRITE(iunit,104,ADVANCE='NO') TIME
      WRITE(iunit,'(A1)',ADVANCE='NO') ','
      WRITE(iunit,'(I12)',ADVANCE='NO') i
      DO iPartState=1,6
        WRITE(iunit,'(A1)',ADVANCE='NO') ','
        WRITE(iunit,104,ADVANCE='NO') PartState(i,iPartState)
      END DO
      WRITE(iunit,'(A1)',ADVANCE='NO') ','
      WRITE(iunit,'(I12)',ADVANCE='NO') PEM%Element(i)
      WRITE(iunit,'(A)') ' '
     END IF
  END DO
  CLOSE(iunit)
END IF
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
104    FORMAT (e25.14)

END SUBROUTINE TrackingParticlePosition

Function CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars ,ONLY: c2, c2_inv
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species
USE MOD_PARTICLE_Vars ,ONLY: PartMPF, usevMPF
USE MOD_Particle_Vars ,ONLY: PartLorentzType
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
REAL                               :: partV2, gamma1, Ekin
!===================================================================================================================================

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv)
  IF(usevMPF)THEN
    Ekin=PartMPF(iPart)*(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2
  ELSE
    Ekin= (gamma1-1.0)* Species(PartSpecies(iPart))%MassIC*Species(PartSpecies(iPart))%MacroParticleFactor*c2
  END IF
ELSE
  partV2 = PartState(iPart,4) * PartState(iPart,4) &
         + PartState(iPart,5) * PartState(iPart,5) &
         + PartState(iPart,6) * PartState(iPart,6)

  IF(usevMPF)THEN
    IF (partV2.LT.1e6)THEN
      Ekin= 0.5 * Species(PartSpecies(iPart))%MassIC * partV2 * PartMPF(iPart)
    ELSE
      gamma1=partV2*c2_inv
      gamma1=1.0/SQRT(1.-gamma1)
      Ekin=PartMPF(iPart)*(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2
    END IF ! ipartV2
  ELSE ! novMPF
    IF (partV2.LT.1e6)THEN
      Ekin= 0.5*Species(PartSpecies(iPart))%MassIC*partV2* Species(PartSpecies(iPart))%MacroParticleFactor
    ELSE
      gamma1=partV2*c2_inv
      gamma1=1.0/SQRT(1.-gamma1)
      Ekin= (gamma1-1.0)* Species(PartSpecies(iPart))%MassIC*Species(PartSpecies(iPart))%MacroParticleFactor*c2
    END IF ! ipartV2
  END IF ! usevMPF
END IF
CalcEkinPart=Ekin
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
IF(NINT(Species(SpeciesID)%ChargeIC/(-1.60217653E-19)).EQ.1) PartIsElectron=.TRUE.

END FUNCTION PARTISELECTRON


SUBROUTINE CalculateElectronDensityCell() 
!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars     ,ONLY:GEO
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell
USE MOD_Particle_Vars          ,ONLY:Species,PartSpecies,PDM,PEM,usevMPF,PartMPF
USE MOD_Preproc                ,ONLY:PP_nElems
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPart, iElem
!===================================================================================================================================

! nullify
ElectronDensityCell=0.

! loop over all particles and count the number of electrons per cell
! CAUTION: we need the number of all real particle instead of simulated particles
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    IF(.NOT.PARTISELECTRON(iPart)) CYCLE
    IF(usevMPF) THEN
      ElectronDensityCell(PEM%Element(iPart))=ElectronDensityCell(PEM%Element(iPart))+PartMPF(iPart)
    ELSE
      ElectronDensityCell(PEM%Element(iPart))=ElectronDensityCell(PEM%Element(iPart)) &
                                             +Species(PartSpecies(iPart))%MacroParticleFactor
    END IF
  END IF ! ParticleInside
END DO ! iPart

! loop over all elements and divide by volume 
DO iElem=1,PP_nElems
  ElectronDensityCell(iElem)=ElectronDensityCell(iElem)/GEO%Volume(iElem)
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateElectronDensityCell


SUBROUTINE CalculateElectronTemperatureCell() 
!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars     ,ONLY:GEO
USE MOD_Preproc                ,ONLY:PP_nElems
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronTemperatureCell
USE MOD_Particle_Vars          ,ONLY:PDM,PEM,BoltzmannConst
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iPart, iElem, nElectronsPerCell(1:PP_nElems), ElemID
!===================================================================================================================================

! nullify
ElectronTemperatureCell=0.
nElectronsPerCell      =0

! loop over all particles and sum-up the electron energy per cell and count the number of electrons per cell
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    IF(.NOT.PARTISELECTRON(iPart)) CYCLE
    ElemID=PEM%Element(iPart)
    nElectronsPerCell(ElemID)=nElectronsPerCell(ElemID)+1
    ElectronTemperatureCell(ElemID)=ElectronTemperatureCell(ElemID)+CalcEkinPart(iPart)
  END IF ! ParticleInside
END DO ! iPart

! loop over all elements and divide by electrons per cell to get average kinetic energy 
DO iElem=1,PP_nElems
  IF(nElectronsPerCell(iElem).EQ.0) CYCLE
  ElectronTemperatureCell(iElem)=2.*ElectronTemperatureCell(iElem)/(3.*REAL(nElectronsPerCell(iElem))*BoltzmannConst)
END DO ! iElem=1,PP_nElems

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
USE MOD_Globals_Vars           ,ONLY:ElectronCharge,ElectronMass
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
  PlasmaFrequencyCell(iElem) = SQRT((ElectronDensityCell(iElem)*ElectronCharge*ElectronCharge)/(ElectronMass*Eps0))
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
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell,ElectronTemperatureCell,DebyeLengthCell
USE MOD_Globals_Vars           ,ONLY:ElectronCharge
USE MOD_Particle_Vars          ,ONLY:BoltzmannConst
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
  IF(ElectronDensityCell(iElem).LE.0.) CYCLE
  DebyeLengthCell(iElem)=SQRT((eps0*BoltzmannConst*ElectronTemperatureCell(iElem))/(ElectronDensityCell(iElem)*ElectronCharge**2))
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateDebyeLengthCell


SUBROUTINE CalculatePartElemData() 
!===================================================================================================================================
! use the plasma frequency per cell to estimate the pic time step
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Analyze_Vars  ,ONLY:CalcPlasmaFrequency,CalcPICTimeStep,CalcElectronDensity&
                                    ,CalcElectronTemperature,CalcDebyeLength
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

! electron density
IF(CalcElectronDensity) CALL CalculateElectronDensityCell()

! electron temperature
IF(CalcElectronTemperature) CALL CalculateElectronTemperatureCell()

! plasma frequency
IF(CalcPlasmaFrequency) CALL CalculatePlasmaFrequencyCell()

! Debye length 
IF(CalcDebyeLength) CALL CalculateDebyeLengthCell()

! PIC time step
IF(CalcPICTimeStep) CALL CalculatePICTimeStepCell()

END SUBROUTINE CalculatePartElemData


SUBROUTINE FinalizeParticleAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars ,ONLY: ParticleAnalyzeInitIsDone,DebyeLengthCell,PICTimeStepCell &
                                    ,ElectronTemperatureCell,ElectronDensityCell,PlasmaFrequencyCell
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
END SUBROUTINE FinalizeParticleAnalyze
#endif /*PARTICLES*/


END MODULE MOD_Particle_Analyze
