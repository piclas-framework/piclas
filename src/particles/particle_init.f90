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

MODULE MOD_ParticleInit
!===================================================================================================================================
! Add comments please!
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

INTERFACE InitParticleGlobals
  MODULE PROCEDURE InitParticleGlobals
END INTERFACE

INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE FinalizeParticles
  MODULE PROCEDURE FinalizeParticles
END INTERFACE

INTERFACE InitialIonization
  MODULE PROCEDURE InitialIonization
END INTERFACE

INTERFACE InitRandomSeed
  MODULE PROCEDURE InitRandomSeed
END INTERFACE

INTERFACE PortabilityGetPID
  FUNCTION GetPID_C() BIND (C, name='getpid')
    !GETPID() is an intrinstic compiler function in gnu. This routine ensures the portability with other compilers.
    USE ISO_C_BINDING,         ONLY: PID_T => C_INT
    IMPLICIT NONE
    INTEGER(KIND=PID_T)        :: GetPID_C
  END FUNCTION GetPID_C
END INTERFACE

PUBLIC::InitParticleGlobals
PUBLIC::InitParticles
PUBLIC::FinalizeParticles
PUBLIC::DefineParametersParticles
PUBLIC::InitialIonization
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersParticles()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
USE MOD_part_RHS    ,ONLY: DefineParametersParticleRHS
IMPLICIT NONE
!==================================================================================================================================
CALL DefineParametersParticleRHS()

CALL prms%SetSection("Particle")

CALL prms%CreateRealOption(     'Particles-ManualTimeStep'  , 'Manual timestep [sec]. This variable is deprecated. '//&
                                                              'Use ManualTimestep instead.', '-1.0')
CALL prms%CreateRealOption(     'Part-AdaptiveWeightingFactor', 'Weighting factor theta for weighting of average'//&
                                                                ' instantaneous values with those of previous iterations.', '0.001')
CALL prms%CreateIntOption(      'Part-nSpecies' ,                 'Number of species used in calculation', '1')
! Ionization
CALL prms%CreateLogicalOption(  'Part-DoInitialIonization'    , 'When restarting from a state, ionize the species to a '//&
                                                                'specific degree', '.FALSE.')
CALL prms%CreateIntOption(      'InitialIonizationSpecies', 'Supply the number of species that are considered for automatic '//&
                                                            'ionization')
CALL prms%CreateIntArrayOption( 'InitialIonizationSpeciesID', 'Supply a vector with the species IDs that are used for the '//&
                                                              'initial ionization.')
CALL prms%CreateRealOption(     'InitialIonizationChargeAverage' , 'Average charge for each atom/molecule in the cell '//&
                                                                   '(corresponds to the ionization degree)')

CALL prms%CreateIntOption(      'Part-MaxParticleNumber', 'Maximum number of Particles per proc (used for array init)'&
                                                                 , '1')
CALL prms%CreateIntOption(      'Part-NumberOfRandomSeeds'    , 'Number of Seeds for Random Number Generator'//&
                                                                'Choose nRandomSeeds \n'//&
                                                                '=-1    Random \n'//&
                                                                '= 0    Debugging-friendly with hard-coded deterministic numbers\n'//&
                                                                '> 0    Debugging-friendly with numbers from ini. ', '0')
CALL prms%CreateIntOption(      'Particles-RandomSeed[$]'     , 'Seed [$] for Random Number Generator', '1', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'Particles-DoPoissonRounding' , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Flag to perform Poisson sampling'//&
                                                                ' instead of random rounding', '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DoTimeDepInflow'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Insertion and SurfaceFlux with'//&
                                                                ' simple random rounding. Linearly ramping of'//&
                                                                ' inflow-number-of-particles is only possible with'//&
                                                                ' PoissonRounding or DoTimeDepInflow', '.FALSE.')
CALL prms%CreateIntOption(      'Part-RotPeriodicAxi'         , 'Axis of rotational periodicity:'//&
                                                                   'x=1, y=2, z=3', '1')
CALL prms%CreateRealOption(     'Part-RotPeriodicAngle'       , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Angle for rotational periodicity [Grad]'&
                                                              , '1.0')
CALL prms%CreateIntOption(      'Part-nPeriodicVectors'       , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Number of the periodic vectors j=1,...,n.'//&
                                                                   ' Value has to be the same as defined in preprog.ini', '0')
CALL prms%CreateRealArrayOption('Part-PeriodicVector[$]'      , 'TODO-DEFINE-PARAMETER\nVector for periodic boundaries.'//&
                                                                   'Has to be the same as defined in preproc.ini in their'//&
                                                                   ' respective order. ', '1. , 0. , 0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-DelayTime'              , "TODO-DEFINE-PARAMETER\n"//&
                                                                "During delay time the particles,"//&
                                                                    " won't be moved so the EM field can be evolved", '0.0')

CALL prms%CreateRealOption(     'Part-SafetyFactor'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Factor to scale the halo region with MPI'&
                                                              , '1.0')
CALL prms%CreateRealOption(     'Particles-HaloEpsVelo'       , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Halo region radius', '0.')

CALL prms%CreateIntOption(      'NbrOfRegions'                , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Number of regions to be mapped to Elements', '0')
CALL prms%CreateRealArrayOption('RegionBounds[$]'                , 'TODO-DEFINE-PARAMETER\nRegionBounds ((xmin,xmax,ymin,...)'//&
                                                                '|1:NbrOfRegions)'&
                                                                , '0. , 0. , 0. , 0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-RegionElectronRef[$]'   , 'rho_ref, phi_ref, and Te[eV] for Region#'&
                                                              , '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-RegionElectronRef[$]-PhiMax'   , 'max. expected phi for Region#\n'//&
                                                                '(linear approx. above! def.: phi_ref)', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'PrintrandomSeeds'            , 'Flag defining if random seeds are written.', '.FALSE.')
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
CALL prms%CreateLogicalOption(  'velocityOutputAtTime'        , 'Flag if leapfrog uses a velocity-output at real time' , '.FALSE.')
#endif

CALL prms%CreateLogicalOption(  'Part-DoFieldIonization'      , 'Do Field Ionization by quantum tunneling.', '.FALSE.')
CALL prms%CreateIntOption(      'FieldIonizationModel'        , 'Field Ionization models. Implemented models are:\n'//&
                                                                ' * Ammosov-Delone-Krainov (ADK) model Bruhwiler 2003\n'//&
                                                                ' * Ammosov-Delone-Krainov (ADK) model Yu 2018')

CALL prms%SetSection("IMD")
! IMD things
CALL prms%CreateRealOption(     'IMDTimeScale'                , 'Time unit of input file.\n The default value is'//&
                                                                ' ~10.18 fs which comes from the unit system in IMD', '10.18e-15')
CALL prms%CreateRealOption(     'IMDLengthScale'              , 'Length unit scale used by IMD which is 1 angstrom'&
                                                              , '1.0e-10')
CALL prms%CreateStringOption(   'IMDAtomFile'                 , 'IMD data file containing the atomic states for PartState(1:6)'&
                                                              , 'no file found')
CALL prms%CreateStringOption(   'IMDCutOff'                   , 'Atom cut-off parameter for reducing the number of improrted '//&
                                                                'IMD particles\n'//&
                                                                '1.) no_cutoff\n'//&
                                                                '2.) Epot\n'//&
                                                                '3.) coordinates\n'//&
                                                                '4.) velocity', 'no_cutoff')
CALL prms%CreateRealOption(     'IMDCutOffxValue'              ,"Cut-off coordinate for"//&
                                                                " IMDCutOff='coordiantes'" &
                                                              , '-999.9')
CALL prms%CreateIntOption(      'IMDnSpecies'                 , 'Count of IMD species', '1')
CALL prms%CreateStringOption(   'IMDInputFile'                , 'Laser data file name containing '//&
                                                                'PartState(1:6) ' &
                                                              , 'no file found')
CALL prms%SetSection("VMPF")

! vmpf stuff
CALL prms%CreateLogicalOption(  'Part-vMPF'                      , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Flag to use variable '//&
                                                                'Macro Particle Factor.', '.FALSE.')
CALL prms%CreateLogicalOption(  'Part-vMPFPartMerge'              , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Enable Particle Merge routines.'&
                                                              , '.FALSE.')
CALL prms%CreateIntOption(      'Part-vMPFMergePolOrder'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Polynomial degree for vMPF particle merge.'&
                                                              , '2')
CALL prms%CreateIntOption(      'Part-vMPFCellSplitOrder'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Order for cell splitting of variable MPF'&
                                                              , '15')
CALL prms%CreateIntOption(      'Part-vMPFMergeParticleTarget', 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Count of particles wanted after merge.', '0')
CALL prms%CreateIntOption(      'Part-vMPFSplitParticleTarget', 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Number of particles wanted after split.','0')
CALL prms%CreateIntOption(      'Part-vMPFMergeParticleIter'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Number of iterations between particle '//&
                                                                'merges.', '100')
CALL prms%CreateStringOption(   'Part-vMPFvelocityDistribution','TODO-DEFINE-PARAMETER\n'//&
                                                                'Velocity distribution for variable '//&
                                                                'MPF.' , 'OVDR')
CALL prms%CreateLogicalOption(  'Part-vMPFrelativistic'              , 'TODO-DEFINE-PARAMETER', '.FALSE.')


! Output of macroscopic values
CALL prms%SetSection("Particle Sampling")
CALL prms%CreateLogicalOption(  'Part-WriteMacroValues'&
  , 'Set [T] to activate ITERATION DEPENDANT h5 output of macroscopic values sampled every [Part-IterationForMacroVal] iterat'//&
  'ions from particles. Sampling starts from simulation start. Can not be enabled together with Part-TimeFracForSampling.\n'//&
  'If Part-WriteMacroValues is true, WriteMacroVolumeValues and WriteMacroSurfaceValues are forced to be true.\n'//&
  '(HALOWIKI:)Write macro values (e.g. rotational Temperature).'&
  , '.FALSE.')
CALL prms%CreateLogicalOption(  'Part-WriteMacroVolumeValues'&
  , 'Similar to Part-WriteMacroValues. Set [T] to activate iteration dependant sampling and h5 output for each element.'//&
  ' Is automatically set true if Part-WriteMacroValues is true.\n'//&
  'Can not be enabled if Part-TimeFracForSampling is set.', '.FALSE.')
CALL prms%CreateLogicalOption(  'Part-WriteMacroSurfaceValues'&
  , 'Similar to Part-WriteMacroValues. Set [T] to activate iteration dependant sampling and h5 output on surfaces.'//&
  ' Is automatically set true if Part-WriteMacroValues is true.\n'//&
  'Can not be enbaled if Part-TimeFracForSampling is set.', '.FALSE.')
CALL prms%CreateIntOption(      'Part-IterationForMacroVal'&
  , 'Set number of iterations used for sampling if Part-WriteMacroValues is set true.', '1')

CALL prms%CreateRealOption(     'Part-TimeFracForSampling'&
  , 'Set value greater 0.0 to enable TIME DEPENDANT sampling. The given simulation time fraction will be sampled. Sampling'//&
  ' starts after TEnd*(1-Part-TimefracForSampling).\n'//&
  'Can not be enabled together with Part-WriteMacroValues.' , '0.0')
CALL prms%CreateIntOption(      'Particles-NumberForDSMCOutputs'&
  , 'Give the number of outputs for time fraction sampling.\n'//&
  'Default value is 1 if Part-TimeFracForSampling is enabled.', '1')

CALL prms%CreateLogicalOption(  'Particles-DSMC-CalcSurfaceVal'&
  , 'Set [T] to activate sampling, analyze and h5 output for surfaces. Therefore either time fraction or iteration sampling'//&
  ' have to be enabled as well.', '.FALSE.')

END SUBROUTINE DefineParametersParticles


!===================================================================================================================================
! Global particle parameters needed for other particle inits
!===================================================================================================================================
SUBROUTINE InitParticleGlobals()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Tracking_Vars,     ONLY: TrackingMethod,TriaTracking,DoRefMapping
USE MOD_Particle_Vars              ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GLOBALS...'

! Find tracking method immediately, a lot of the later variables depend on it
TrackingMethod = GETINTFROMSTR('TrackingMethod')
SELECT CASE(TrackingMethod)
CASE(REFMAPPING)
  DoRefMapping=.TRUE.
  TriaTracking=.FALSE.
CASE(TRACING)
  DoRefMapping=.FALSE.
  TriaTracking=.FALSE.
CASE(TRIATRACKING)
  DoRefMapping=.FALSE.
  TriaTracking=.TRUE.
CASE DEFAULT
  SWRITE(UNIT_stdOut,'(A)')' TrackingMethod not implemented! Select refmapping (1), tracing (2) or triatracking (3).'
  CALL abort(&
  __STAMP__&
  ,'TrackingMethod not implemented! TrackingMethod=',IntInfoOpt=TrackingMethod)
END SELECT
IF (Symmetry%Order.LE.2) THEN
  DoRefMapping=.FALSE.
  TriaTracking=.TRUE.
  SWRITE(UNIT_stdOut,'(A)') "TrackingMethod set to TriaTracking due to Symmetry2D."
END IF

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GLOBALS DONE'

END SUBROUTINE InitParticleGlobals


SUBROUTINE InitParticles()
!===================================================================================================================================
! Glue Subroutine for particle initialization
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_DSMC_Init                  ,ONLY: InitDSMC
USE MOD_DSMC_Vars                  ,ONLY: useDSMC,DSMC,DSMC_Solution
USE MOD_IO_HDF5                    ,ONLY: AddToElemData,ElementOut
USE MOD_LoadBalance_Vars           ,ONLY: nPartsPerElem
USE MOD_Mesh_Vars                  ,ONLY: nElems
USE MOD_SurfaceModel_Porous        ,ONLY: InitPorousBoundaryCondition
USE MOD_Particle_Boundary_Sampling ,ONLY: InitParticleBoundarySampling
USE MOD_SurfaceModel_Vars          ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound
USE MOD_Particle_Tracking_Vars     ,ONLY: TrackingMethod
USE MOD_Particle_Vars              ,ONLY: ParticlesInitIsDone,WriteMacroVolumeValues,WriteMacroSurfaceValues,nSpecies
USE MOD_Restart_Vars               ,ONLY: DoRestart
USE MOD_Particle_Emission_Init     ,ONLY: InitialParticleInserting
USE MOD_Particle_SurfFlux_Init     ,ONLY: InitializeParticleSurfaceflux
USE MOD_SurfaceModel_Init          ,ONLY: InitSurfaceModel
USE MOD_Particle_Surfaces          ,ONLY: InitParticleSurfaces
USE MOD_Particle_Mesh_Vars         ,ONLY: GEO
USE MOD_Part_Emission              ,ONLY: AdaptiveBCAnalyze
USE MOD_Particle_Boundary_Init     ,ONLY: InitParticleBoundaryRotPeriodic, InitAdaptiveWallTemp
#if USE_MPI
USE MOD_Particle_MPI               ,ONLY: InitParticleCommSize
!USE MOD_Particle_MPI_Emission      ,ONLY: InitEmissionParticlesToProcs
#endif
#if (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Init                ,ONLY: InitFPFlow
USE MOD_Particle_Vars              ,ONLY: Symmetry
#endif
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Init                   ,ONLY: InitBGK
USE MOD_Particle_Vars              ,ONLY: Symmetry
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
IF(ParticlesInitIsDone)THEN
   SWRITE(*,*) "InitParticles already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES ...'

IF(TrackingMethod.NE.TRIATRACKING) THEN
  CALL InitParticleSurfaces()
END IF

IF(.NOT.ALLOCATED(nPartsPerElem))THEN
  ALLOCATE(nPartsPerElem(1:nElems))
  nPartsPerElem=0
  CALL AddToElemData(ElementOut,'nPartsPerElem',LongIntArray=nPartsPerElem(:))
END IF

CALL InitializeVariables()

! Insert the initial particles
CALL InitialParticleInserting()
! Initialize particle surface flux to be performed per iteration
CALL InitializeParticleSurfaceflux()

! Initialize volume sampling
IF(useDSMC .OR. WriteMacroVolumeValues) THEN
! definition of DSMC sampling values
  DSMC%SampNum = 0
  ALLOCATE(DSMC_Solution(1:11,1:nElems,1:nSpecies))
  DSMC_Solution = 0.0
END IF

! Initialize surface sampling / rotational periodic mapping
IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal.OR.(ANY(PartBound%Reactive)).OR.(nPorousBC.GT.0).OR.GEO%RotPeriodicBC) THEN
  CALL InitParticleBoundarySampling()
  CALL InitParticleBoundaryRotPeriodic()
  CALL InitAdaptiveWallTemp()
END IF

! Initialize porous boundary condition (requires BCdata_auxSF and SurfMesh from InitParticleBoundarySampling)
IF(nPorousBC.GT.0) CALL InitPorousBoundaryCondition()

IF (useDSMC) THEN
  CALL InitDSMC()
  CALL InitSurfaceModel()
#if (PP_TimeDiscMethod==300)
IF (Symmetry%Order.EQ.1) CALL abort(__STAMP__&
  ,'ERROR: 1D Fokker-Planck flow is not implemented yet')
  CALL InitFPFlow()
#endif
#if (PP_TimeDiscMethod==400)
IF (Symmetry%Order.EQ.1) CALL abort(__STAMP__&
  ,'ERROR: 1D BGK is not implemented yet')
  CALL InitBGK()
#endif
ELSE IF (WriteMacroVolumeValues.OR.WriteMacroSurfaceValues) THEN
  DSMC%ElectronicModel = 0
END IF

#if USE_MPI
! has to be called AFTER InitializeVariables and InitDSMC
CALL InitParticleCommSize()
#endif

! sampling of near adaptive boundary element values in the first time step to get initial distribution for porous BC
IF(.NOT.DoRestart) THEN
  IF(nPorousBC.GT.0) CALL AdaptiveBCAnalyze(initSampling_opt=.TRUE.)
END IF

ParticlesInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticles


SUBROUTINE InitializeVariables()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_DSMC_Symmetry          ,ONLY: DSMC_1D_InitVolumes, DSMC_2D_InitVolumes, DSMC_2D_InitRadialWeighting
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_Part_RHS               ,ONLY: InitPartRHS
USE MOD_Particle_Mesh          ,ONLY: InitParticleMesh
USE MOD_Particle_Emission_Init ,ONLY: InitializeVariablesSpeciesInits
USE MOD_Particle_Boundary_Init ,ONLY: InitializeVariablesPartBoundary, InitializeVariablesAuxBC
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_Particle_Surfaces_Vars ,ONLY: TriaSurfaceFlux
USE MOD_PICInit                ,ONLY: InitPIC
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: InitEmissionComm
USE MOD_Particle_MPI_Halo      ,ONLY: IdentifyPartExchangeProcs
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolationAnalytic
#endif /*CODE_ANALYZE*/
USE MOD_DSMC_AmbipolarDiffusion,ONLY: InitializeVariablesAmbipolarDiff
USE MOD_TimeDisc_Vars          ,ONLY: ManualTimeStep,useManualTimeStep
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: ManualTimeStepParticle ! temporary variable
!===================================================================================================================================
! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')
CALL AllocateParticleArrays()
CALL InitializeVariablesRandomNumbers()

! initialization of surface model flags
DoPoissonRounding = GETLOGICAL('Particles-DoPoissonRounding','.FALSE.')
DoTimeDepInflow   = GETLOGICAL('Particles-DoTimeDepInflow','.FALSE.')
DelayTime = GETREAL('Part-DelayTime','0.')
!--- Read Manual Time Step: Old variable name still supported
ManualTimeStepParticle = GETREAL('Particles-ManualTimeStep')
IF(ManualTimeStepParticle.GT.0.0)THEN
  ManualTimeStep = ManualTimeStepParticle
  IF (ManualTimeStep.GT.0.0) useManualTimeStep=.True.
END IF ! ManualTimeStepParticle.GT.0.0

nSpecies = GETINT('Part-nSpecies','1')
IF (nSpecies.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF
ALLOCATE(Species(1:nSpecies))

CALL InitializeVariablesSpeciesInits()
! Which Lorentz boost method should be used?
CALL InitPartRHS()
CALL InitializeVariablesPartBoundary()


!-- AuxBCs
CALL InitializeVariablesAuxBC()

!-- Get PIC deposition (skip DSMC, FP-Flow and BGS-Flow related timediscs)
#if (PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400)
DoDeposition    = .FALSE.
!DoInterpolation = .FALSE.
CALL PrintOption('No PIC-related Time discretization, turning deposition off. DoDeposition','*CHANGE',LogOpt=DoDeposition)
!CALL PrintOption('No PIC-related Time discretization, turning interpolation off. DoInterpolation','*CHANGE',LogOpt=DoDeposition)
#else
DoDeposition    = GETLOGICAL('PIC-DoDeposition')
#endif /*(PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400)*/

!-- Get PIC interpolation (could be skipped above, but DSMC octree requires some interpolation variables, which are allocated before
! init DSMC determines whether DSMC%UseOctree is true or false)
DoInterpolation = GETLOGICAL('PIC-DoInterpolation')
#ifdef CODE_ANALYZE
! Check if an analytic function is to be used for interpolation
DoInterpolationAnalytic   = GETLOGICAL('PIC-DoInterpolationAnalytic')
IF(DoInterpolationAnalytic) DoInterpolation = DoInterpolationAnalytic
#endif /*CODE_ANALYZE*/

! Build BGM and initialize particle mesh
CALL InitParticleMesh()
#if USE_MPI
!-- Build MPI communication
CALL IdentifyPartExchangeProcs()
#endif
CALL InitPIC()

! === 2D/1D/Axisymmetric initialization
! Calculate the volumes for 2D simulation (requires the GEO%zminglob/GEO%zmaxglob from InitFIBGM)
IF(Symmetry%Order.EQ.2) CALL DSMC_2D_InitVolumes()
IF(Symmetry%Order.EQ.1) CALL DSMC_1D_InitVolumes()
IF(Symmetry%Axisymmetric) THEN
  IF(RadialWeighting%DoRadialWeighting) THEN
    ! Initialization of RadialWeighting in 2D axisymmetric simulations
    RadialWeighting%PerformCloning = .TRUE.
    CALL DSMC_2D_InitRadialWeighting()
  END IF
  IF(.NOT.TriaTracking) CALL abort(&
    __STAMP__&
    ,'ERROR: Axisymmetric simulation only supported with TriaTracking = T')
  IF(.NOT.TriaSurfaceFlux) CALL abort(&
    __STAMP__&
    ,'ERROR: Axisymmetric simulation only supported with TriaSurfaceFlux = T')
END IF

#if USE_MPI
CALL InitEmissionComm()
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*USE_MPI*/

CALL InitializeVariablesElectronFluidRegions()
CALL InitializeVariablesIMD()
CALL InitializeVariablesWriteMacroValues()
CALL InitializeVariablesvMPF()
CALL InitializeVariablesIonization()
CALL InitializeVariablesVarTimeStep()
CALL InitializeVariablesAmbipolarDiff()

END SUBROUTINE InitializeVariables


SUBROUTINE AllocateParticleArrays()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefMapping
USE MOD_Mesh_Vars               ,ONLY: nElems
USE MOD_DSMC_Vars               ,ONLY: useDSMC
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ALLOCSTAT
!===================================================================================================================================
IF(DoRefMapping)THEN
  ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,' Cannot allocate partposref!')
  PartPosRef=-888.
END IF
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
velocityOutputAtTime = GETLOGICAL('velocityOutputAtTime','.FALSE.')
IF (velocityOutputAtTime) THEN
  ALLOCATE(velocityAtTime(1:3,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) &
    CALL abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate velocityAtTime array!')
  velocityAtTime=0.
END IF
#endif
#if defined(LSERK)
ALLOCATE(Pt_temp(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
Pt_temp=0.
#endif
#if defined(IMPA) || defined(ROS)
CALL InitializeVariablesImplicit()
#endif

ALLOCATE(PartState(1:6,1:PDM%maxParticleNumber)       , &
         LastPartPos(1:3,1:PDM%maxParticleNumber)     , &
         Pt(1:3,1:PDM%maxParticleNumber)              , &
         PartSpecies(1:PDM%maxParticleNumber)         , &
         PDM%ParticleInside(1:PDM%maxParticleNumber)  , &
         PDM%nextFreePosition(1:PDM%maxParticleNumber), &
         PDM%dtFracPush(1:PDM%maxParticleNumber)      , &
         PDM%IsNewPart(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
PDM%ParticleInside(1:PDM%maxParticleNumber) = .FALSE.
PDM%dtFracPush(1:PDM%maxParticleNumber)     = .FALSE.
PDM%IsNewPart(1:PDM%maxParticleNumber)      = .FALSE.
LastPartPos(1:3,1:PDM%maxParticleNumber)    = 0.
PartState=0.
Pt=0.
PartSpecies        = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)=0

ALLOCATE(PEM%GlobalElemID(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
  ,' Cannot allocate PEM%GlobalElemID(1:PDM%maxParticleNumber) array!')

ALLOCATE(PEM%LastGlobalElemID(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
  ,' Cannot allocate PEM%LastGlobalElemID(1:PDM%maxParticleNumber) array!')

IF (useDSMC) THEN
  ALLOCATE(PEM%pStart(1:nElems)                         , &
           PEM%pNumber(1:nElems)                        , &
           PEM%pEnd(1:nElems)                           , &
           PEM%pNext(1:PDM%maxParticleNumber)           , STAT=ALLOCSTAT)
           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
__STAMP__&
    , ' Cannot allocate DSMC PEM arrays!')
  END IF
END IF
IF (useDSMC) THEN
  ALLOCATE(PDM%PartInit(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
__STAMP__&
    ,' Cannot allocate DSMC PEM arrays!')
  END IF
END IF

END SUBROUTINE AllocateParticleArrays

SUBROUTINE InitializeVariablesIonization()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Initial Ionization of species
DoInitialIonization = GETLOGICAL('Part-DoInitialIonization','.FALSE.')
IF(DoInitialIonization)THEN
  ! Supply the number of species that are considered for automatic ionization
  InitialIonizationSpecies = GETINT('InitialIonizationSpecies')
  ALLOCATE(InitialIonizationSpeciesID(1:InitialIonizationSpecies))

  ! Supply a vector with the species IDs
  InitialIonizationSpeciesID = GETINTARRAY('InitialIonizationSpeciesID',InitialIonizationSpecies)

  ! Average charge for each atom/molecule in the cell (corresponds to the ionization degree)
  InitialIonizationChargeAverage = GETREAL('InitialIonizationChargeAverage')
END IF
DoFieldIonization = GETLOGICAL('Part-DoFieldIonization')
IF(DoFieldIonization)THEN
  FieldIonizationModel = GETINT('FieldIonizationModel')
END IF

END SUBROUTINE InitializeVariablesIonization


SUBROUTINE InitializeVariablesElectronFluidRegions()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Tools ,ONLY: MapRegionToElem
USE MOD_Particle_Mesh_Vars  ,ONLY: NbrOfRegions,RegionBounds
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf, hilf2
INTEGER               :: iRegions
REAL                  :: phimax_tmp
!===================================================================================================================================
!-- Read parameters for region mapping
NbrOfRegions = GETINT('NbrOfRegions','0')
IF (NbrOfRegions .GT. 0) THEN
  ALLOCATE(RegionBounds(1:6,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I0)') iRegions
    RegionBounds(1:6,iRegions) = GETREALARRAY('RegionBounds'//TRIM(hilf2),6,'0. , 0. , 0. , 0. , 0. , 0.')
  END DO

  CALL MapRegionToElem()
  ALLOCATE(RegionElectronRef(1:3,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I0)') iRegions
    ! 1:3 - rho_ref, phi_ref, and Te[eV]
    RegionElectronRef(1:3,iRegions) = GETREALARRAY('Part-RegionElectronRef'//TRIM(hilf2),3,'0. , 0. , 1.')
    WRITE(UNIT=hilf,FMT='(G0)') RegionElectronRef(2,iRegions)
    phimax_tmp = GETREAL('Part-RegionElectronRef'//TRIM(hilf2)//'-PhiMax',TRIM(hilf))
    IF (phimax_tmp.NE.RegionElectronRef(2,iRegions)) THEN !shift reference point (rho_ref, phi_ref) to phi_max:
      RegionElectronRef(1,iRegions) = RegionElectronRef(1,iRegions) &
        * EXP((phimax_tmp-RegionElectronRef(2,iRegions))/RegionElectronRef(3,iRegions))
      RegionElectronRef(2,iRegions) = phimax_tmp
      SWRITE(*,*) 'WARNING: BR-reference point is shifted to:', RegionElectronRef(1:2,iRegions)
    END IF
  END DO
END IF

END SUBROUTINE InitializeVariablesElectronFluidRegions

SUBROUTINE InitializeVariablesVarTimeStep()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Mesh_Vars               ,ONLY: nElems
USE MOD_Particle_Tracking_Vars  ,ONLY: TriaTracking
USE MOD_Particle_VarTimeStep    ,ONLY: VarTimeStep_CalcElemFacs
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! ------- Variable Time Step Initialization (parts requiring completed particle_init and readMesh)
IF(VarTimeStep%UseVariableTimeStep) THEN
  ! Initializing the particle time step array used during calculation for the distribution (after maxParticleNumber was read-in)
  ALLOCATE(VarTimeStep%ParticleTimeStep(1:PDM%maxParticleNumber))
  VarTimeStep%ParticleTimeStep = 1.
  IF(.NOT.TriaTracking) THEN
    CALL abort(&
      __STAMP__&
      ,'ERROR: Variable time step is only supported with TriaTracking = T')
  END IF
  IF(VarTimeStep%UseLinearScaling) THEN
    IF(Symmetry%Order.LE.2) THEN
      ! 2D: particle-wise scaling in the radial direction, ElemFac array only utilized for the output of the time step
      ALLOCATE(VarTimeStep%ElemFac(nElems))
      VarTimeStep%ElemFac = 1.0
    ELSE
      ! 3D: The time step for each cell is precomputed, ElemFac is allocated in the routine
      CALL VarTimeStep_CalcElemFacs()
    END IF
  END IF
  IF(VarTimeStep%UseDistribution) THEN
    ! ! Apply a min-mean filter combo if the distribution was adapted
    ! ! (is performed here to have the element neighbours already defined)
    ! IF(VarTimeStep%AdaptDistribution) CALL VarTimeStep_SmoothDistribution()
    ! Disable AdaptDistribution to avoid adapting during a load balance restart
    IF(VarTimeStep%AdaptDistribution) VarTimeStep%AdaptDistribution = .FALSE.
  END IF
END IF

END SUBROUTINE InitializeVariablesVarTimeStep


SUBROUTINE InitializeVariablesRandomNumbers()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSeed, nRandomSeeds, SeedSize
CHARACTER(32)         :: hilf
!===================================================================================================================================
!--- initialize randomization
! Read print flags
printRandomSeeds = GETLOGICAL('printRandomSeeds','.FALSE.')
nRandomSeeds = GETINT('Part-NumberOfRandomSeeds','0')
CALL RANDOM_SEED(Size = SeedSize)    ! specifies compiler specific minimum number of seeds
ALLOCATE(Seeds(SeedSize))
Seeds(:)=1 ! to ensure a solid run when an unfitting number of seeds is provided in ini
IF(nRandomSeeds.EQ.-1) THEN
  ! ensures different random numbers through irreproducable random seeds (via System_clock)
  CALL InitRandomSeed(nRandomSeeds,SeedSize,Seeds)
ELSE IF(nRandomSeeds.EQ.0) THEN
 !   IF (Restart) THEN
 !   CALL !numbers from state file
 ! ELSE IF (.NOT.Restart) THEN
  CALL InitRandomSeed(nRandomSeeds,SeedSize,Seeds)
ELSE IF(nRandomSeeds.GT.0) THEN
  ! read in numbers from ini
  IF(nRandomSeeds.GT.SeedSize) THEN
    SWRITE (*,*) 'Expected ',SeedSize,'seeds. Provided ',nRandomSeeds,'. Computer uses default value for all unset values.'
  ELSE IF(nRandomSeeds.LT.SeedSize) THEN
    SWRITE (*,*) 'Expected ',SeedSize,'seeds. Provided ',nRandomSeeds,'. Computer uses default value for all unset values.'
  END IF
  DO iSeed=1,MIN(SeedSize,nRandomSeeds)
    WRITE(UNIT=hilf,FMT='(I0)') iSeed
    Seeds(iSeed)= GETINT('Particles-RandomSeed'//TRIM(hilf))
  END DO
  IF (ALL(Seeds(:).EQ.0)) THEN
    CALL ABORT(&
     __STAMP__&
     ,'Not all seeds can be set to zero ')
  END IF
  CALL InitRandomSeed(nRandomSeeds,SeedSize,Seeds)
ELSE
  SWRITE (*,*) 'Error: nRandomSeeds not defined.'//&
  'Choose nRandomSeeds'//&
  '=-1    pseudo random'//&
  '= 0    hard-coded deterministic numbers'//&
  '> 0    numbers from ini. Expected ',SeedSize,'seeds.'
END IF

END SUBROUTINE InitializeVariablesRandomNumbers


SUBROUTINE InitializeVariablesWriteMacroValues()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_TimeDisc_Vars          ,ONLY: TEnd
USE MOD_Particle_Boundary_Vars ,ONLY: AdaptWallTemp
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Include surface values in the macroscopic output
DSMC%CalcSurfaceVal = GETLOGICAL('Particles-DSMC-CalcSurfaceVal','.FALSE.')
! Sampling for and output every given number of iterations (sample is reset after an output)
WriteMacroValues = GETLOGICAL('Part-WriteMacroValues','.FALSE.')
IF(WriteMacroValues)THEN
  WriteMacroVolumeValues = GETLOGICAL('Part-WriteMacroVolumeValues','.TRUE.')
  WriteMacroSurfaceValues = GETLOGICAL('Part-WriteMacroSurfaceValues','.TRUE.')
  IF(.NOT.(WriteMacroVolumeValues.AND.WriteMacroSurfaceValues))THEN
     CALL abort(&
__STAMP__&
    ,'ERROR in particle_init.f90: Part-WriteMacroValues=T => WriteMacroVolumeValues and WriteMacroSurfaceValues must be T!')
  END IF
ELSE
  WriteMacroVolumeValues = GETLOGICAL('Part-WriteMacroVolumeValues','.FALSE.')
  WriteMacroSurfaceValues = GETLOGICAL('Part-WriteMacroSurfaceValues','.FALSE.')
  IF(WriteMacroVolumeValues.AND.WriteMacroSurfaceValues)THEN
    WriteMacroValues = .TRUE.
  END IF
END IF
MacroValSamplIterNum = GETINT('Part-IterationForMacroVal','1')
IF(WriteMacroSurfaceValues.AND.(.NOT.DSMC%CalcSurfaceVal)) DSMC%CalcSurfaceVal = .TRUE.

! Continuous sampling with multiple outputs (sample is not reset after an output)
DSMC%TimeFracSamp = GETREAL('Part-TimeFracForSampling')
IF(DSMC%TimeFracSamp.GT.0.0) THEN
  IF(WriteMacroVolumeValues.OR.WriteMacroSurfaceValues)THEN
    CALL abort(__STAMP__, &
      'ERROR Init Macrosampling: WriteMacroValues and Time fraction sampling enabled at the same time')
  END IF
  DSMC%NumOutput = GETINT('Particles-NumberForDSMCOutputs')
  IF(DSMC%NumOutput.NE.0) THEN
    DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
  ELSE
    CALL abort(__STAMP__,&
      'ERROR Init Macrosampling: Please define a number of outputs (Particles-NumberForDSMCOutputs) greater than zero!')
  END IF
ELSE
  IF(DSMC%NumOutput.NE.0) THEN
    SWRITE(UNIT_STDOUT,*)'WARNING: NumberForDSMCOutputs was set to 0 because TimeFracForSampling is 0.0'
  END IF
  DSMC%NumOutput = 0
END IF

! Adaptive wall temperature should not be used with continuous sampling with multiple outputs as the sample is not reset
IF(AdaptWallTemp) THEN
  IF (DSMC%NumOutput.GT.1) THEN
    CALL abort(__STAMP__, &
      'ERROR: Enabled adaptation of the wall temperature and multiple outputs during a continuous sample is not supported!')
  END IF
END IF

END SUBROUTINE InitializeVariablesWriteMacroValues

SUBROUTINE InitializeVariablesvMPF()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Part_MPFtools          ,ONLY: DefinePolyVec, DefineSplitVec
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ALLOCSTAT
!===================================================================================================================================
! init varibale MPF per particle
IF (usevMPF) THEN
  enableParticleMerge = GETLOGICAL('Part-vMPFPartMerge','.FALSE.')
  IF (enableParticleMerge) THEN
    vMPFMergePolyOrder = GETINT('Part-vMPFMergePolOrder','2')
    vMPFMergeCellSplitOrder = GETINT('Part-vMPFCellSplitOrder','15')
    vMPFMergeParticleTarget = GETINT('Part-vMPFMergeParticleTarget','0')
    IF (vMPFMergeParticleTarget.EQ.0) WRITE(*,*) 'vMPFMergeParticleTarget equals zero: no merging is performed!'
    vMPFSplitParticleTarget = GETINT('Part-vMPFSplitParticleTarget','0')
    IF (vMPFSplitParticleTarget.EQ.0) WRITE(*,*) 'vMPFSplitParticleTarget equals zero: no split is performed!'
    vMPFMergeParticleIter = GETINT('Part-vMPFMergeParticleIter','100')
    vMPF_velocityDistribution = TRIM(GETSTR('Part-vMPFvelocityDistribution','OVDR'))
    vMPF_relativistic = GETLOGICAL('Part-vMPFrelativistic','.FALSE.')
    IF(vMPF_relativistic.AND.(vMPF_velocityDistribution.EQ.'MBDR')) THEN
      CALL abort(&
__STAMP__&
      ,'Relativistic handling of vMPF is not possible using MBDR velocity distribution!')
    END IF
    ALLOCATE(vMPF_SpecNumElem(1:nElems,1:nSpecies))
    CALL DefinePolyVec(vMPFMergePolyOrder)
    CALL DefineSplitVec(vMPFMergeCellSplitOrder)
  END IF
  ALLOCATE(PartMPF(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
__STAMP__&
    ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
  END IF
END IF
END SUBROUTINE InitializeVariablesvMPF


SUBROUTINE InitializeVariablesIMD()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: ElementaryCharge
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefMapping
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iIMDSpec, iSpec
LOGICAL               :: IsIMDSpecies
CHARACTER(32)         :: hilf
!===================================================================================================================================
! IMD data import from *.chkpt file
DoImportIMDFile=.FALSE. ! default
IMDLengthScale=0.0

IMDTimeScale          = GETREAL('IMDTimeScale','10.18e-15')
IMDLengthScale        = GETREAL('IMDLengthScale','1.0E-10')
IMDAtomFile           = GETSTR( 'IMDAtomFile','no file found')
IMDCutOff             = GETSTR( 'IMDCutOff','no_cutoff')
IMDCutOffxValue       = GETREAL('IMDCutOffxValue','-999.9')

IF(TRIM(IMDAtomFile).NE.'no file found')DoImportIMDFile=.TRUE.
IF(DoImportIMDFile)THEN
  DoRefMapping=.FALSE. ! for faster init don't use DoRefMapping!
  CALL PrintOption('DoImportIMDFile=T. Setting DoRefMapping =','*CHANGE',LogOpt=DoRefMapping)
END IF

! get information for IMD atom/ion charge determination and distribution
IMDnSpecies         = GETINT('IMDnSpecies','1')
IMDInputFile        = GETSTR('IMDInputFile','no file found')
ALLOCATE(IMDSpeciesID(IMDnSpecies))
ALLOCATE(IMDSpeciesCharge(IMDnSpecies))
iIMDSpec=1
DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  IsIMDSpecies = GETLOGICAL('Part-Species'//TRIM(hilf)//'-IsIMDSpecies','.FALSE.')
  IF(IsIMDSpecies)THEN
    IMDSpeciesID(iIMDSpec)=iSpec
    IMDSpeciesCharge(iIMDSpec)=NINT(Species(iSpec)%ChargeIC/ElementaryCharge)
    iIMDSpec=iIMDSpec+1
  END IF
END DO
END SUBROUTINE InitializeVariablesIMD


#if defined(IMPA) || defined(ROS)
SUBROUTINE InitializeVariablesImplicit()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_TimeDisc_Vars ,ONLY: nRKStages
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ALLOCSTAT
!===================================================================================================================================
#ifdef IMPA
ALLOCATE(PartStage(1:6,1:nRKStages-1,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartStage arrays!')
END IF
ALLOCATE(PartStateN(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartStateN arrays!')
END IF
ALLOCATE(PartQ(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate PartQ arrays!')
END IF
! particle function values at X0
ALLOCATE(F_PartX0(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate F_PartX0 arrays!')
END IF
! particle function values at Xk
ALLOCATE(F_PartXk(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate F_PartXk arrays!')
END IF
! and the required norms
ALLOCATE(Norm_F_PartX0(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate Norm_F_PartX0 arrays!')
END IF
ALLOCATE(Norm_F_PartXk(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate Norm_F_PartXk arrays!')
END IF
ALLOCATE(Norm_F_PartXk_old(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate Norm_F_PartXk_old arrays!')
END IF
ALLOCATE(PartDeltaX(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate PartDeltaX arrays!')
END IF
ALLOCATE(PartLambdaAccept(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate PartLambdaAccept arrays!')
END IF
ALLOCATE(DoPartInNewton(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
      ,'Cannot allocate DoPartInNewton arrays!')
END IF
#endif /* IMPA */
#ifdef ROS
ALLOCATE(PartStage(1:6,1:nRKStages-1,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartStage arrays!')
END IF
ALLOCATE(PartStateN(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartStateN arrays!')
END IF
ALLOCATE(PartQ(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate PartQ arrays!')
END IF
ALLOCATE(PartDtFrac(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartDtFrac arrays!')
END IF
PartDtFrac=1.
ALLOCATE(PEM%NormVec(1:3,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
   CALL abort(&
 __STAMP__&
   ,' Cannot allocate the normal vector for reflections!')
END IF
PEM%NormVec=0
ALLOCATE(PEM%PeriodicMoved(1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
   CALL abort(&
 __STAMP__&
   ,' Cannot allocate the stage position and element arrays!')
END IF
PEM%PeriodicMoved=.FALSE.
#endif /* ROSENBROCK */

#if IMPA
ALLOCATE(PartIsImplicit(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartIsImplicit arrays!')
END IF
PartIsImplicit=.FALSE.
ALLOCATE(PartDtFrac(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartDtFrac arrays!')
END IF
PartDtFrac=1.
ALLOCATE(PEM%NormVec(1:3,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
   CALL abort(&
 __STAMP__&
   ,' Cannot allocate the normal vector for reflections!')
END IF
PEM%NormVec=0
ALLOCATE(PEM%PeriodicMoved(1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
   CALL abort(&
 __STAMP__&
   ,' Cannot allocate the stage position and element arrays!')
END IF
PEM%PeriodicMoved=.FALSE.
#endif

END SUBROUTINE InitializeVariablesImplicit
#endif

SUBROUTINE InitialIonization()
!----------------------------------------------------------------------------------------------------------------------------------!
! 1.) assign charges to each atom/molecule using the charge supplied by the user
! 2.) reconstruct the electron phase space using the summed charges in each cell for which an electron is
!     created to achieve an ionization degree supplied by the user
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY:  ElementaryCharge
USE MOD_PreProc
USE MOD_Particle_Vars       ,ONLY: PDM,PEM,PartState,nSpecies,Species,PartSpecies
USE MOD_Particle_Vars       ,ONLY: InitialIonizationChargeAverage,InitialIonizationSpeciesID,InitialIonizationSpecies
USE MOD_Mesh_Vars           ,ONLY: NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo,offsetElem
USE MOD_DSMC_Vars           ,ONLY: CollisMode,DSMC,PartStateIntEn
USE MOD_part_emission_tools ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_Vars           ,ONLY: useDSMC
USE MOD_Eval_xyz            ,ONLY: TensorProductInterpolation
USE MOD_Particle_Analyze    ,ONLY: PARTISELECTRON
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: ChargeLower,ChargeUpper
INTEGER       :: ElemCharge(1:PP_nElems)
INTEGER       :: ElecSpecIndx,iSpec,location,iElem,iPart,ParticleIndexNbr
REAL          :: ChargeProbability
REAL          :: iRan
REAL          :: PartPosRef(1:3)
REAL          :: CellElectronTemperature=300
!INTEGER       :: SpeciesID(3)
!INTEGER       :: SpeciesCharge(3)
CHARACTER(32) :: hilf
INTEGER       :: SpeciesCharge(1:InitialIonizationSpecies)
INTEGER       :: I
!===================================================================================================================================
! Determine the charge number of each species
DO I = 1, InitialIonizationSpecies
  iSpec = InitialIonizationSpeciesID(I)
  SpeciesCharge(I) = NINT(Species(iSpec)%ChargeIC/(ElementaryCharge))
END DO ! I = 1, InitialIonizationSpecies

! ---------------------------------------------------------------------------------------------------------------------------------
! 1.) reconstruct ions and determine charge
! ---------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,*)'InitialIonization:'
SWRITE(UNIT_stdOut,*)'  1.) Reconstructing ions and determining the charge of each particle'

! Initialize the element charge with zero
ElemCharge(1:PP_nElems)=0

! Loop over all particles in the vector list
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN

    ! If the current particle is part of the species list, then a charge can be assigned
    IF(ANY(PartSpecies(iPart).EQ.InitialIonizationSpeciesID(:)))THEN !

      IF(PARTISELECTRON(iPart))THEN
        CYCLE
      END IF
      ! Get the cell charge average value and select and upper and lower charge number
      ChargeLower       = INT(InitialIonizationChargeAverage)
      ChargeUpper       = ChargeLower+1
      ChargeProbability = REAL(ChargeUpper)-InitialIonizationChargeAverage !e.g. 2-1.4 = 0.6 -> 60% probability to get lower charge

      ! Compare the random number with the charge difference
      CALL RANDOM_NUMBER(iRan)
      IF(iRan.LT.ChargeProbability)THEN ! Select the lower charge number
        ! Determines the location of the element in the array with min value: get the index of the corresponding charged ion
        ! species
        location                            = MINLOC(ABS(SpeciesCharge-ChargeLower),1)
        ElemCharge(PEM%LocalElemID(iPart))      = ElemCharge(PEM%LocalElemID(iPart))+ChargeLower
      ELSE ! Select the upper charge number
        ! Determines the location of the element in the array with min value: get the index of the corresponding charged ion
        ! species
        location                            = MINLOC(ABS(SpeciesCharge-ChargeUpper),1)
        ElemCharge(PEM%LocalElemID(iPart))      = ElemCharge(PEM%LocalElemID(iPart))+ChargeUpper
      END IF

      ! Set the species ID to atom/singly charged ion/doubly charged ... and so on
      PartSpecies(iPart)=InitialIonizationSpeciesID(location)
    END IF
  END IF
END DO

! ---------------------------------------------------------------------------------------------------------------------------------
! 2.) reconstruct electrons
! ---------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,*)'  2.) Reconstructing electrons'

! Initialize the species index for the electron species with -1
ElecSpecIndx = -1

! Loop over all species and find the index corresponding to the electron species: take the first electron species that is
! encountered
DO iSpec = 1, nSpecies
  IF (Species(iSpec)%ChargeIC.GE.0.0) CYCLE
    IF(NINT(Species(iSpec)%ChargeIC/(-ElementaryCharge)).EQ.1)THEN
      ElecSpecIndx = iSpec
    EXIT
  END IF
END DO
IF (ElecSpecIndx.EQ.-1) CALL abort(&
  __STAMP__&
  ,'Electron species not found. Cannot create electrons without the defined species!')

WRITE(UNIT=hilf,FMT='(I0)') iSpec
SWRITE(UNIT_stdOut,'(A)')'  Using iSpec='//TRIM(hilf)//' as electron species index.'
WRITE(UNIT=hilf,FMT=WRITEFORMAT) CellElectronTemperature
SWRITE(UNIT_stdOut,'(A)')'  Using T='//TRIM(hilf)//' K for the initial electron temperatture (maxwell_lpn) in each cell.'
! Loop over all elements and the sum of charges in each element (for each charge assigned in an element, an electron is created)
DO iElem=1,PP_nElems
  DO iPart=1,ElemCharge(iElem) ! 1 electron for each charge of each element

    ! Set the next free position in the particle vector list
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
    ParticleIndexNbr            = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    PDM%ParticleVecLength       = PDM%ParticleVecLength + 1

    !Set new SpeciesID of new particle (electron)
    PDM%ParticleInside(ParticleIndexNbr) = .true.
    PartSpecies(ParticleIndexNbr) = ElecSpecIndx

    ! Place the electron randomly in the reference cell
    CALL RANDOM_NUMBER(PartPosRef(1:3)) ! get random reference space
    PartPosRef(1:3)=PartPosRef(1:3)*2. - 1. ! map (0,1) -> (-1,1)

    ! Get the physical coordinates that correspond to the reference coordinates
    CALL TensorProductInterpolation(PartPosRef(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem) &
                      ,PartState(1:3,ParticleIndexNbr)) !Map into phys. space

    ! Set the internal energies (vb, rot and electronic) to zero if needed
    IF ((useDSMC).AND.(CollisMode.GT.1)) THEN
      PartStateIntEn(1,ParticleIndexNbr) = 0.
      PartStateIntEn(2,ParticleIndexNbr) = 0.
      IF (DSMC%ElectronicModel.GT.0)  PartStateIntEn(3,ParticleIndexNbr) = 0.
    END IF

    ! Set the element ID of the electron to the current element ID
    PEM%GlobalElemID(ParticleIndexNbr) = iElem + offsetElem

    ! Set the electron velocity using the Maxwellian distribution (use the function that is suitable for small numbers)
    CALL CalcVelocity_maxwell_lpn(ElecSpecIndx, PartState(4:6,ParticleIndexNbr),&
                                  Temperature=CellElectronTemperature)
  END DO
END DO


END SUBROUTINE InitialIonization


SUBROUTINE FinalizeParticles()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize particle variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Halo      ,ONLY: FinalizePartExchangeProcs
USE MOD_PICDepo_Vars           ,ONLY: SendShapeElemID,SendElemShapeID,ShapeMapping,CNShapeMapping
#endif /*USE_MPI*/
!#if USE_MPI
!USE MOD_Particle_MPI_Emission      ,ONLY: FinalizeEmissionParticlesToProcs
!#endif
!USE MOD_DSMC_Vars,                  ONLY: SampDSMC
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!#if USE_MPI
!CALL FinalizeEmissionParticlesToProcs()
!#endif
#if defined(LSERK)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
SDEALLOCATE( Pt_temp)
#endif
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
IF (velocityOutputAtTime) THEN
  SDEALLOCATE(velocityAtTime)
END IF
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
#if defined(ROS) || defined(IMPA)
SDEALLOCATE(PartStage)
SDEALLOCATE(PartStateN)
SDEALLOCATE(PartQ)
SDEALLOCATE(PartDtFrac)
SDEALLOCATE(PEM%GlobalElemID)
SDEALLOCATE(PEM%NormVec)
SDEALLOCATE(PEM%PeriodicMoved)
#endif /*defined(ROS) || defined(IMPA)*/
#if defined(IMPA)
SDEALLOCATE(F_PartXk)
SDEALLOCATE(F_PartX0)
SDEALLOCATE(Norm_F_PartXk_old)
SDEALLOCATE(Norm_F_PartXk)
SDEALLOCATE(Norm_F_PartX0)
SDEALLOCATE(PartDeltaX)
SDEALLOCATE(PartLambdaAccept)
SDEALLOCATE(DoPartInNewton)
SDEALLOCATE(PartIsImplicit)
#endif /*defined(IMPA)*/
SDEALLOCATE(PartPosRef)
SDEALLOCATE(PartState)
SDEALLOCATE(LastPartPos)
SDEALLOCATE(PartSpecies)
SDEALLOCATE(Pt)
SDEALLOCATE(PDM%ParticleInside)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%dtFracPush)
SDEALLOCATE(PDM%IsNewPart)
SDEALLOCATE(vMPF_SpecNumElem)
SDEALLOCATE(PartMPF)
SDEALLOCATE(Species)
SDEALLOCATE(Adaptive_MacroVal)
SDEALLOCATE(SpecReset)
SDEALLOCATE(IMDSpeciesID)
SDEALLOCATE(IMDSpeciesCharge)
SDEALLOCATE(VarTimeStep%ParticleTimeStep)
SDEALLOCATE(VarTimeStep%ElemFac)
SDEALLOCATE(VarTimeStep%ElemWeight)
SDEALLOCATE(PEM%GlobalElemID)
SDEALLOCATE(PEM%LastGlobalElemID)
SDEALLOCATE(PEM%pStart)
SDEALLOCATE(PEM%pNumber)
SDEALLOCATE(PEM%pEnd)
SDEALLOCATE(PEM%pNext)
SDEALLOCATE(seeds)
SDEALLOCATE(RegionBounds)
SDEALLOCATE(RegionElectronRef)
SDEALLOCATE(PartPosLandmark)
#if USE_MPI
SDEALLOCATE(SendShapeElemID)
SDEALLOCATE(SendElemShapeID)
SDEALLOCATE(ShapeMapping)
SDEALLOCATE(CNShapeMapping)
! particle MPI halo exchange
CALL FinalizePartExchangeProcs()
#endif
END SUBROUTINE FinalizeParticles


SUBROUTINE InitRandomSeed(nRandomSeeds,SeedSize,Seeds)
!===================================================================================================================================
!> Initialize pseudo random numbers: Create Random_seed array
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
!===================================================================================================================================
IMPLICIT NONE
! VARIABLES
INTEGER,INTENT(IN)             :: nRandomSeeds
INTEGER,INTENT(IN)             :: SeedSize
INTEGER,INTENT(INOUT)          :: Seeds(SeedSize)
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                        :: iSeed,DateTime(8),ProcessID,iStat,OpenFileID,GoodSeeds
INTEGER(KIND=8)                :: Clock,AuxilaryClock
LOGICAL                        :: uRandomExists
!==================================================================================================================================

uRandomExists=.FALSE.
IF (nRandomSeeds.NE.-1) THEN
  Clock     = 1536679165842_8
  ProcessID = 3671
ELSE
! First try if the OS provides a random number generator
  OPEN(NEWUNIT=OpenFileID, FILE="/dev/urandom", ACCESS="stream", &
       FORM="unformatted", ACTION="read", STATUS="old", IOSTAT=iStat)
  IF (iStat.EQ.0) THEN
    READ(OpenFileID) Seeds
    CLOSE(OpenFileID)
    uRandomExists=.TRUE.
  ELSE
    ! Fallback to XOR:ing the current time and pid. The PID is
    ! useful in case one launches multiple instances of the same
    ! program in parallel.
    CALL SYSTEM_CLOCK(COUNT=Clock)
    IF (Clock .EQ. 0) THEN
      CALL DATE_AND_TIME(values=DateTime)
      Clock =(DateTime(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
      + DateTime(2) * 31_8 * 24 * 60 * 60 * 1000 &
      + DateTime(3) * 24_8 * 60 * 60 * 1000 &
      + DateTime(5) * 60 * 60 * 1000 &
      + DateTime(6) * 60 * 1000 &
      + DateTime(7) * 1000 &
      + DateTime(8)
    END IF
    ProcessID = GetPID_C()
  END IF
END IF
IF(.NOT. uRandomExists) THEN
  Clock = IEOR(Clock, INT(ProcessID, KIND(Clock)))
  AuxilaryClock=Clock
  DO iSeed = 1, SeedSize
#if USE_MPI
    IF (nRandomSeeds.EQ.0) THEN
      AuxilaryClock=AuxilaryClock+PartMPI%MyRank
    ELSE IF(nRandomSeeds.GT.0) THEN
      AuxilaryClock=AuxilaryClock+(PartMPI%MyRank+1)*Seeds(iSeed)*37
    END IF
#else
    IF (nRandomSeeds.GT.0) THEN
      AuxilaryClock=AuxilaryClock+Seeds(iSeed)*37
    END IF
#endif
    IF (AuxilaryClock .EQ. 0) THEN
      AuxilaryClock = 104729
    ELSE
      AuxilaryClock = MOD(AuxilaryClock, 4294967296_8)
    END IF
    AuxilaryClock = MOD(AuxilaryClock * 279470273_8, 4294967291_8)
    GoodSeeds = INT(MOD(AuxilaryClock, INT(HUGE(0),KIND=8)), KIND(0))
    Seeds(iSeed) = GoodSeeds
  END DO
END IF
CALL RANDOM_SEED(PUT=Seeds)

END SUBROUTINE InitRandomSeed


END MODULE MOD_ParticleInit
