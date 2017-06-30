#include "boltzplatz.h"

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

INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE FinalizeParticles
  MODULE PROCEDURE FinalizeParticles
END INTERFACE

PUBLIC::InitParticles,FinalizeParticles
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticles()
!===================================================================================================================================
! Glue Subroutine for particle initialization 
!===================================================================================================================================
! MODULES
USE MOD_Globals!,       ONLY: MPIRoot,UNIT_STDOUT
USE MOD_ReadInTools
USE MOD_Particle_Vars,              ONLY: ParticlesInitIsDone, WriteMacroValues, nSpecies
USE MOD_part_emission,              ONLY: InitializeParticleEmission, InitializeParticleSurfaceflux
USE MOD_DSMC_Analyze,               ONLY: InitHODSMC
USE MOD_DSMC_Init,                  ONLY: InitDSMC
USE MOD_LD_Init,                    ONLY: InitLD
USE MOD_LD_Vars,                    ONLY: useLD
USE MOD_DSMC_Vars,                  ONLY: useDSMC, DSMC, DSMC_HOSolution,HODSMC
USE MOD_Mesh_Vars,                  ONLY : nElems
USE MOD_InitializeBackgroundField,  ONLY:InitializeBackgroundField
USE MOD_PICInterpolation_Vars,      ONLY: useBGField
#ifdef MPI
USE MOD_Particle_MPI,               ONLY:InitParticleCommSize
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef MPI
#endif
!===================================================================================================================================
IF(ParticlesInitIsDone)THEN
   SWRITE(*,*) "InitParticles already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES ...'

CALL InitializeVariables()
IF(useBGField) CALL InitializeBackgroundField()

CALL InitializeParticleEmission()
CALL InitializeParticleSurfaceflux()

IF (useDSMC) THEN
  CALL  InitDSMC()
  IF (useLD) CALL InitLD
ELSE IF (WriteMacroValues) THEN
  DSMC%CalcSurfaceVal  = .FALSE.
  DSMC%ElectronicModel = .FALSE.
  DSMC%OutputMeshInit  = .FALSE.
  DSMC%OutputMeshSamp  = .FALSE.
END IF

#ifdef MPI
! has to be called AFTER InitializeVariables and InitDSMC 
CALL InitParticleCommSize()
#endif

IF(useDSMC .OR. WriteMacroValues) THEN
! definition of DSMC sampling values
  DSMC%SampNum = 0
  HODSMC%SampleType = TRIM(GETSTR('DSMC-HOSampling-Type','cell_mean'))
  IF (TRIM(HODSMC%SampleType).EQ.'cell_mean') THEN
    HODSMC%nOutputDSMC = 1
    SWRITE(*,*) 'DSMCHO output order is set to 1 for sampling type cell_mean!'
  ELSE
    HODSMC%nOutputDSMC = GETINT('Particles-DSMC-OutputOrder','1')
  END IF     
  ALLOCATE(DSMC_HOSolution(1:11,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,1:nElems,1:nSpecies))         
  DSMC_HOSolution = 0.0
  CALL InitHODSMC()
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
USE MOD_Globals!, ONLY:MPIRoot,UNIT_STDOUT,myRank,nProcessors
USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_Particle_Vars!, ONLY: 
USE MOD_Particle_Boundary_Vars,ONLY:PartBound,nPartBound,nAdaptiveBC
USE MOD_Particle_Mesh_Vars    ,ONLY:NbrOfRegions,RegionBounds
USE MOD_Mesh_Vars,             ONLY:nElems, BoundaryName,BoundaryType, nBCs
USE MOD_Particle_Surfaces_Vars,ONLY:BCdata_auxSF
USE MOD_DSMC_Vars,             ONLY:useDSMC
USE MOD_Particle_Output_Vars,  ONLY:WriteFieldsToVTK, OutputMesh
USE MOD_part_MPFtools,         ONLY:DefinePolyVec, DefineSplitVec
USE MOD_PICInterpolation,      ONLY:InitializeInterpolation
USE MOD_PICInit,               ONLY:InitPIC
USE MOD_Particle_Mesh,         ONLY:InitFIBGM,MapRegionToElem
USE MOD_Particle_Tracking_Vars,ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,     ONLY:SafetyFactor,halo_eps_velo,PartMPI
USE MOD_part_pressure,         ONLY:ParticlePressureIni,ParticlePressureCellIni
#if defined(IMEX) || defined (IMPA)
USE MOD_TimeDisc_Vars,         ONLY: nRKStages
#endif /*IMEX*/
#ifdef MPI
USE MOD_Particle_MPI,          ONLY: InitEmissionComm
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iInit, iPartBound, iSeed, iCC
INTEGER               :: SeedSize, iPBC, iBC, iSwaps, iRegions, iExclude
INTEGER               :: ALLOCSTAT
CHARACTER(32)         :: hilf , hilf2, hilf3
CHARACTER(200)        :: tmpString
LOGICAL               :: TrueRandom, PartDens_OnlyInit                                                  !
INTEGER,ALLOCATABLE   :: iSeeds(:)
REAL                  :: iRan, aVec, bVec   ! random numbers for random vectors
REAL                  :: lineVector(3), v_drift_line, A_ins
INTEGER               :: iVec, MaxNbrOfSpeciesSwaps
LOGICAL                       :: exitTrue
#ifdef MPI
#endif
!===================================================================================================================================
! Read print flags
printRandomSeeds = GETLOGICAL('printRandomSeeds','.TRUE.')
! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
#if defined(LSERK)
!print*, "SFSDRWE#"
ALLOCATE(Pt_temp(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT)  
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
Pt_temp=0.
#endif 

#if defined(IMEX) 
ALLOCATE(PartStage(1:PDM%maxParticleNumber,1:6,1:nRKStages-1), STAT=ALLOCSTAT)  ! save memory
!ALLOCATE(PartStage(1:PDM%maxParticleNumber,1:6,1:nRKStages), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate ParStage arrays!')
END IF
ALLOCATE(PartStateN(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT)  
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate PartStateN arrays!')
END IF
#endif /* IMEX */

#ifdef IMPA
#if (PP_TimeDiscMethod!=110)
ALLOCATE(PartStage(1:PDM%maxParticleNumber,1:6,1:nRKStages-1), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate ParStage arrays!')
END IF
ALLOCATE(PartStateN(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT)  
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartStateN arrays!')
END IF
#endif
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
ALLOCATE(Norm2_F_PartX0(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate Norm2_F_PartX0 arrays!')
END IF
ALLOCATE(Norm2_F_PartXk(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate Norm2_F_PartXk arrays!')
END IF
ALLOCATE(Norm2_F_PartXk_old(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'Cannot allocate Norm2_F_PartXk_old arrays!')
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

#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
ALLOCATE(PartIsImplicit(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)  ! save memory
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate PartIsImplicit arrays!')
END IF
PartIsImplicit=.FALSE.
! ALLOCATE(StagePartPos(1:PDM%maxParticleNumber,1:3) &
!         ,PEM%StageElement(1:PDM%maxParticleNumber) ,  STAT=ALLOCSTAT) 
! IF (ALLOCSTAT.NE.0) THEN
!   CALL abort(&
! __STAMP__&
!   ,' Cannot allocate the stage position and element arrays!')
! END IF
! StagePartPos=0
! PEM%StageElement=0
#endif

IF(DoRefMapping)THEN
  ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,' Cannot allocate partposref!')
  PartPosRef=-888.
END IF

! predefine random vectors
NumRanVec = GETINT('Particles-NumberOfRandomVectors','100000')
IF ((usevMPF).OR.(useDSMC)) THEN
  ALLOCATE(RandomVec(NumRanVec, 3))
  RandomVec = 0
  DO iVec = 1, NumRanVec  ! calculation of NumRanVec different Vectors
    CALL RANDOM_NUMBER(iRan)
    bVec              = 1 - 2*iRan
    aVec              = SQRT(1 - bVec**2)
    RandomVec(iVec,1) = bVec
    CALL RANDOM_NUMBER(iRan)
    bVec              = Pi *2 * iRan
    RandomVec(iVec,2) = aVec * COS(bVec)
    RandomVec(iVec,3) = aVec * SIN(bVec)
  END DO
END IF

ALLOCATE(PartState(1:PDM%maxParticleNumber,1:6)       , &
         LastPartPos(1:PDM%maxParticleNumber,1:3)     , &
         Pt(1:PDM%maxParticleNumber,1:3)              , &
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
PDM%dtFracPush(1:PDM%maxParticleNumber) = .FALSE.
PDM%IsNewPart(1:PDM%maxParticleNumber) = .FALSE.
LastPartPos(1:PDM%maxParticleNumber,1:3)    = 0.
PartState=0.
Pt=0.
PartSpecies        = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)=0

! initialize of adsorbed particle tracking 
KeepWallParticles = .FALSE.

nSpecies = GETINT('Part-nSpecies','1')


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
  END IF
  ALLOCATE(PartMPF(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
__STAMP__&
    ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
  END IF
END IF

!WriteOutputMesh (=vtk mesh at start, seperate mesh for each proc
OutputMesh = GETLOGICAL('Part-WriteOutputMesh','.FALSE.')
           
! output of macroscopic values
WriteMacroValues = GETLOGICAL('Part-WriteMacroValues','.FALSE.')
MacroValSamplIterNum = GETINT('Part-IterationForMacroVal','1')
!ParticlePushMethod = TRIM(GETSTR('Part-ParticlePushMethod','boris_leap_frog_scheme')
WriteFieldsToVTK = GETLOGICAL('Part-WriteFieldsToVTK','.FALSE.')

!!!! Logicals for Constant Pressure in Cells
! are particles to be ADDED to cells in order to reach constant pressure? Default YES
PartPressAddParts = GETLOGICAL('Part-ConstPressAddParts','.TRUE.')
! are particles to be REMOVED from cells in order to reach constant pressure? Default NO
PartPressRemParts = GETLOGICAL('Part-ConstPressRemParts','.FALSE.')

! Read particle species data
!nSpecies = CNTSTR('Part-Species-SpaceIC')

IF (nSpecies.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF
PartPressureCell = .FALSE.
ALLOCATE(Species(1:nSpecies))

DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(hilf)//'-nInits','0')
  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits)) 
  DO iInit = 0, Species(iSpec)%NumberOfInits
    ! set help characters
    IF(iInit.EQ.0)THEN
      hilf2=TRIM(hilf)
    ELSE ! iInit >0
      WRITE(UNIT=hilf2,FMT='(I2)') iInit
      hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
    END IF ! iInit
    ! get species values // only once
    IF(iInit.EQ.0)THEN
      !General Species Values
      Species(iSpec)%ChargeIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-ChargeIC','0.')
      Species(iSpec)%MassIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-MassIC','0.')
      Species(iSpec)%MacroParticleFactor   = GETREAL('Part-Species'//TRIM(hilf2)//'-MacroParticleFactor','1.')
#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
      Species(iSpec)%IsImplicit            = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-IsImplicit','.FALSE.')
#endif
    END IF ! iInit
    ! get emission and init data
    Species(iSpec)%Init(iInit)%UseForInit           = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForInit','.TRUE.')
    Species(iSpec)%Init(iInit)%UseForEmission       = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForEmission','.TRUE.')
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-SpaceIC','cuboid'))
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution','constant'))
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'tangential_constant')THEN
      Species(iSpec)%Init(iInit)%Rotation        = GETINT('Part-Species'//TRIM(hilf2)//'-rotation','1')
      Species(iSpec)%Init(iInit)%VelocitySpread  = GETREAL('Part-Species'//TRIM(hilf2)//'-velocityspread','0.')
      IF(Species(iSpec)%Init(iInit)%VelocitySpread.LT.0. .OR. Species(iSpec)%Init(iInit)%VelocitySpread.GT.1.) CALL abort(&
__STAMP__&
          ,' Wrong input parameter for VelocitySpread in [0;1].')
      Species(iSpec)%Init(iInit)%VelocitySpreadMethod  = GETINT('Part-Species'//TRIM(hilf2)//'-velocityspreadmethod','0')
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
    Species(iSpec)%Init(iInit)%ConstantPressure      = GETREAL('Part-Species'//TRIM(hilf2)//'-ConstantPressure','0.')
    Species(iSpec)%Init(iInit)%ConstPressureRelaxFac = GETREAL('Part-Species'//TRIM(hilf2)//'-ConstPressureRelaxFac','1.')
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity','0.')
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(hilf2)//'-ParticleEmissionType','2')
    Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL('Part-Species'//TRIM(hilf2)//'-ParticleEmission','0.')
    Species(iSpec)%Init(iInit)%NSigma                = GETREAL('Part-Species'//TRIM(hilf2)//'-NSigma','10.')
    Species(iSpec)%Init(iInit)%NumberOfExcludeRegions= GETINT('Part-Species'//TRIM(hilf2)//'-NumberOfExcludeRegions','0')
    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus = 0
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell-juettner') THEN
      Species(iSpec)%Init(iInit)%MJxRatio       = GETREAL('Part-Species'//TRIM(hilf2)//'-MJxRatio','0')
      Species(iSpec)%Init(iInit)%MJyRatio       = GETREAL('Part-Species'//TRIM(hilf2)//'-MJyRatio','0')
      Species(iSpec)%Init(iInit)%MJzRatio       = GETREAL('Part-Species'//TRIM(hilf2)//'-MJzRatio','0')
    END IF
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'weibel') THEN
      Species(iSpec)%Init(iInit)%WeibelVeloPar       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPar','0')
      Species(iSpec)%Init(iInit)%WeibelVeloPer       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPer','0')
    END IF
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ. 'OneD-twostreaminstabilty') THEN
      Species(iSpec)%Init(iInit)%OneDTwoStreamVelo   = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamVelo','0')
      Species(iSpec)%Init(iInit)%OneDTwoStreamTransRatio = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamTransRatio','0')
    END IF

    !----------- various checks/calculations after read-in of Species(i)%Init(iInit)%-data ----------------------------------!
    !--- Check if Initial ParticleInserting is really used
    IF ( ((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.1).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2)) &
      .AND.(Species(iSpec)%Init(iInit)%UseForInit) ) THEN
      IF ( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0) &
      .AND. (ABS(Species(iSpec)%Init(iInit)%PartDensity).LE.0.) ) THEN
        Species(iSpec)%Init(iInit)%UseForInit=.FALSE.
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled as neither ParticleNumber"
        SWRITE(*,*) "nor PartDensity detected for Species, Init ", iSpec, iInit
      END IF
    END IF
    !--- cuboid-/cylinder-height calculation from v and dt
    IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
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
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cuboid') &
          .AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cylinder')) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is only supported for cuboid or cylinder!')
      IF (Species(iSpec)%Init(iInit)%UseForInit) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is not supported for initial ParticleInserting!')
    END IF
    !--- virtual pre-insertion (vpi) checks and calculations
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
      IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) .AND. (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2) ) &
        CALL abort(&
__STAMP__&
        ,' Wrong emission-type for virtual Pre-Inserting region!')
      IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell_lpn') &
        CALL abort(&
__STAMP__&
        ,' Only maxwell_lpn is implemened as velocity-distribution for virtual Pre-Inserting region!')
      IF (Species(iSpec)%Init(iInit)%UseForInit) &
        CALL abort(&
__STAMP__&
          ,' virtual Pre-Inserting is not supported for initial ParticleInserting. Use additional Init!')
      !-- Virtual Pre-Inserting is used correctly !
      Species(iSpec)%Init(iInit)%VirtPreInsert = .TRUE.
      SWRITE(*,*) "Virtual Pre-Inserting is used for Species, Init ", iSpec, iInit
      IF (Species(iSpec)%Init(iInit)%PartDensity .EQ. 0.) THEN
        SWRITE(*,*) "WARNING: If VPI-BC is open, a backflow might not be compensated"
        SWRITE(*,*) "         (use PartDensity instead of ParticleEmission)!"
      END IF
      Species(iSpec)%Init(iInit)%vpiDomainType = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-vpiDomainType','perpendicular_extrusion'))
      SELECT CASE ( TRIM(Species(iSpec)%Init(iInit)%vpiDomainType) )
      CASE ( 'freestream' )
        IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .NE. 'cuboid_vpi' ) THEN
          CALL abort(&
__STAMP__&
            ,' Only cuboid_vpi is supported for a freestream vpiDomainType! (Use default vpiDomainType for cylinder.)')
        ELSE
          Species(iSpec)%Init(iInit)%vpiBVBuffer(1) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV1BufferNeg','.TRUE.')
          Species(iSpec)%Init(iInit)%vpiBVBuffer(2) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV1BufferPos','.TRUE.')
          Species(iSpec)%Init(iInit)%vpiBVBuffer(3) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV2BufferNeg','.TRUE.')
          Species(iSpec)%Init(iInit)%vpiBVBuffer(4) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV2BufferPos','.TRUE.')
        END IF
      CASE ( 'orifice' )
        Species(iSpec)%Init(iInit)%vpiBVBuffer = .TRUE.
        IF ( ABS(Species(iSpec)%Init(iInit)%Radius2IC) .GT. 0. ) THEN
          CALL abort(&
__STAMP__&
            ,' Annular orifice is not implemented yet!')
        END IF
      CASE ( 'perpendicular_extrusion' )
        Species(iSpec)%Init(iInit)%vpiBVBuffer = .TRUE. !dummy
      CASE DEFAULT
        CALL abort(&
__STAMP__&
,'vpiDomainType is not implemented!')
      END SELECT
      !--
    ELSE
      Species(iSpec)%Init(iInit)%VirtPreInsert = .FALSE.
    END IF
    !--- integer check for ParticleEmissionType 2
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         ((Species(iSpec)%Init(iInit)%ParticleEmission-INT(Species(iSpec)%Init(iInit)%ParticleEmission)).NE.0)) THEN
       CALL abort(&
__STAMP__&
       ,' If ParticleEmissionType = 2 (parts per iteration), ParticleEmission has to be an integer number')
    END IF
    !--- flag for cell-based constant-pressure-EmiTypes
    IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 4).OR. &
        (Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6)) PartPressureCell = .TRUE.
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
        .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
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
    !--- read stuff for ExcludeRegions and normalize/calculate corresponding vectors
    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      ALLOCATE(Species(iSpec)%Init(iInit)%ExcludeRegion(1:Species(iSpec)%Init(iInit)%NumberOfExcludeRegions)) 
      IF (((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') &
       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) &
      .OR.((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi'))) THEN
        DO iExclude=1,Species(iSpec)%Init(iInit)%NumberOfExcludeRegions
          WRITE(UNIT=hilf3,FMT='(I2)') iExclude
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
          IF ((TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).EQ.'cuboid') .OR. &
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
          ELSE IF ( (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cuboid') .AND. &
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
,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
      END IF
    END IF
    !--- stuff for calculating ParticleEmission/InitialParticleNumber from PartDensity when this value is not used for LD-stuff
    !                                                                                  (additional not-LD checks might be necassary)
    PartDens_OnlyInit=.FALSE.
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.).AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'LD_insert')) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) THEN
        IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2) .AND. (Species(iSpec)%Init(iInit)%UseForInit) ) THEN
          PartDens_OnlyInit=.TRUE.
        ELSE
          CALL abort(&
__STAMP__&
            , 'PartDensity without LD is only supported for EmiType1 or initial ParticleInserting with EmiType1/2!')
        END IF
      END IF
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
        IF  ((((TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'emmert') ) THEN
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
                  ,'PartDensity without LD is only supported for CalcHeightFromDt, vpi, or initial ParticleInserting!')
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
              IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * Species(iSpec)%Init(iInit)%CuboidHeightIC * A_ins)
              ELSE !cylinder
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * Species(iSpec)%Init(iInit)%CylinderHeightIC * A_ins)
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
          ,'Only const. or maxwell(_lpn) is supported as velocityDistr. for PartDensity without LD!')
        END IF
      ELSE IF (Species(iSpec)%Init(iInit)%VirtPreInsert) THEN
        IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
               CALL abort(&
__STAMP__&
          ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
        ELSE
          SWRITE(*,*) "PartDensity is used for VPI of Species, Init ", iSpec, iInit !Value is calculated inside SetParticlePostion!
        END IF
      ELSE
        CALL abort(&
__STAMP__&
        ,'PartDensity without LD is only supported for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
      END IF
    END IF
    !--- determine StartnumberOfInits (start loop index, e.g., for emission loops)
    IF(iInit.EQ.0)THEN
      !!!for new case: check if to be included!!!
      IF((( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0)&
        .AND.(Species(iSpec)%Init(iInit)%ParticleEmission.EQ.0.) )  &
        .AND.(Species(iSpec)%Init(iInit)%PartDensity.EQ.0.) )       &
        .AND.(Species(iSpec)%Init(iInit)%ConstantPressure.EQ.0.)    &
        .AND.(Species(iSpec)%NumberOfInits.GT.0))       THEN 
        Species(iSpec)%StartnumberOfInits = 1 ! only new style paramaters defined (Part-Species(i)-Init(iInit)-***)
      ELSE
        Species(iSpec)%StartnumberOfInits = 0 ! old style parameters has been defined for inits/emissions (Part-Species(i)-***)
      END IF
      SWRITE(*,*) "StartnumberOfInits of Species ", iSpec, " = ", Species(iSpec)%StartnumberOfInits
    END IF ! iInit .EQ.0

  END DO ! iInit
END DO ! iSpec 

! Which Lorentz boost method should be used?
PartLorentzType = GETINT('Part-LorentzType','3')

! Read in boundary parameters
nPartBound = GETINT('Part-nBounds','1.')
IF (nPartBound.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nPartBound .LE. 0:', nPartBound)
END IF
ALLOCATE(PartBound%SourceBoundName(1:nPartBound))
ALLOCATE(PartBound%TargetBoundCond(1:nPartBound))
ALLOCATE(PartBound%MomentumACC(1:nPartBound))
ALLOCATE(PartBound%WallTemp(1:nPartBound))
ALLOCATE(PartBound%TransACC(1:nPartBound))
ALLOCATE(PartBound%VibACC(1:nPartBound))
ALLOCATE(PartBound%RotACC(1:nPartBound))
ALLOCATE(PartBound%Resample(1:nPartBound))
ALLOCATE(PartBound%WallVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientCondition(1:nPartBound))
ALLOCATE(PartBound%AmbientTemp(1:nPartBound))
ALLOCATE(PartBound%AmbientMeanPartMass(1:nPartBound))
ALLOCATE(PartBound%AmbientBeta(1:nPartBound))
ALLOCATE(PartBound%AmbientVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientDens(1:nPartBound))
ALLOCATE(PartBound%AmbientDynamicVisc(1:nPartBound))
ALLOCATE(PartBound%AmbientThermalCond(1:nPartBound))

ALLOCATE(PartBound%Adaptive(1:nPartBound))
ALLOCATE(PartBound%AdaptiveType(1:nPartBound))
ALLOCATE(PartBound%AdaptiveTemp(1:nPartBound))
ALLOCATE(PartBound%AdaptivePressure(1:nPartBound))
nAdaptiveBC = 0
PartBound%Adaptive(:) = .FALSE.
PartBound%AdaptiveType(:) = -1
PartBound%AdaptiveTemp(:) = -1.
PartBound%AdaptivePressure(:) = -1.

ALLOCATE(PartBound%Voltage(1:nPartBound))
ALLOCATE(PartBound%Voltage_CollectCharges(1:nPartBound))
PartBound%Voltage_CollectCharges(:)=0.
ALLOCATE(PartBound%NbrOfSpeciesSwaps(1:nPartBound))
!--determine MaxNbrOfSpeciesSwaps for correct allocation
MaxNbrOfSpeciesSwaps=0
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I2)') iPartBound
  PartBound%NbrOfSpeciesSwaps(iPartBound)= GETINT('Part-Boundary'//TRIM(hilf)//'-NbrOfSpeciesSwaps','0')
  MaxNbrOfSpeciesSwaps=max(PartBound%NbrOfSpeciesSwaps(iPartBound),MaxNbrOfSpeciesSwaps)
END DO
IF (MaxNbrOfSpeciesSwaps.gt.0) THEN
  ALLOCATE(PartBound%ProbOfSpeciesSwaps(1:nPartBound))
  ALLOCATE(PartBound%SpeciesSwaps(1:2,1:MaxNbrOfSpeciesSwaps,1:nPartBound))
END IF
!--
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I2)') iPartBound
  tmpString = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-Condition','open'))
  SELECT CASE (TRIM(tmpString))
  CASE('open')
     PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
     PartBound%AmbientCondition(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-AmbientCondition','.FALSE.')
     IF(PartBound%AmbientCondition(iPartBound)) THEN
       PartBound%AmbientTemp(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientTemp','0')
       PartBound%AmbientMeanPartMass(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientMeanPartMass','0')
       PartBound%AmbientBeta(iPartBound) = &
       SQRT(PartBound%AmbientMeanPartMass(iPartBound)/(2*BoltzmannConst*PartBound%AmbientTemp(iPartBound)))
       PartBound%AmbientVelo(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-AmbientVelo',3,'0. , 0. , 0.')
       PartBound%AmbientDens(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientDens','0')
       PartBound%AmbientDynamicVisc(iPartBound)=&
           GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientDynamicVisc','1.72326582572253E-5') ! N2:T=288K
       PartBound%AmbientThermalCond(iPartBound)=&
           GETREAL('Part-Boundary'//TRIM(hilf)//'-AmbientThermalCond','2.42948500556027E-2') ! N2:T=288K
     END IF
     PartBound%Adaptive(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Adaptive','.FALSE.')
     IF(PartBound%Adaptive(iPartBound)) THEN
       nAdaptiveBC = nAdaptiveBC + 1
       PartBound%AdaptiveType(iPartBound) = GETINT('Part-Boundary'//TRIM(hilf)//'-AdaptiveType','2')
       PartBound%AdaptiveTemp(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AdaptiveTemp','0.')
       PartBound%AdaptivePressure(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilf)//'-AdaptivePressure','0.')
       IF (PartBound%AdaptiveTemp(iPartBound)*PartBound%AdaptivePressure(iPartBound).EQ.0.) THEN
         CALL abort(&
__STAMP__&
,'Error during ParticleBoundary init: Part-Boundary'//TRIM(hilf)//'-AdaptiveTemp or -AdaptivePressure not defined')
       END IF
     END IF
     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
  CASE('reflective')
     PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
     PartBound%MomentumACC(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-MomentumACC','0')
     PartBound%WallTemp(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp','0')
     PartBound%TransACC(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-TransACC','0')
     PartBound%VibACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-VibACC','0')
     PartBound%RotACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotACC','0')
     PartBound%Resample(iPartBound)        = GETLOGICAL('Part-Boundary'//TRIM(hilf)//'-Resample','.FALSE.')
     PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(hilf)//'-WallVelo',3,'0. , 0. , 0.')
     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
     IF (PartBound%NbrOfSpeciesSwaps(iPartBound).gt.0) THEN  
       !read Species to be changed at wall (in, out), out=0: delete
       PartBound%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-Boundary'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
       DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(iPartBound)
         WRITE(UNIT=hilf2,FMT='(I2)') iSwaps
         PartBound%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
             GETINTARRAY('Part-Boundary'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
       END DO
     END IF
  CASE('periodic')
     PartBound%TargetBoundCond(iPartBound) = PartBound%PeriodicBC
  CASE('simple_anode')
     PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleAnodeBC
  CASE('simple_cathode')
     PartBound%TargetBoundCond(iPartBound) = PartBound%SimpleCathodeBC
  CASE('symmetric')
     PartBound%TargetBoundCond(iPartBound) = PartBound%SymmetryBC
     PartBound%WallVelo(1:3,iPartBound)    = (/0.,0.,0./)
  CASE DEFAULT
     SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
     CALL abort(&
__STAMP__&
         ,'Particle Boundary Condition does not exist')
  END SELECT
  PartBound%SourceBoundName(iPartBound) = TRIM(GETSTR('Part-Boundary'//TRIM(hilf)//'-SourceName'))
END DO

DEALLOCATE(PartBound%AmbientMeanPartMass)
DEALLOCATE(PartBound%AmbientTemp)
! Set mapping from field boundary to particle boundary index
ALLOCATE(PartBound%MapToPartBC(1:nBCs))
PartBound%MapToPartBC(:)=-10
DO iPBC=1,nPartBound
  DO iBC = 1, nBCs
    IF (BoundaryType(iBC,1).EQ.0) THEN
      PartBound%MapToPartBC(iBC) = -1 !there are no internal BCs in the mesh, they are just in the name list!
      SWRITE(*,*)"... PartBound",iPBC,"is internal bound, no mapping needed"
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

ALLOCATE(PEM%Element(1:PDM%maxParticleNumber), PEM%lastElement(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(&
__STAMP__&
  ,' Cannot allocate PEM arrays!')
END IF
IF (useDSMC.OR.PartPressureCell) THEN
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

!--- Read Manual Time Step
useManualTimeStep = .FALSE.
ManualTimeStep = GETREAL('Particles-ManualTimeStep', '0.0')
IF (ManualTimeStep.GT.0.0) THEN
  useManualTimeStep=.True.
END IF
#if (PP_TimeDiscMethod==201||PP_TimeDiscMethod==200)
  dt_part_ratio = GETREAL('Particles-dt_part_ratio', '3.8')
  overrelax_factor = GETREAL('Particles-overrelax_factor', '1.0')
#if (PP_TimeDiscMethod==200)
IF ( ALMOSTEQUAL(overrelax_factor,1.0) .AND. .NOT.ALMOSTEQUAL(dt_part_ratio,3.8) ) THEN
  overrelax_factor = dt_part_ratio !compatibility
END IF
#endif
#endif

!--- initialize randomization (= random if one or more seeds are 0 or random is wanted)
nrSeeds = GETINT('Part-NumberOfRandomSeeds','0')
IF (nrSeeds.GT.0) THEN
   ALLOCATE(seeds(1:nrSeeds))
   DO iSeed = 1, nrSeeds
      WRITE(UNIT=hilf,FMT='(I2)') iSeed
      seeds(iSeed) = GETINT('Particles-RandomSeed'//TRIM(hilf),'0')
   END DO
END IF

CALL RANDOM_SEED(Size = SeedSize)                       ! Check for number of needed Seeds
TrueRandom = .FALSE.                             ! FALSE for defined random seed

! to be stored in HDF5-state file!!!
IF (nrSeeds.GT.0) THEN
   IF (nrSeeds.NE.SeedSize) THEN
      IF(printRandomSeeds)THEN
        IPWRITE(UNIT_StdOut,*) 'Error: Number of seeds for RNG must be ',SeedSize
        IPWRITE(UNIT_StdOut,*) 'Random RNG seeds are used'
      END IF
      TrueRandom = .TRUE.
   END IF
   DO iSeed = 1, nrSeeds
      IF (Seeds(iSeed).EQ.0) THEN
         IF(printRandomSeeds)THEN
           IPWRITE(UNIT_StdOut,*) 'Error: ',SeedSize,' seeds for RNG must be defined'
           IPWRITE(UNIT_StdOut,*) 'Random RNG seeds are used'
         END IF
         TrueRandom = .TRUE.
      END IF
   END DO
ELSE
   TrueRandom = .TRUE.
END IF

IF (TrueRandom) THEN
   CALL RANDOM_SEED()
ELSE
#ifdef MPI
   Seeds(1:SeedSize) = Seeds(1:SeedSize)+PartMPI%MyRank
#endif
   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
END IF

ALLOCATE(iseeds(SeedSize))
iseeds(:)=0
CALL RANDOM_SEED(GET = iseeds(1:SeedSize))
! to be stored in HDF5-state file!!!
IF(printRandomSeeds)THEN
  IPWRITE(UNIT_stdOut,*) 'Random seeds in PIC_init:'
  DO iSeed = 1,SeedSize
     IPWRITE(UNIT_stdOut,*) iseeds(iSeed)
  END DO
END IF
DEALLOCATE(iseeds)
!DoZigguratSampling = GETLOGICAL('Particles-DoZigguratSampling','.FALSE.')
DoPoissonRounding = GETLOGICAL('Particles-DoPoissonRounding','.FALSE.')
DoTimeDepInflow   = GETLOGICAL('Particles-DoTimeDepInflow','.FALSE.')

DO iSpec = 1, nSpecies
  DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
    IF(Species(iSpec)%Init(iInit)%InflowRiseTime.GT.0.)THEN
      IF(.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow)  CALL CollectiveStop(&
__STAMP__, &
' Linearly ramping of inflow-number-of-particles is only possible with PoissonRounding or DoTimeDepInflow!')
    END IF      
  END DO ! iInit = 0, Species(iSpec)%NumberOfInits
END DO ! iSpec = 1, nSpecies


DelayTime = GETREAL('Part-DelayTime','0.')

!-- Read Flag for BGG via VTK-File (needed for particles in general as asked for in set velocity!)
useVTKFileBGG = GETLOGICAL('Particles-useVTKFileBGG','.FALSE.')
!-- Read Flag if warnings to be displayed for rejected velocities when virtual Pre-Inserting region (vpi) is used with PartDensity
OutputVpiWarnings = GETLOGICAL('Particles-OutputVpiWarnings','.FALSE.')


! init interpolation
CALL InitializeInterpolation() ! not any more required ! has to be called earliear
CALL InitPIC()
! always, because you have to construct a halo_eps region around each bc element

#ifdef MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT FIBGM...' 
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')
!-- Finalizing InitializeVariables
CALL InitFIBGM()
!CALL InitSFIBGM()
#ifdef MPI
CALL InitEmissionComm()
#endif /*MPI*/
#ifdef MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/

SWRITE(UNIT_StdOut,'(132("-"))')

!-- Read parameters for particle-data on region mapping

!-- Read parameters for region mapping
NbrOfRegions = GETINT('NbrOfRegions','0')
IF (NbrOfRegions .GT. 0) THEN
  ALLOCATE(RegionBounds(1:6,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I2)') iRegions
    RegionBounds(1:6,iRegions) = GETREALARRAY('RegionBounds'//TRIM(hilf2),6,'0. , 0. , 0. , 0. , 0. , 0.')
  END DO
END IF

IF (NbrOfRegions .GT. 0) THEN
  CALL MapRegionToElem()
  ALLOCATE(RegionElectronRef(1:3,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I2)') iRegions
    RegionElectronRef(1:3,iRegions) = GETREALARRAY('Part-RegionElectronRef'//TRIM(hilf2),3,'0. , 0. , 1.')
  END DO
END IF

exitTrue=.false.
DO iSpec = 1,nSpecies
  DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 3).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 5)) THEN
      CALL ParticlePressureIni()
      exitTrue=.true.
      EXIT
    END IF
  END DO
  IF (exitTrue) EXIT
END DO

exitTrue=.false.
DO iSpec = 1,nSpecies
  DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
    IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 4).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6)) THEN
      CALL ParticlePressureCellIni()
      exitTrue=.true.
      EXIT
    END IF
  END DO
  IF (exitTrue) EXIT
END DO




IF(enableParticleMerge) THEN
 CALL DefinePolyVec(vMPFMergePolyOrder) 
 CALL DefineSplitVec(vMPFMergeCellSplitOrder)
END IF

!-- Floating Potential
ALLOCATE(BCdata_auxSF(1:nPartBound))
DO iPartBound=1,nPartBound
  BCdata_auxSF(iPartBound)%SideNumber=-1 !init value when not used
  BCdata_auxSF(iPartBound)%GlobalArea=0.
  BCdata_auxSF(iPartBound)%LocalArea=0.
END DO
nDataBC_CollectCharges=0
nCollectChargesBCs = GETINT('PIC-nCollectChargesBCs','0')
IF (nCollectChargesBCs .GT. 0) THEN
#if !(defined (PP_HDG) && (PP_nVar==1))
  CALL abort(__STAMP__&
    , 'CollectCharges only implemented for electrostatic HDG!')
#endif
  ALLOCATE(CollectCharges(1:nCollectChargesBCs))
  DO iCC=1,nCollectChargesBCs
    WRITE(UNIT=hilf,FMT='(I2)') iCC
    CollectCharges(iCC)%BC = GETINT('PIC-CollectCharges'//TRIM(hilf)//'-BC','0')
    IF (CollectCharges(iCC)%BC.LT.1 .OR. CollectCharges(iCC)%BC.GT.nPartBound) THEN
      CALL abort(__STAMP__&
      , 'nCollectChargesBCs must be between 1 and nPartBound!')
    ELSE IF (BCdata_auxSF(CollectCharges(iCC)%BC)%SideNumber.EQ. -1) THEN !not set yet
      BCdata_auxSF(CollectCharges(iCC)%BC)%SideNumber=0
      nDataBC_CollectCharges=nDataBC_CollectCharges+1 !side-data will be set in InitializeParticleSurfaceflux!!!
    END IF
    CollectCharges(iCC)%NumOfRealCharges = GETREAL('PIC-CollectCharges'//TRIM(hilf)//'-NumOfRealCharges','0.')
    CollectCharges(iCC)%NumOfNewRealCharges = 0.
    CollectCharges(iCC)%ChargeDist = GETREAL('PIC-CollectCharges'//TRIM(hilf)//'-ChargeDist','0.')
  END DO !iCC
END IF !nCollectChargesBCs .GT. 0


END SUBROUTINE InitializeVariables


SUBROUTINE FinalizeParticles() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize particle variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars
!USE MOD_DSMC_Vars,                  ONLY: SampDSMC
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if defined(LSERK)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
SDEALLOCATE( Pt_temp)
#endif
!SDEALLOCATE(SampDSMC)
SDEALLOCATE(PartPosRef)
SDEALLOCATE(RandomVec)
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
!SDEALLOCATE(Species%Init)
SDEALLOCATE(Species)
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%MomentumACC)
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%TransACC)
SDEALLOCATE(PartBound%VibACC)
SDEALLOCATE(PartBound%RotACC)
SDEALLOCATE(PartBound%Resample)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%AmbientCondition)
SDEALLOCATE(PartBound%AmbientTemp)
SDEALLOCATE(PartBound%AmbientMeanPartMass)
SDEALLOCATE(PartBound%AmbientBeta)
SDEALLOCATE(PartBound%AmbientVelo)
SDEALLOCATE(PartBound%AmbientDens)
SDEALLOCATE(PartBound%AmbientDynamicVisc)
SDEALLOCATE(PartBound%AmbientThermalCond)
SDEALLOCATE(PartBound%Adaptive)
SDEALLOCATE(PartBound%AdaptiveType)
SDEALLOCATE(PartBound%AdaptiveTemp)
SDEALLOCATE(PartBound%AdaptivePressure)
SDEALLOCATE(PartBound%Voltage)
SDEALLOCATE(PartBound%Voltage_CollectCharges)
SDEALLOCATE(PartBound%NbrOfSpeciesSwaps)
SDEALLOCATE(PartBound%ProbOfSpeciesSwaps)
SDEALLOCATE(PartBound%SpeciesSwaps)
SDEALLOCATE(PartBound%MapToPartBC)
SDEALLOCATE(PEM%Element)
SDEALLOCATE(PEM%lastElement)
SDEALLOCATE(PEM%pStart)
SDEALLOCATE(PEM%pNumber)
SDEALLOCATE(PEM%pEnd)
SDEALLOCATE(PEM%pNext)
SDEALLOCATE(seeds)
SDEALLOCATE(RegionBounds)
SDEALLOCATE(RegionElectronRef)
END SUBROUTINE FinalizeParticles

END MODULE MOD_ParticleInit
