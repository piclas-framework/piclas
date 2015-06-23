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

PUBLIC::InitParticles
!===================================================================================================================================

CONTAINS


SUBROUTINE InitParticles()
!===================================================================================================================================
! Glue Subroutine for particle initialization 
!===================================================================================================================================
! MODULES
USE MOD_Globals!,       ONLY: MPIRoot,UNIT_STDOUT
USE MOD_ReadInTools
USE MOD_Particle_Vars, ONLY: ParticlesInitIsDone, WriteMacroValues, nSpecies
USE MOD_part_emission, ONLY: InitializeParticleEmission
USE MOD_DSMC_Init,     ONLY: InitDSMC
!USE MOD_LD_Init,       ONLY: InitLD
!USE MOD_LD_Vars,       ONLY: useLD
USE MOD_DSMC_Vars,     ONLY: useDSMC, DSMC, SampDSMC
USE MOD_Mesh_Vars,     ONLY : nElems
USE MOD_InitializeBackgroundField
USE MOD_PICInterpolation_Vars, ONLY: useBGField
#ifdef MPI
USE MOD_Particle_MPI,     ONLY:InitParticleCommSize
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
#ifdef MPI
CALL InitParticleCommSize()
#endif
IF(useBGField) CALL InitializeBackgroundField()
CALL InitializeParticleEmission()

IF (useDSMC) THEN
  CALL  InitDSMC()
  !IF (useLD) CALL InitLD
ELSE IF (WriteMacroValues) THEN
  DSMC%SampNum = 0
  ALLOCATE(SampDSMC(nElems,nSpecies))
  SampDSMC(1:nElems,1:nSpecies)%PartV(1)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV(2)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV(3)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(1) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(2) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(3) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartNum   = 0
  SampDSMC(1:nElems,1:nSpecies)%SimPartNum   = 0
  SampDSMC(1:nElems,1:nSpecies)%ERot      = 0
  SampDSMC(1:nElems,1:nSpecies)%EVib      = 0
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
USE MOD_Particle_Mesh_Vars,    ONLY: PartBound,nPartBound
USE MOD_Mesh_Vars,             ONLY: nElems, BoundaryName,BoundaryType, nBCs,nSides,BC
USE MOD_Restart_Vars,          ONLY: DoRestart
USE MOD_DSMC_Vars,             ONLY: useDSMC
USE MOD_Particle_Output_Vars,  ONLY: WriteFieldsToVTK, OutputMesh
USE MOD_part_MPFtools,         ONLY: DefinePolyVec, DefineSplitVec
USE MOD_PICInterpolation_Vars, ONLY: InterpolationType
USE MOD_PICInterpolation,      ONLY: InitializeInterpolation
USE MOD_PICInit,               ONLY: InitPIC
USE MOD_Particle_Mesh,         ONLY: InitFIBGM
USE MOD_Particle_Surfaces_Vars,ONLY: DoRefMapping
!USE MOD_Particle_Mesh_Vars,    ONLY:Geo
#ifdef MPI
!USE MOD_part_MPI_Vars,         ONLY: PMPIVAR
USE MOD_Particle_MPI_Vars,     ONLY: SafetyFactor,halo_eps_velo,PartMPI
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
  INTEGER               :: iSpec, iInit, iPartBound, iSeed,iSide
  INTEGER               :: SeedSize, iPBC, iBC, iSwaps
  INTEGER               :: ALLOCSTAT
  CHARACTER(32)         :: hilf , hilf2
  CHARACTER(200)        :: tmpString
  LOGICAL               :: TrueRandom                                                  !
  INTEGER,ALLOCATABLE   :: iSeeds(:)
  REAL                  :: iRan, aVec, bVec   ! random numbers for random vectors
  REAL                  :: lineVector(3), v_drift_line, A_ins
  INTEGER               :: iVec, MaxNbrOfSpeciesSwaps
#ifdef MPI
#endif
!===================================================================================================================================
!#ifdef MPI
!   PMPIVAR%COMM   = MPI_COMM_WORLD
!   PMPIVAR%iProc  = myRank
!   PMPIVAR%nProcs = nProcessors
!   CALL MPI_COMM_GROUP(PMPIVAR%COMM,PMPIVAR%GROUP,IERROR)
!   PMPIVAR%GROUPWORLD=PMPIVAR%GROUP
!!IPWRITE(UNIT_stdOut,*)'INIT: GROUPWORLD',PMPIVAR%GROUPWORLD
!#endif

! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* RK3 and RK4 only */
!print*, "SFSDRWE#"
ALLOCATE(Pt_temp(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT)  
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
Pt_temp=0.
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
         PDM%nextFreePosition(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
! always zero
PDM%ParticleInside(1:PDM%maxParticleNumber) = .FALSE.
LastPartPos(1:PDM%maxParticleNumber,1:3)    = 0.
PartState=0.
Pt=0.
PartSpecies        = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)=0

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
    vMPF_velocityDistribution = GETSTR('Part-vMPFvelocityDistribution','OVDR')
    ALLOCATE(vMPF_SpecNumElem(1:nElems,1:nSpecies))
  END IF
  ALLOCATE(PartMPF(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
  END IF
END IF

!WriteOutputMesh (=vtk mesh at start, seperate mesh for each proc
OutputMesh = GETLOGICAL('Part-WriteOutputMesh','.FALSE.')
           
! output of macroscopic values
WriteMacroValues = GETLOGICAL('Part-WriteMacroValues','.FALSE.')
MacroValSamplIterNum = GETINT('Part-IterationForMacroVal','1')
!ParticlePushMethod = GETSTR('Part-ParticlePushMethod','boris_leap_frog_scheme')
WriteFieldsToVTK = GETLOGICAL('Part-WriteFieldsToVTK','.FALSE.')

!!!! Logicals for Constant Pressure in Cells
! are particles to be ADDED to cells in order to reach constant pressure? Default YES
PartPressAddParts = GETLOGICAL('Part-ConstPressAddParts','.TRUE.')
! are particles to be REMOVED from cells in order to reach constant pressure? Default NO
PartPressRemParts = GETLOGICAL('Part-ConstPressRemParts','.FALSE.')

! Read particle species data
!nSpecies = CNTSTR('Part-Species-SpaceIC')

IF (nSpecies.LE.0) THEN
  CALL abort(__STAMP__&
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
    END IF ! iInit
    ! get emission and init data
    Species(iSpec)%Init(iInit)%UseForInit           = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForInit','.TRUE.')
    Species(iSpec)%Init(iInit)%UseForEmission       = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-UseForEmission','.TRUE.')
    Species(iSpec)%Init(iInit)%SpaceIC               = GETSTR('Part-Species'//TRIM(hilf2)//'-SpaceIC','cuboid')
    Species(iSpec)%Init(iInit)%velocityDistribution  = GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution','constant')
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
    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'weibel') THEN
      Species(iSpec)%Init(iInit)%WeibelVeloPar       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPar','0')
      Species(iSpec)%Init(iInit)%WeibelVeloPer       = GETREAL('Part-Species'//TRIM(hilf2)//'-WeibelVeloPer','0')
    END IF
    IF(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ. 'OneD-twostreaminstabilty') THEN
      Species(iSpec)%Init(iInit)%OneDTwoStreamVelo   = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamVelo','0')
      Species(iSpec)%Init(iInit)%OneDTwoStreamTransRatio = GETREAL('Part-Species'//TRIM(hilf2)//'-OneDTwoStreamTransRatio','0')
    END IF

    !-- various checks/calculations after read-in of Species(i)%Init(iInit)%-data
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
      IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1).AND.(Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2)) &
        CALL abort(__STAMP__&
        ,' Wrong emission-type for virtual Pre-Inserting region!')
      IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell_lpn') &
        CALL abort(__STAMP__&
        ,' Only maxwell_lpn is implemened as velocity-distribution for virtual Pre-Inserting region!')
      IF (Species(iSpec)%Init(iInit)%initialParticleNumber.GT.0) THEN
        CALL abort(__STAMP__&
          ,' initialParticleNumber does not work for virtual Pre-Inserting. Use additional Init!')
      ELSE
        Species(iSpec)%Init(iInit)%UseForInit = .FALSE. !Do not make particle_emission stuff during initialization
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled, as both VPI and"
        SWRITE(*,*) "initialParticleNumber=0 detected for Species, Init ", iSpec, iInit
      END IF
      Species(iSpec)%Init(iInit)%VirtPreInsert = .TRUE.
      SWRITE(*,*) "Virtual Pre-Inserting is used for Species, Init ", iSpec, iInit
      IF (Species(iSpec)%Init(iInit)%PartDensity .EQ. 0.) THEN
        SWRITE(*,*) "WARNING: If VPI-BC is open, a backflow might not be compensated"
        SWRITE(*,*) "         (use PartDensity instead of ParticleEmission)!"
      END IF
    ELSE
      Species(iSpec)%Init(iInit)%VirtPreInsert = .FALSE.
    END IF
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         ((Species(iSpec)%Init(iInit)%ParticleEmission-INT(Species(iSpec)%Init(iInit)%ParticleEmission)).NE.0)) THEN
       CALL abort(__STAMP__&
       ,' If ParticleEmissionType = 2 (parts per iteration), ParticleEmission has to be an integer number')
    END IF
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

    !--- stuff for calculating ParticleEmission from PartDensity
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.).AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'LD_insert')) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) &
        CALL abort(__STAMP__&
        , 'Only emission-type 1 is supported for PartDensity without LD!')
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
        IF  (((TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') ) THEN
          IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
            CALL abort(__STAMP__&
            ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
          ELSE
            !---calculation of Base-Area and corresponding component of VeloVecIC
            lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
              Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
            lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
              Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
            lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
              Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
            IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
               CALL abort(__STAMP__&
               ,'BaseVectors are parallel!')
            ELSE
              A_ins = SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + lineVector(3) * lineVector(3))
              lineVector = lineVector / A_ins
              v_drift_line = Species(iSpec)%Init(iInit)%VeloIC * &
                ( Species(iSpec)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(iSpec)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
                + Species(iSpec)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
              IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .EQ. 'cylinder' ) &
                A_ins = Pi * (Species(iSpec)%Init(iInit)%RadiusIC**2-Species(iSpec)%Init(iInit)%Radius2IC**2)
              !---calculation of particle flow (macroparticles/s) through boundary
              Species(iSpec)%Init(iInit)%ParticleEmission = Species(iSpec)%Init(iInit)%PartDensity &
                                                          / Species(iSpec)%MacroParticleFactor * v_drift_line * A_ins
              END IF
          END IF
        ELSE
          CALL abort(__STAMP__&
          ,'Only const. or maxwell(_lpn) is supported as velocityDistr. for PartDensity without LD!')
        END IF
      ELSE IF (Species(iSpec)%Init(iInit)%VirtPreInsert) THEN
        IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
               CALL abort(__STAMP__&
          ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
        ELSE
          SWRITE(*,*) "PartDensity is used for VPI of Species, Init ", iSpec, iInit
        END IF
      ELSE
        CALL abort(__STAMP__&
        ,'PartDensity without LD is only supported for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
      END IF
    END IF

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

!! Read parameter for FastInitBackgroundMesh (FIBGM)
!GEO%FIBGMdeltas(1:3)              = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
!GEO%FactorFIBGM(1:3)              = GETREALARRAy('Part-FactorFIBGM',3,'1. , 1. , 1.')
!GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

! Read in boundary parameters
nPartBound = GETINT('Part-nBounds','1.')
IF (nPartBound.LE.0) THEN
  CALL abort(__STAMP__&
  ,'ERROR: nPartBound .LE. 0:', nPartBound)
END IF
ALLOCATE(PartBound%SourceBoundName(1:nPartBound))
ALLOCATE(PartBound%TargetBoundCond(1:nPartBound))
ALLOCATE(PartBound%MomentumACC(1:nPartBound))
ALLOCATE(PartBound%WallTemp(1:nPartBound))
ALLOCATE(PartBound%TransACC(1:nPartBound))
ALLOCATE(PartBound%VibACC(1:nPartBound))
ALLOCATE(PartBound%RotACC(1:nPartBound))
ALLOCATE(PartBound%WallVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientCondition(1:nPartBound))
ALLOCATE(PartBound%AmbientTemp(1:nPartBound))
ALLOCATE(PartBound%AmbientMeanPartMass(1:nPartBound))
ALLOCATE(PartBound%AmbientBeta(1:nPartBound))
ALLOCATE(PartBound%AmbientVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientDens(1:nPartBound))
ALLOCATE(PartBound%AmbientDynamicVisc(1:nPartBound))
ALLOCATE(PartBound%AmbientThermalCond(1:nPartBound))

ALLOCATE(PartBound%Voltage(1:nPartBound))
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
  tmpString = GETSTR('Part-Boundary'//TRIM(hilf)//'-Condition','open')
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
     PartBound%Voltage(iPartBound)         = GETREAL('Part-Boundary'//TRIM(hilf)//'-Voltage','0')
  CASE('reflective')
     PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
     PartBound%MomentumACC(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf)//'-MomentumACC','0')
     PartBound%WallTemp(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-WallTemp','0')
     PartBound%TransACC(iPartBound)        = GETREAL('Part-Boundary'//TRIM(hilf)//'-TransACC','0')
     PartBound%VibACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-VibACC','0')
     PartBound%RotACC(iPartBound)          = GETREAL('Part-Boundary'//TRIM(hilf)//'-RotACC','0')
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
  CASE DEFAULT
     SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
     CALL abort(__STAMP__&
         ,'Particle Boundary Condition does not exist')
  END SELECT
  PartBound%SourceBoundName(iPartBound) = GETSTR('Part-Boundary'//TRIM(hilf)//'-SourceName')
END DO
DEALLOCATE(PartBound%AmbientMeanPartMass)
DEALLOCATE(PartBound%AmbientTemp)
! Set mapping from field boundary to particle boundary index
ALLOCATE(PartBound%Map(1:nBCs))
PartBound%Map(:)=-10
DO iPBC=1,nPartBound
  DO iBC = 1, nBCs
    IF (BoundaryType(iBC,1).EQ.0) THEN
      PartBound%Map(iBC) = 1
      SWRITE(*,*)"PartBound",iPBC,"is internal bound, no mapping needed"
    END IF
    IF (TRIM(BoundaryName(iBC)).EQ.TRIM(PartBound%SourceBoundName(iPBC))) THEN
      PartBound%Map(iBC) = PartBound%TargetBoundCond(iPBC) 
      SWRITE(*,*)"Mapped PartBound",iPBC,"on FieldBound",BoundaryType(iBC,1),",i.e.:",TRIM(BoundaryName(iBC))
    END IF
  END DO
END DO
! Errorhandler for PartBound-Types that could not be mapped to the 
! FieldBound-Types.
DO iBC = 1,nBCs
  IF (PartBound%Map(iBC).EQ.-10) THEN
    CALL abort(__STAMP__&
    ,' PartBound%Map for Boundary is not set. iBC: :',iBC)
  END IF
END DO

!ALLOCATE(PartBound%SideBCType(1:nSides) )
!DO iSide=1,nSides
!  IF(iSide.LE.nBCSides)THEN
!    SELECT CASE( PartBound%Map(BC(SideID)))
!    CASE(PartBound%OpenBC)
!      PartBound%SideBCType(iSide)=PartBound%OpenBC
!    CASE(PartBound%ReflectiveBC)
!      PartBound%SideBCType(iSide)=PartBound%ReflectiveBC
!    CASE(PartBound%SimpleCathodeBC)
!      PartBound%SideBCType(iSide)=PartBound%SimpleCathodeBC
!    CASE(PartBound%SimpleCathodeBC)
!      PartBound%SideBCType(iSide)=PartBound%
!  ELSE
!    IF(BC(SideID).EQ.1) PartBound%SideBCType(iSide)=PartBound%PeriodicBC
!  END IF
!END DO ! iSide=1,nSides

ALLOCATE(PEM%Element(1:PDM%maxParticleNumber), PEM%lastElement(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(__STAMP__&
  ,' Cannot allocate PEM arrays!')
END IF
IF (useDSMC.OR.PartPressureCell) THEN
  ALLOCATE(PEM%pStart(1:nElems)                         , &
           PEM%pNumber(1:nElems)                        , &
           PEM%pEnd(1:nElems)                           , &
           PEM%pNext(1:PDM%maxParticleNumber)           , STAT=ALLOCSTAT) 
           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    , ' Cannot allocate DSMC PEM arrays!')
  END IF
END IF
IF (useDSMC) THEN
  ALLOCATE(PDM%PartInit(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
           !PDM%nextUsedPosition(1:PDM%maxParticleNumber)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__&
    ,' Cannot allocate DSMC PEM arrays!')
  END IF
END IF

!--- Read Manual Time Step
useManualTimeStep = .FALSE.
ManualTimeStep = GETREAL('Particles-ManualTimeStep', '0.0')
IF (ManualTimeStep.GT.0.0) THEN
  useManualTimeStep=.True.
END IF

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

IF (nrSeeds.GT.0) THEN
   IF (nrSeeds.NE.SeedSize) THEN
      IPWRITE(UNIT_stdOut,*) 'Error: Number of seeds for RNG must be ',SeedSize
      IPWRITE(UNIT_stdOut,*) 'Random RNG seeds are used'
      TrueRandom = .TRUE.
   END IF
   DO iSeed = 1, nrSeeds
      IF (Seeds(iSeed).EQ.0) THEN
         IPWRITE(UNIT_stdOut,*) 'Error: ',SeedSize,' seeds for RNG must be defined'
         IPWRITE(UNIT_stdOut,*) 'Random RNG seeds are used'
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
IPWRITE(UNIT_stdOut,*) 'Random seeds in PIC_init:'
DO iSeed = 1,SeedSize
   IPWRITE(UNIT_stdOut,*) iseeds(iSeed)
END DO
DEALLOCATE(iseeds)

DelayTime = GETREAL('Part-DelayTime','0.')
! init interpolation
CALL InitializeInterpolation() ! not any more required ! has to be called earliear
CALL InitPIC()
#ifdef MPI
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')
#endif /*MPI*/
!-- Finalizing InitializeVariables
CALL InitFIBGM()
#ifdef MPI
CALL InitEmissionComm()
#endif /*MPI*/
!CALL DomainUpdate()
IF(enableParticleMerge) THEN
 !IF (TRIM(InterpolationType).NE.'particle_position') CALL DefineElemT_inv()
 CALL DefinePolyVec(vMPFMergePolyOrder) 
 CALL DefineSplitVec(vMPFMergeCellSplitOrder)
END IF

END SUBROUTINE InitializeVariables

END MODULE MOD_ParticleInit
