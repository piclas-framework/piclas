MODULE MOD_Particle_Vars
!===================================================================================================================================
! Contains the Particles' variables (general for all modules: PIC, DSMC, FP)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL, PARAMETER       :: BoltzmannConst=1.380648813E-23                      ! Boltzmann constant [J/K] SI-Unit!
REAL                  :: ManualTimeStep                                      ! Manual TimeStep
LOGICAL               :: useManualTimeStep                                   ! Logical Flag for manual timestep. For consistency
                                                                             ! with IAG programming style
REAL                  :: dt_max_particles                                    ! Maximum timestep for particles (for static fields!)
REAL                  :: dt_maxwell                                          ! timestep for field solver (for static fields only!)
REAL                  :: dt_adapt_maxwell                                    ! adapted timestep for field solver dependent  
                                                                             ! on particle velocity (for static fields only!)
INTEGER               :: NextTimeStepAdjustmentIter                          ! iteration of next timestep change
INTEGER               :: MaxwellIterNum                                      ! number of iterations for the maxwell solver
INTEGER               :: WeirdElems                                          ! Number of Weird Elements (=Elements which are folded
                                                                             ! into themselves)
REAL    , ALLOCATABLE :: PartState(:,:)                                      ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE :: PartPosRef(:,:)                                     ! (1:3,1:NParts) particles pos mapped to -1|1 space
INTEGER , ALLOCATABLE :: PartPosGauss(:,:)                                   ! (1:NParts,1:3) Gauss point localization of particles
REAL    , ALLOCATABLE :: Pt(:,:)                                             ! Derivative of PartState (vx,xy,vz) only
                                                                             ! since temporal derivative of position
                                                                             ! is the velocity. Thus we can take 
                                                                             ! PartState(:,4:6) as Pt(1:3)
                                                                             ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE :: Pt_temp(:,:)                                        ! LSERK4 additional derivative of PartState
                                                                             ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE :: LastPartPos(:,:)                                    ! (1:NParts,1:3) with 2nd index: x,y,z
INTEGER , ALLOCATABLE :: PartSpecies(:)                                      ! (1:NParts) 
REAL    , ALLOCATABLE :: PartMPF(:)                                          ! (1:NParts) MacroParticleFactor by variable MPF
INTEGER               :: PartLorentzType
CHARACTER(LEN=256)    :: ParticlePushMethod                                  ! Type of PP-Method
INTEGER               :: nrSeeds                                             ! Number of Seeds for Random Number Generator
INTEGER , POINTER     :: seeds(:)                                 =>NULL()   ! Seeds for Random Number Generator

TYPE tConstPressure
  INTEGER                                :: nElemTotalInside                  ! Number of elements totally in Emission Particle  
  INTEGER                                :: nElemPartlyInside                 ! Number of elements partly in Emission Particle 
  INTEGER, ALLOCATABLE                   :: ElemTotalInside(:)                ! List of elements totally in Emission Particle 
                                                                              ! ElemTotalInside(1:nElemTotalInside)
  INTEGER, ALLOCATABLE                   :: ElemPartlyInside(:)               ! List of elements partly in Emission Particle 
                                                                              ! ElemTotalInside(1:nElemPartlyInside)
  INTEGER(2), ALLOCATABLE                :: ElemStat(:)                       ! Status of Element to Emission Particle Space
                                                                              ! ElemStat(nElem) = 1  -->  Element is totally insid
                                                                              !                 = 2  -->  Element is partly  insid
                                                                              !                 = 3  -->  Element is totally outsi
  REAL                                   :: OrthoVector(3)                    ! Vector othogonal on BaseVector1IC and BaseVector2
  REAL                                   :: Determinant                       ! Determinant for solving a 3x3 system of equations
                                                                              ! to see whether a point is inside a cuboid
  REAL                                   :: EkinInside                        ! Kinetic Energy in Emission-Area
  REAL                                   :: InitialTemp                       ! Initial MWTemerature
  REAL, ALLOCATABLE                      :: ConstPressureSamp(:,:)            ! ElemTotalInside(1:nElemTotalInside,1 = v_x
                                                                              !                                    2 = v_y
                                                                              !                                    3 = v_z
                                                                              !                                    4 = dens. [1/m3]
                                                                              !                                    5 = pressure
                                                                              !                                    6 = v of sound**2
END TYPE

TYPE tExcludeRegion
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  REAL                                   :: RadiusIC                         ! Radius for IC circle
  REAL                                   :: Radius2IC                        ! Radius2 for IC cylinder (ring)
  REAL                                   :: NormalIC(3)                      ! Normal / Orientation of cylinder (altern. to BV1/2)
  REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                                   :: CuboidHeightIC                   ! third measure of cuboid
                                                                             ! (set 0 for flat rectangle),
                                                                             ! negative value = opposite direction
  REAL                                   :: CylinderHeightIC                 ! third measure of cylinder
                                                                             ! (set 0 for flat circle),
                                                                             ! negative value = opposite direction
  REAL                                   :: ExcludeBV_lenghts(2)                    ! lenghts of BV1/2 (to be calculated)
END TYPE

TYPE tInit                                                                   ! Particle Data for each init emission for each species
  !Specific Emission/Init values
  LOGICAL                                :: UseForInit                       ! Use Init/Emission for init.?
  LOGICAL                                :: UseForEmission                   ! Use Init/Emission for emission?
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  INTEGER(8)                             :: initialParticleNumber            ! Number of Particles at time 0.0
  REAL                                   :: RadiusIC                         ! Radius for IC circle
  REAL                                   :: Radius2IC                        ! Radius2 for IC cylinder (ring)
  REAL                                   :: RadiusICGyro                     ! Radius for Gyrotron gyro radius
  REAL                                   :: NormalIC(3)                      ! Normal / Orientation of circle
  REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                                   :: CuboidHeightIC                   ! third measure of cuboid
                                                                             ! (set 0 for flat rectangle),
                                                                             ! negative value = opposite direction
  REAL                                   :: CylinderHeightIC                 ! third measure of cylinder
                                                                             ! (set 0 for flat rectangle),
                                                                             ! negative value = opposite direction
  LOGICAL                                :: CalcHeightFromDt                 ! Calc. cuboid/cylinder height from v and dt?
  REAL                                   :: VeloIC                           ! velocity for inital Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: Amplitude                        ! Amplitude for sin-deviation initiation.
  REAL                                   :: WaveNumber                       ! WaveNumber for sin-deviation initiation.
  INTEGER                                :: maxParticleNumberX               ! Maximum Number of all Particles in x direction
  INTEGER                                :: maxParticleNumberY               ! Maximum Number of all Particles in y direction
  INTEGER                                :: maxParticleNumberZ               ! Maximum Number of all Particles in z direction 
  REAL                                   :: WeibelVeloPar                    ! Parrallel velocity component for Weibel
  REAL                                   :: WeibelVeloPer                    ! Perpendicular velocity component for Weibel
  REAL                                   :: OneDTwoStreamVelo                ! Stream Velocity for the Two Stream Instability
  REAL                                   :: OneDTwoStreamTransRatio          ! Ratio between perpendicular and parallel velocity
  REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: ConstantPressure                 ! Pressure for an Area with a Constant Pressure
  REAL                                   :: ConstPressureRelaxFac            ! RelaxFac. for ConstPressureSamp
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3) for LD_insert or
                                                                             ! (vpi_)cub./cyl. as alternative to Part.Emis. in Type1
  INTEGER                                :: ParticleEmissionType             ! Emission Type 1 = emission rate in 1/s,
                                                                             !               2 = emission rate 1/iteration
                                                                             !               3 = user def. emission rate
                                                                             !               4 = const. cell pressure
                                                                             !               5 = cell pres. w. complete part removal
                                                                             !               6 = outflow BC (characteristics method)
  REAL                                   :: ParticleEmission                 ! Emission in [1/s] or [1/Iteration]
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  REAL                                   :: Nsigma                           ! sigma multiple of maxwell for virtual insert length
  LOGICAL                                :: VirtPreInsert                    ! virtual Pre-Inserting region (adapted SetPos/Velo)?
  CHARACTER(40)                          :: vpiDomainType                    ! specifying Keyword for virtual Pre-Inserting region
                                                                             ! implemented: - perpendicular_extrusion (default)
                                                                             !              - freestream
                                                                             !              - orifice
                                                                             !              - ...more following...
  LOGICAL                                :: vpiBVBuffer(4)                   ! incl. buffer region in -BV1/+BV1/-BV2/+BV2 direction?
  TYPE(tConstPressure)                   :: ConstPress!(:)           =>NULL() !
  INTEGER                                :: NumberOfExcludeRegions           ! Number of different regions to be excluded
  TYPE(tExcludeRegion), ALLOCATABLE      :: ExcludeRegion(:)
#ifdef MPI
  INTEGER                                :: InitComm                          ! number of init-communicator
#endif /*MPI*/
END TYPE tInit

TYPE tSpecies                                                                ! Particle Data for each Species
  !General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  REAL                                   :: ChargeIC                         ! Particle Charge (without MPF)
  REAL                                   :: MassIC                           ! Particle Mass (without MPF)
  REAL                                   :: MacroParticleFactor              ! Number of Microparticle per Macroparticle
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  INTEGER                                :: StartnumberOfInits               ! 0 if old emit defined (array is copied into 0. entry)
END TYPE

INTEGER                                  :: nSpecies                         ! number of species
TYPE(tSpecies), ALLOCATABLE              :: Species(:)  !           => NULL() ! Species Data Vector

TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleAnodeBC           = 4      ! = 4 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleCathodeBC         = 5      ! = 5 (s.u.) Boundary Condition Integer Definition
  !INTEGER                                :: MPINeighborhoodBC       = 6      ! = 6 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)           , POINTER :: SourceBoundName(:) =>NULL() ! Link part 1 for mapping Boltzplatz BCs to Particle BC
  INTEGER                      , POINTER :: TargetBoundCond(:) =>NULL() ! Link part 2 for mapping Boltzplatz BCs to Particle BC
  INTEGER                      , POINTER :: Map(:)             =>NULL() ! Map from Boltzplatz BCindex to Particle BC
  INTEGER                      , POINTER :: MapToPartBC(:)     =>NULL() ! Map from Boltzplatz BCindex to Particle BC (NOT TO TYPE!)
  REAL    , ALLOCATABLE                  :: MomentumACC(:)      
  REAL    , ALLOCATABLE                  :: WallTemp(:)     
  REAL    , ALLOCATABLE                  :: TransACC(:)     
  REAL    , ALLOCATABLE                  :: VibACC(:) 
  REAL    , ALLOCATABLE                  :: RotACC(:) 
  REAL    , ALLOCATABLE                  :: WallVelo(:,:) 
  REAL    , ALLOCATABLE                  :: Voltage(:)
  INTEGER , ALLOCATABLE                  :: NbrOfSpeciesSwaps(:)          !Number of Species to be changed at wall
  REAL    , ALLOCATABLE                  :: ProbOfSpeciesSwaps(:)         !Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE                  :: SpeciesSwaps(:,:,:)           !Species to be changed at wall (in, out), out=0: delete
  LOGICAL , ALLOCATABLE                  :: AmbientCondition(:)
  REAL    , ALLOCATABLE                  :: AmbientTemp(:)
  REAL    , ALLOCATABLE                  :: AmbientMeanPartMass(:)
  REAL    , ALLOCATABLE                  :: AmbientBeta(:)
  REAL    , ALLOCATABLE                  :: AmbientVelo(:,:)
  REAL    , ALLOCATABLE                  :: AmbientDens(:)
  REAL    , ALLOCATABLE                  :: AmbientDynamicVisc(:)               ! dynamic viscousity
  REAL    , ALLOCATABLE                  :: AmbientThermalCond(:)               ! thermal conuctivity

END TYPE

INTEGER                                  :: nPartBound                       ! number of particle boundaries
TYPE(tPartBoundary)                      :: PartBound                         ! Boundary Data for Particles

TYPE tParticleElementMapping
  INTEGER                      , POINTER :: Element(:)             =>NULL()  ! Element number allocated to each Particle
  INTEGER                      , POINTER :: lastElement(:)         =>NULL()  ! Element number allocated
                                                                             ! to each Particle at previous timestep
!----------------------------------------------------------------------------!----------------------------------
                                                                             ! Following vectors are assigned in
                                                                             ! SUBROUTINE UpdateNextFreePosition
                                                                             ! IF (PIC%withDSMC .OR. PIC%withFP)
  INTEGER                      , POINTER :: pStart(:)              =>NULL()  ! Start of Linked List for Particles in Element
                                                                             ! pStart(1:PIC%nElem)
  INTEGER                      , POINTER :: pNumber(:)             =>NULL()  ! Number of Particles in Element
                                                                             ! pStart(1:PIC%nElem)
  INTEGER                      , POINTER :: pEnd(:)                =>NULL()  ! End of Linked List for Particles in Element
                                                                             ! pEnd(1:PIC%nElem)
  INTEGER                      , POINTER :: pNext(:)               =>NULL()  ! Next Particle in same Element (Linked List)
                                                                             ! pStart(1:PIC%maxParticleNumber)
END TYPE

TYPE(tParticleElementMapping)            :: PEM

TYPE tParticleDataManagement
  INTEGER                                :: CurrentNextFreePosition           ! Index of nextfree index in nextFreePosition-Array
  INTEGER                                :: maxParticleNumber                 ! Maximum Number of all Particles
  INTEGER                                :: ParticleVecLength                 ! Vector Length for Particle Push Calculation
  INTEGER                                :: insideParticleNumber              ! Number of all recent Particles inside
  INTEGER , ALLOCATABLE                  :: PartInit(:)                       ! (1:NParts), initial emission condition number
                                                                              ! the calculation area
  INTEGER                      , POINTER :: nextFreePosition(:)    =>NULL()   ! next_free_Position(1:max_Particle_Number)
                                                                              ! List of free Positon
  LOGICAL                      , POINTER :: ParticleInside(:)      =>NULL()   ! Particle_inside(1:Particle_Number)
END TYPE

TYPE (tParticleDataManagement)           :: PDM

TYPE tFastInitBGM
  INTEGER                                :: nElem                             ! Number of elements in background mesh cell
  INTEGER, ALLOCATABLE                   :: Element(:)                        ! List of elements in BGM cell
#ifdef MPI     
  INTEGER, ALLOCATABLE                   :: ShapeProcs(:)                     ! first Entry: Number of Shapeprocs, 
                                                                              ! following: ShapeProcs
  INTEGER, ALLOCATABLE                   :: PaddingProcs(:)                   ! first Entry: Number of Paddingprocs, 
                                                                              ! following: PaddingProcs
  INTEGER, ALLOCATABLE                   :: SharedProcs(:)                    ! first Entry: Number of Sharedprocs, 
                                                                              ! following: SharedProcs
  INTEGER                                :: nBCSides                          ! number BC sides in BGM cell
#endif                     
END TYPE

TYPE tGeometry
  INTEGER, ALLOCATABLE                   :: ElemToNodeID(:,:)                 ! ElemToNodeID(1:nElemNodes,1:nElems)
  INTEGER, ALLOCATABLE                   :: ElemSideNodeID(:,:,:)             ! ElemSideNodeID(1:nSideNodes,1:nLocSides,1:nElems)
                                                                              ! From element sides to node IDs
  INTEGER, ALLOCATABLE                   :: PeriodicElemSide(:,:)             ! 0=not periodic side, others=PeriodicVectorsNum
  INTEGER, ALLOCATABLE                   :: PeriodicBGMVectors(:,:)           ! = periodic vectors in backgroundmesh coords
  LOGICAL, ALLOCATABLE                   :: ConcaveElemSide(:,:)              ! Whether LocalSide of Element is concave side
  REAL, ALLOCATABLE                      :: NodeCoords(:,:)                   ! Node Coordinates (1:nDim,1:nNodes)
  REAL, ALLOCATABLE                      :: Volume(:)                         ! Volume(nElems) for nearest_blurrycenter
  REAL, ALLOCATABLE                      :: DeltaEvMPF(:)                     ! Energy difference due to particle merge
  REAL, ALLOCATABLE                      :: PeriodicVectors(:,:)              ! PeriodicVectors(1:3,1:nPeriodicVectors), 1:3=x,y,z
  TYPE (tFastInitBGM)          , POINTER :: FIBGM(:,:,:)           =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: FIBGMimin                         ! smallest index of FastInitBGM (x)
  INTEGER                                :: FIBGMimax                         ! biggest index of FastInitBGM (x)
  INTEGER                                :: FIBGMkmin                         ! smallest index of FastInitBGM (y)
  INTEGER                                :: FIBGMkmax                         ! biggest index of FastInitBGM (y)
  INTEGER                                :: FIBGMlmin                         ! smallest index of FastInitBGM (z)
  INTEGER                                :: FIBGMlmax                         ! biggest index of FastInitBGM (z)
  INTEGER                                :: nPeriodicVectors                  ! Number of periodic Vectors
  REAL                                   :: FIBGMdeltas(3)                    ! size of background mesh cell for particle init
  REAL                                   :: FactorFIBGM(3)                    ! scaling factor for FIBGM
  REAL                                   :: xminglob                          ! global minimum x coord of all nodes
  REAL                                   :: yminglob                          ! global minimum y coord of all nodes
  REAL                                   :: zminglob                          ! global minimum z coord of all nodes
  REAL                                   :: xmaxglob                          ! global max x coord of all nodes
  REAL                                   :: ymaxglob                          ! global max y coord of all nodes
  REAL                                   :: zmaxglob                          ! global max z coord of all nodes
  REAL                                   :: xmin                              ! minimum x coord of all nodes
  REAL                                   :: xmax                              ! maximum x coord of all nodes
  REAL                                   :: ymin                              ! minimum y coord of all nodes
  REAL                                   :: ymax                              ! maximum y coord of all nodes
  REAL                                   :: zmin                              ! minimum z coord of all nodes
  REAL                                   :: zmax                              ! maximum z coord of all nodes
  REAL                                   :: nnx(3),nny(3),nnz(3)              ! periodic vectors
  LOGICAL                                :: SelfPeriodic                      ! does process have periodic bounds with itself?
END TYPE

TYPE (tGeometry)                         :: GEO

REAL                                     :: DelayTime
REAL                                     :: Time                              ! Simulation Time

LOGICAL           :: ParticlesInitIsDone=.FALSE.

LOGICAL                                  :: WriteMacroValues                  ! Output of macroscopic values
INTEGER                                  :: MacroValSamplIterNum              ! Number of iterations for sampling   
                                                                              ! macroscopic values
REAL                                     :: MacroValSampTime                  ! Sampling time for WriteMacroVal. (e.g., for td201)
LOGICAL                                  :: usevMPF                           ! use the vMPF per particle
LOGICAL                                  :: enableParticleMerge               ! enables the particle merge routines
LOGICAL                                  :: doParticleMerge=.false.           ! flag for particle merge
INTEGER                                  :: vMPFMergeParticleTarget           ! number of particles wanted after merge
INTEGER                                  :: vMPFSplitParticleTarget           ! number of particles wanted after split
INTEGER                                  :: vMPFMergeParticleIter             ! iterations between particle merges
INTEGER                                  :: vMPFMergePolyOrder                ! order of polynom for vMPF merge
INTEGER                                  :: vMPFMergeCellSplitOrder           ! order of cell splitting (vMPF)
INTEGER, ALLOCATABLE                     :: vMPF_OrderVec(:,:)                ! Vec of vMPF poynom orders
INTEGER, ALLOCATABLE                     :: vMPF_SplitVec(:,:)                ! Vec of vMPF cell split orders
INTEGER, ALLOCATABLE                     :: vMPF_SplitVecBack(:,:,:)          ! Vec of vMPF cell split orders backward
REAL, ALLOCATABLE                        :: PartStateMap(:,:)                 ! part pos mapped on the -1,1 cube  
INTEGER, ALLOCATABLE                     :: PartStatevMPFSpec(:)              ! part state indx of spec to merge
REAL, ALLOCATABLE                        :: vMPFPolyPoint(:,:)                ! Points of Polynom in LM 
REAL, ALLOCATABLE                        :: vMPFPolySol(:)                    ! Solution of Polynom in LM
REAL                                     :: vMPF_oldMPFSum                    ! Sum of all old MPF in cell
REAL                                     :: vMPF_oldEngSum                    ! Sum of all old energies in cell
REAL                                     :: vMPF_oldMomSum(3)                 ! Sum of all old momentums in cell
REAL, ALLOCATABLE                        :: vMPFOldVelo(:,:)                  ! Old Particle Velo for Polynom
REAL, ALLOCATABLE                        :: vMPFOldBrownVelo(:,:)             ! Old brownian Velo
REAL, ALLOCATABLE                        :: vMPFOldPos(:,:)                   ! Old Particle Pos for Polynom
REAL, ALLOCATABLE                        :: vMPFOldMPF(:)                     ! Old Particle MPF
REAL, ALLOCATABLE                        :: vMPFNewPosNum(:)
INTEGER, ALLOCATABLE                     :: vMPF_SpecNumElem(:,:)             ! number of particles of spec (:,i) in element (j,:)
CHARACTER(30)                            :: vMPF_velocityDistribution         ! specifying keyword for velocity distribution
REAL, ALLOCATABLE                        :: vMPF_NewPosRefElem(:,:)          ! new positions in ref elem
LOGICAL                                  :: vMPF_relativistic
LOGICAL                                  :: PartPressureCell                  ! Flag: constant pressure in cells emission (type4)
LOGICAL                                  :: PartPressAddParts                 ! Should Parts be added to reach wanted pressure?
LOGICAL                                  :: PartPressRemParts                 ! Should Parts be removed to reach wanted pressure?
INTEGER                                  :: NumRanVec      ! Number of predefined random vectors
REAL  , ALLOCATABLE                      :: RandomVec(:,:) ! Random Vectos (NumRanVec, direction)
REAL, ALLOCATABLE                        :: RegionElectronRef(:,:)          ! RegionElectronRef((rho0,phi0,Te[eV])|1:NbrOfRegions)
LOGICAL                                  :: useVTKFileBGG                     ! Flag for BGG via VTK-File
REAL, ALLOCATABLE                        :: BGGdataAtElem(:,:)                ! data for BGG via VTK-File
LOGICAL                                  :: OutputVpiWarnings                 ! Flag for warnings for rejected v if VPI+PartDensity

INTEGER(8)                               :: nTotalPart
INTEGER(8)                               :: nTotalHalfPart

!===================================================================================================================================
END MODULE MOD_Particle_Vars
