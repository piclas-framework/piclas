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
REAL                  :: ManualTimeStep                                      ! Manual TimeStep
LOGICAL               :: useManualTimeStep                                   ! Logical Flag for manual timestep. For consistency
                                                                             ! with IAG programming style
LOGICAL               :: Symmetry2D                                          ! Enables 2D simulation: symmetry in xy-Plane
LOGICAL               :: Symmetry2DAxisymmetric                              ! Enables axisymmetric simulation around z-axis
LOGICAL               :: DoFieldIonization                                   ! Do Field Ionization by quantum tunneling
INTEGER               :: FieldIonizationModel                                !'Field Ionization models. Implemented models are:
!                                                                            ! * Ammosov-Delone-Krainov (ADK) model
!                                                                            ! * Ammosov-Delone-Krainov (ADK) model Yu 2018
LOGICAL,ALLOCATABLE   :: SpecReset(:)                                        ! Flag for resetting species distribution with init
                                                                             ! during restart
LOGICAL               :: KeepWallParticles                                   ! Flag for tracking of adsorbed Particles
LOGICAL               :: SolidSimFlag                                        ! Flag telling if Solid boundary is existing
LOGICAL               :: LiquidSimFlag                                       ! Flag telling if Liquid boundary is existing
INTEGER               :: PartSurfaceModel                                    ! Model used for wall interaction
                                                                             ! 0 perfect/diffusive reflection
                                                                             ! 1 adsorption (Kisluik) / desorption (Polanyi Wigner)
                                                                             ! 2 Recombination coefficient (Laux model)
                                                                             ! 3 adsorption/desorption + chemical interaction
                                                                             !   (SMCR with UBI-QEP, TST and TCE)
                                                                             ! 4 TODO
                                                                             ! 5 SEE (secondary e- emission) by Levko2015
                                                                             ! 6 SEE (secondary e- emission) by Pagonakis2016
                                                                             !   (orignally from Harrower1956)
LOGICAL               :: printRandomSeeds                                    ! print random seeds or not
! IMD: Molecular Dynamics Model - ion distribution info
LOGICAL               :: DoInitialIonization                                 ! When restarting from a state, ionize the species to a
                                                                             ! specific degree
INTEGER               :: InitialIonizationSpecies                            ! Supply the number of species that are considered for
                                                                             ! automatic ionization
INTEGER , ALLOCATABLE :: InitialIonizationSpeciesID(:)                       ! Supply a vector with the species IDs that are used
                                                                             ! for the initial ionization.
REAL                  :: InitialIonizationChargeAverage                      ! Average charge for each atom/molecule in the cell
LOGICAL               :: DoImportIMDFile                                     ! read IMD (MD-Simulation) data from *.chkpt file
REAL                  :: IMDTimeScale                                        ! Time unit of input file
REAL                  :: IMDLengthScale                                      ! global IMD length scale
INTEGER               :: IMDNumber                                           ! Output number IMD Data file
CHARACTER(255)        :: IMDInputFile                                        ! Laser data file name containing PartState(1:6)
INTEGER               :: IMDnSpecies                                         ! number of IMD species
INTEGER , ALLOCATABLE :: IMDSpeciesID(:)                                     ! species ID for distributing the IMD atoms/ions
INTEGER , ALLOCATABLE :: IMDSpeciesCharge(:)                                 ! charge number of IMD atoms/ions
CHARACTER(255)        :: IMDAtomFile                                         ! Laser data file name containing PartState(1:6)
REAL                  :: IMDCutOffxValue                                     ! cut-off coordinate for IMDCutOff='coordiantes'
CHARACTER(255)        :: IMDCutOff                                           ! cut-off type for IMD data reduction: 1.) no_cutoff
                                                                             !                                      2.) Epot
                                                                             !                                      3.) coordinates
                                                                             !                                      4.) velocity
REAL                  :: dt_max_particles                                    ! Maximum timestep for particles (for static fields!)
REAL                  :: dt_maxwell                                          ! timestep for field solver (for static fields only!)
REAL                  :: dt_adapt_maxwell                                    ! adapted timestep for field solver dependent
                                                                             ! on particle velocity (for static fields only!)
REAL                  :: dt_part_ratio, overrelax_factor                     ! factors for td200/201 overrelaxation/subcycling
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
LOGICAL               :: DoForceFreeSurfaceFlux                              ! switch if the stage reconstruction uses a force
#if (PP_TimeDiscMethod==509)
LOGICAL               :: velocityOutputAtTime
REAL    , ALLOCATABLE :: velocityAtTime(:,:)
#endif /*(PP_TimeDiscMethod==509)*/
#if defined(ROS) || defined(IMPA)
REAL    , ALLOCATABLE :: PartStage (:,:,:)                                   ! ERK4 additional function values
REAL    , ALLOCATABLE :: PartStateN(:,:)                                     ! ParticleState at t^n
REAL    , ALLOCATABLE :: PartdtFrac(:)                                       ! dual use variable:
REAL    , ALLOCATABLE :: PartQ(:,:)                                          ! PartilceState at t^n or state at RK-level 0
                                                                             ! 1) time fraction of domain entering (surface flux)
                                                                             ! 2) fraction of time step for push (surface flux)
#endif /*IMPA || ROS*/
#if defined(IMPA)
LOGICAL , ALLOCATABLE :: PartIsImplicit(:)                                   ! select, if specific particle is explicit or implicit
REAL    , ALLOCATABLE :: PartDeltaX(:,:)                                     ! Change of particle during Newton step
LOGICAL , ALLOCATABLE :: PartLambdaAccept(:)                                 ! Accept particle search direction
! Newton iteration
REAL    , ALLOCATABLE :: F_PartX0(:,:)                                       ! Particle function evaluated at t^0
REAL    , ALLOCATABLE :: F_PartXK(:,:)                                       ! Particle function evaluated at iteration step k
REAL    , ALLOCATABLE :: Norm_F_PartX0    (:)                               ! and the corresponding L2 norm
REAL    , ALLOCATABLE :: Norm_F_PartXK    (:)                               ! and the corresponding L2 norm
REAL    , ALLOCATABLE :: Norm_F_PartXK_Old(:)                               ! and the corresponding L2 norm
LOGICAL , ALLOCATABLE :: DoPartInNewton(:)                                   ! particle is treated implicitly && Newtons method
                                                                             ! is performed on it
#endif
REAL    , ALLOCATABLE :: Pt_temp(:,:)                                        ! LSERK4 additional derivative of PartState

                                                                             ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE :: LastPartPos(:,:)                                    ! (1:NParts,1:3) with 2nd index: x,y,z
INTEGER , ALLOCATABLE :: PartSpecies(:)                                      ! (1:NParts)
REAL    , ALLOCATABLE :: PartMPF(:)                                          ! (1:NParts) MacroParticleFactor by variable MPF
INTEGER               :: PartLorentzType
CHARACTER(LEN=256)    :: ParticlePushMethod                                  ! Type of PP-Method
INTEGER               :: nrSeeds                                             ! Number of Seeds for Random Number Generator
INTEGER , ALLOCATABLE :: seeds(:)                        !        =>NULL()   ! Seeds for Random Number Generator

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
  INTEGER                                :: Rotation                         ! direction of rotation, similar to TE-mode
  INTEGER                                :: VelocitySpreadMethod             ! method to compute the velocity spread
  REAL                                   :: InflowRiseTime                   ! time to ramp the number of inflow particles
                                                                             ! linearly from zero to unity
  REAL                                   :: VelocitySpread                   ! velocity spread in percent
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
  REAL                                   :: MJxRatio                         ! x direction portion of velocity for Maxwell-Juettner
  REAL                                   :: MJyRatio                         ! y direction portion of velocity for Maxwell-Juettner
  REAL                                   :: MJzRatio                         ! z direction portion of velocity for Maxwell-Juettner
  REAL                                   :: WeibelVeloPar                    ! Parallel velocity component for Weibel
  REAL                                   :: WeibelVeloPer                    ! Perpendicular velocity component for Weibel
  REAL                                   :: OneDTwoStreamVelo                ! Stream Velocity for the Two Stream Instability
  REAL                                   :: OneDTwoStreamTransRatio          ! Ratio between perpendicular and parallel velocity
  REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: ConstantPressure                 ! Pressure for an Area with a Constant Pressure
  REAL                                   :: ConstPressureRelaxFac            ! RelaxFac. for ConstPressureSamp
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3) for LD_insert or
                                                                             ! (vpi_)cub./cyl. as alternative to Part.Emis. in Type1
  INTEGER                                :: ElemTemperatureFileID
  REAL , ALLOCATABLE                     :: ElemTemperatureIC(:,:)           ! Temperature from macrorestart [1:3,1:nElems)
  INTEGER                                :: ElemPartDensityFileID
  REAL , ALLOCATABLE                     :: ElemPartDensity(:)
  INTEGER                                :: ElemVelocityICFileID
  REAL , ALLOCATABLE                     :: ElemVelocityIC(:,:)
  INTEGER                                :: ElemTVibFileID
  REAL , ALLOCATABLE                     :: ElemTVib(:)                      ! vibrational temperature [nElems]
  INTEGER                                :: ElemTRotFileID
  REAL , ALLOCATABLE                     :: ElemTRot(:)                      ! rotational temperature [nElems]
  INTEGER                                :: ElemTelecFileID
  REAL , ALLOCATABLE                     :: ElemTelec(:)                     ! electronic temperature [nElems]
  INTEGER                                :: ParticleEmissionType             ! Emission Type 1 = emission rate in 1/s,
                                                                             !               2 = emission rate 1/iteration
                                                                             !               3 = user def. emission rate
                                                                             !               4 = const. cell pressure
                                                                             !               5 = cell pres. w. complete part removal
                                                                             !               6 = outflow BC (characteristics method)
  REAL                                   :: ParticleEmission                 ! Emission in [1/s] or [1/Iteration]
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=4)                        :: InsertedParticleMisMatch=0       ! error in number of inserted particles of last step
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
#if USE_MPI
  INTEGER                                :: InitComm                          ! number of init-communicator
#endif /*USE_MPI*/
END TYPE tInit

TYPE tSurfFluxSubSideData
  REAL                                   :: projFak                          ! VeloVecIC projected to inwards normal
  REAL                                   :: a_nIn                            ! speed ratio projected to inwards normal
  REAL                                   :: Velo_t1                          ! Velo comp. of first orth. vector
  REAL                                   :: Velo_t2                          ! Velo comp. of second orth. vector
  REAL                                   :: nVFR                             ! normal volume flow rate through subside
  REAL                                   :: Dmax                             ! maximum Jacobian determinant of subside for opt. ARM
  REAL,ALLOCATABLE                       :: BezierControlPoints2D(:,:,:)     ! BCP of SubSide projected to VeloVecIC
                                                                             ! (1:2,0:NGeo,0:NGeo)
END TYPE tSurfFluxSubSideData

TYPE tSurfFluxLink
  INTEGER                                :: PartIdx
  INTEGER,ALLOCATABLE                    :: SideInfo(:)
  TYPE(tSurfFluxLink), POINTER           :: next => null()
END TYPE tSurfFluxLink

TYPE typeSurfaceflux
  INTEGER                                :: BC                               ! PartBound to be emitted from
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  REAL                                   :: VeloIC                           ! velocity for inital Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  LOGICAL                                :: VeloIsNormal                     ! VeloIC is in Surf-Normal instead of VeloVecIC
  LOGICAL                                :: ReduceNoise                      ! reduce stat. noise by global calc. of PartIns
  LOGICAL                                :: AcceptReject                     ! perform ARM for skewness of RefMap-positioning
  INTEGER                                :: ARM_DmaxSampleN                  ! number of sample intervals in xi/eta for Dmax-calc.
  REAL                                   :: VFR_total                        ! Total Volumetric flow rate through surface
  REAL                     , ALLOCATABLE :: VFR_total_allProcs(:)            ! -''-, all values for root in ReduceNoise-case
  REAL                                   :: VFR_total_allProcsTotal          !     -''-, total
  REAL                                   :: totalAreaSF                      ! Total area of the respective surface flux
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=8)                        :: tmpInsertedParticle              ! tmp Number of all already inserted Particles
  INTEGER(KIND=8)                        :: tmpInsertedParticleSurplus       ! tmp accumulated "negative" number of inserted Particles
  TYPE(tSurfFluxSubSideData), ALLOCATABLE :: SurfFluxSubSideData(:,:,:)      ! SF-specific Data of Sides (1:N,1:N,1:SideNumber)
  LOGICAL                                :: SimpleRadialVeloFit              ! fit of veloR/veloTot=-r*(A*exp(B*r)+C)
  REAL                                   :: preFac !A
  REAL                                   :: powerFac !B
  REAL                                   :: shiftFac !C
  LOGICAL                                :: CircularInflow                   ! Circular region, which can be used to define small
                                                                             ! geometry features on large boundaries
  INTEGER                                :: dir(3)                           ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                                   :: origin(2)                        ! origin in orth. coordinates of polar system
  REAL                                   :: rmax                             ! max radius of to-be inserted particles
  REAL                                   :: rmin                             ! min radius of to-be inserted particles
  INTEGER, ALLOCATABLE                   :: SurfFluxSideRejectType(:)        ! Type if parts in side can be rejected (1:SideNumber)
  REAL                                   :: PressureFraction
  TYPE(tSurfFluxLink), POINTER           :: firstSurfFluxPart => null()      ! pointer to first particle inserted for iSurfaceFlux
                                                                             ! used for linked list during sampling
  TYPE(tSurfFluxLink), POINTER           :: lastSurfFluxPart => null()       ! pointer to last particle inserted for iSurfaceFlux
                                                                             ! used for abort criterion in do while during sampling
  LOGICAL                                :: Adaptive                         ! Is the surface flux an adaptive boundary?
  INTEGER                                :: AdaptiveType                     ! Chose the adaptive type, description in DefineParams
  REAL                                   :: AdaptiveMassflow                 ! Mass flow [kg/s], which is held constant
  REAL                                   :: AdaptivePressure                 ! Static pressure [Pa], which is held constant
  REAL, ALLOCATABLE                      :: AdaptivePreviousVelocity(:,:)    ! A velocity is stored in case of negative values
  REAL, ALLOCATABLE                      :: ConstMassflowWeight(:,:,:)       ! Adaptive, Type 4: Weighting factor for SF-sides to
                                                                             ! insert the right amount of particles
  REAL, ALLOCATABLE                      :: CircleAreaPerTriaSide(:,:,:)     ! Adaptive, Type 4: Area within a triangle, determined
                                                                             ! through Monte Carlo integration (initially)
  INTEGER                                :: AdaptivePartNumOut               ! Adaptive, Type 4: Number of particles exiting through
                                                                             ! the adaptive boundary condition
END TYPE

TYPE tSpecies                                                                ! Particle Data for each Species
  !General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  REAL                                   :: ChargeIC                         ! Particle Charge (without MPF)
  REAL                                   :: MassIC                           ! Particle Mass (without MPF)
  REAL                                   :: MacroParticleFactor              ! Number of Microparticle per Macroparticle
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  INTEGER                                :: StartnumberOfInits               ! 0 if old emit defined (array is copied into 0. entry)
  TYPE(typeSurfaceflux),ALLOCATABLE      :: Surfaceflux(:)                   ! Particle Data for each SurfaceFlux emission
  INTEGER                                :: nSurfacefluxBCs                  ! Number of SF emissions
#if IMPA
  LOGICAL                                :: IsImplicit
#endif
END TYPE

LOGICAL                                  :: UseCircularInflow                !
LOGICAL                                  :: UseAdaptive                 !
REAL                                     :: AdaptiveWeightFac                ! weighting factor theta for weighting of average
                                                                             ! instantaneous values with those
                                                                             ! of previous iterations
REAL, ALLOCATABLE                        :: Adaptive_MacroVal(:,:,:)         ! Macroscopic value near boundaries
                                                                             ! (1:14,1:nElems,1:nSpecies)
                                                                             !  1:  VELOX
                                                                             !  2:  VELOY
                                                                             !  3:  VELOZ
                                                                             !  4:  TEMPX
                                                                             !  5:  TEMPY
                                                                             !  6:  TEMPZ
                                                                             !  7:  NUMBER DENSITY
                                                                             !  8:  TVIB
                                                                             !  9:  TROT
                                                                             ! 10:  TELEC
                                                                             ! 11:  Pumping capacity [m3/s]
                                                                             ! 12:  Static pressure [Pa]
                                                                             ! 13:  Integral pressure difference [Pa]
REAL,ALLOCATABLE                         :: MacroRestartData_tmp(:,:,:,:)    ! Array of macrovalues read from macrorestartfile

INTEGER                                  :: nSpecies                         ! number of species
INTEGER                                  :: nMacroRestartFiles                ! number of macroscopic restart files used for particles
TYPE(tSpecies), ALLOCATABLE              :: Species(:)  !           => NULL() ! Species Data Vector

LOGICAL                                  :: PartMeshHasPeriodicBCs
#if defined(IMPA) || defined(ROS)
LOGICAL                                  :: PartMeshHasReflectiveBCs
#endif
TYPE tParticleElementMapping
  INTEGER                , ALLOCATABLE   :: Element(:)      !      =>NULL()  ! Element number allocated to each Particle
  INTEGER                , ALLOCATABLE   :: lastElement(:)  !      =>NULL()  ! Element number allocated
#if defined(IMPA) || defined(ROS)
  INTEGER                , ALLOCATABLE   :: ElementN(:)  !      =>NULL()  ! Element number allocated
  REAL                   , ALLOCATABLE   :: NormVec(:,:)  !      =>NULL()  ! Element number allocated
  LOGICAL                , ALLOCATABLE   :: PeriodicMoved(:)                 ! flag, if the particle moved over periodic bcs
#endif
                                                                             ! to each Particle at previous timestep
!----------------------------------------------------------------------------!----------------------------------
                                                                             ! Following vectors are assigned in
                                                                             ! SUBROUTINE UpdateNextFreePosition
                                                                             ! IF (PIC%withDSMC .OR. PIC%withFP)
  INTEGER                , ALLOCATABLE    :: pStart(:)         !     =>NULL()  ! Start of Linked List for Particles in Element
                                                               !               ! pStart(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pNumber(:)        !     =>NULL()  ! Number of Particles in Element
                                                               !               ! pStart(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: wNumber(:)        !     =>NULL()  ! Number of Wall-Particles in Element
                                                                               ! pStart(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pEnd(:)           !     =>NULL()  ! End of Linked List for Particles in Element
                                                               !               ! pEnd(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pNext(:)          !     =>NULL()  ! Next Particle in same Element (Linked List)
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
  INTEGER ,ALLOCATABLE                   :: nextFreePosition(:)  !  =>NULL()  ! next_free_Position(1:max_Particle_Number)
                                                                              ! List of free Positon
  LOGICAL ,ALLOCATABLE                   :: ParticleInside(:)    !  =>NULL()  ! Particle_inside(1:Particle_Number)
  LOGICAL , ALLOCATABLE                  :: ParticleAtWall(:)                 ! Particle_adsorbed_on_to_wall(1:Particle_number)
  INTEGER , ALLOCATABLE                  :: PartAdsorbSideIndx(:,:)           ! Surface index on which Particle i adsorbed
                                                                              ! (1:3,1:PDM%maxParticleNumber)
                                                                              ! 1: surface index ElemToSide(i,localsideID,ElementID)
                                                                              ! 2: p
                                                                              ! 3: q
  LOGICAL ,ALLOCATABLE                   :: dtFracPush(:)                     ! Push random fraction only
  LOGICAL ,ALLOCATABLE                   :: IsNewPart(:)                      ! Reconstruct RK-scheme in next stage
END TYPE

TYPE (tParticleDataManagement)           :: PDM

REAL                                     :: DelayTime

LOGICAL                                  :: ParticlesInitIsDone=.FALSE.

LOGICAL                                  :: WRITEMacroValues = .FALSE.
LOGICAL                                  :: WriteMacroVolumeValues =.FALSE.   ! Output of macroscopic values in volume
LOGICAL                                  :: WriteMacroSurfaceValues=.FALSE.   ! Output of macroscopic values on surface
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
INTEGER, ALLOCATABLE                     :: vMPFNewPosNum(:)
INTEGER, ALLOCATABLE                     :: vMPF_SpecNumElem(:,:)             ! number of particles of spec (:,i) in element (j,:)
CHARACTER(30)                            :: vMPF_velocityDistribution         ! specifying keyword for velocity distribution
REAL, ALLOCATABLE                        :: vMPF_NewPosRefElem(:,:)          ! new positions in ref elem
LOGICAL                                  :: vMPF_relativistic
LOGICAL                                  :: PartPressureCell                  ! Flag: constant pressure in cells emission (type4)
LOGICAL                                  :: PartPressAddParts                 ! Should Parts be added to reach wanted pressure?
LOGICAL                                  :: PartPressRemParts                 ! Should Parts be removed to reach wanted pressure?
REAL, ALLOCATABLE                        :: RegionElectronRef(:,:)          ! RegionElectronRef((rho0,phi0,Te[eV])|1:NbrOfRegions)
LOGICAL                                  :: OutputVpiWarnings                 ! Flag for warnings for rejected v if VPI+PartDensity
LOGICAL                                  :: DoSurfaceFlux                     ! Flag for emitting by SurfaceFluxBCs
LOGICAL                                  :: DoPoissonRounding                 ! Perform Poisson sampling instead of random rounding
LOGICAL                                  :: DoTimeDepInflow                   ! Insertion and SurfaceFlux w simple random rounding
LOGICAL                                  :: DoZigguratSampling                ! Sample normal randoms with Ziggurat method

INTEGER(8)                               :: nTotalPart
INTEGER(8)                               :: nTotalHalfPart

INTEGER :: nCollectChargesBCs
INTEGER :: nDataBC_CollectCharges
TYPE tCollectCharges
  INTEGER                              :: BC
  REAL                                 :: NumOfRealCharges
  REAL                                 :: NumOfNewRealCharges
  REAL                                 :: ChargeDist
END TYPE
TYPE(tCollectCharges), ALLOCATABLE     :: CollectCharges(:)

TYPE tVariableTimeStep
  LOGICAL                              :: UseVariableTimeStep
  LOGICAL                              :: UseLinearScaling
  LOGICAL                              :: UseDistribution
  REAL, ALLOCATABLE                    :: ParticleTimeStep(:)
  REAL, ALLOCATABLE                    :: ElemFac(:)
  REAL, ALLOCATABLE                    :: ElemWeight(:)
  REAL                                 :: StartPoint(3)
  REAL                                 :: EndPoint(3)
  REAL                                 :: Direction(3)
  REAL                                 :: ScaleFac
  LOGICAL                              :: Use2DTimeFunc
  REAL                                 :: StagnationPoint
  REAL                                 :: TimeScaleFac2DFront
  REAL                                 :: TimeScaleFac2DBack
  REAL                                 :: DistributionMaxTimeFactor
  REAL                                 :: DistributionMinTimeFactor
  INTEGER                              :: DistributionMinPartNum
  LOGICAL                              :: AdaptDistribution
  REAL                                 :: TargetMCSoverMFP
  REAL                                 :: TargetMaxCollProb
END TYPE
TYPE(tVariableTimeStep)                :: VarTimeStep


!===================================================================================================================================
END MODULE MOD_Particle_Vars
