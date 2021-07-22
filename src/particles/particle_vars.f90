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
TYPE tSymmetry
  INTEGER             :: Order                                               ! 1-3 D
  LOGICAL             :: Axisymmetric
END TYPE tSymmetry

TYPE(tSymmetry)       :: Symmetry

LOGICAL               :: DoFieldIonization                                   ! Do Field Ionization by quantum tunneling
INTEGER               :: FieldIonizationModel                                !'Field Ionization models. Implemented models are:
!                                                                            ! * Ammosov-Delone-Krainov (ADK) model
!                                                                            ! * Ammosov-Delone-Krainov (ADK) model Yu 2018
LOGICAL,ALLOCATABLE   :: SpecReset(:)                                        ! Flag for resetting species distribution with init
                                                                             ! during restart
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
REAL    , ALLOCATABLE :: PartState(:,:)                                      ! 1st index: x,y,z,vx,vy,vz
!                                                                            ! 2nd index: 1:NParts
REAL    , ALLOCATABLE :: PartPosRef(:,:)                                     ! (1:3,1:NParts) particles pos mapped to -1|1 space
INTEGER , ALLOCATABLE :: PartPosGauss(:,:)                                   ! (1:NParts,1:3) Gauss point localization of particles
REAL    , ALLOCATABLE :: Pt(:,:)                                             ! Derivative of PartState (vx,xy,vz) only
                                                                             ! since temporal derivative of position
                                                                             ! is the velocity. Thus we can take
                                                                             ! PartState(4:6,:) as Pt(1:3)
                                                                             ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
LOGICAL               :: DoForceFreeSurfaceFlux                              ! switch if the stage reconstruction uses a force
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
LOGICAL               :: velocityOutputAtTime
REAL    , ALLOCATABLE :: velocityAtTime(:,:)
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
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
REAL    , ALLOCATABLE :: LastPartPos(:,:)                                    ! 1st index: x,y,z
!                                                                            ! 2nd index: 1:NParts with 2nd index
INTEGER , ALLOCATABLE :: PartSpecies(:)                                      ! (1:NParts)
REAL    , ALLOCATABLE :: PartMPF(:)                                          ! (1:NParts) MacroParticleFactor by variable MPF
INTEGER               :: PartLorentzType
CHARACTER(LEN=256)    :: ParticlePushMethod                                  ! Type of PP-Method
INTEGER               :: nrSeeds                                             ! Number of Seeds for Random Number Generator
INTEGER , ALLOCATABLE :: seeds(:)                        !        =>NULL()   ! Seeds for Random Number Generator

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
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  REAL                                   :: RadiusIC                         ! Radius for IC circle
  REAL                                   :: Radius2IC                        ! Radius2 for IC cylinder (ring)
  REAL                                   :: RadiusICGyro                     ! Radius for Gyrotron gyro radius
  REAL                                   :: InflowRiseTime                   ! time to ramp the number of inflow particles
                                                                             ! linearly from zero to unity
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
  REAL                                   :: VeloIC                           ! velocity for inital Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: Amplitude                        ! Amplitude for sin-deviation initiation.
  REAL                                   :: WaveNumber                       ! WaveNumber for sin-deviation initiation.
  INTEGER                                :: maxParticleNumberX               ! Maximum Number of all Particles in x direction
  INTEGER                                :: maxParticleNumberY               ! Maximum Number of all Particles in y direction
  INTEGER                                :: maxParticleNumberZ               ! Maximum Number of all Particles in z direction
  REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  INTEGER                                :: ParticleEmissionType             ! Emission Type 0 = only initial,
                                                                             !               1 = emission rate in 1/s,
                                                                             !               2 = emission rate 1/iteration
  REAL                                   :: ParticleNumber                   ! Initial, Emission in [1/s] or [1/Iteration]
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=4)                        :: InsertedParticleMisMatch=0       ! error in number of inserted particles of last step
  INTEGER                                :: NumberOfExcludeRegions           ! Number of different regions to be excluded
  TYPE(tExcludeRegion), ALLOCATABLE      :: ExcludeRegion(:)
#if USE_MPI
  INTEGER                                :: InitComm                          ! number of init-communicator
#endif /*USE_MPI*/
!====================================photo ionization =======================================================
  LOGICAL                            :: FirstQuadrantOnly  ! Only insert particles in the first quadrant that is spanned by the
                                                           ! vectors x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R]
  REAL                               :: PulseDuration      ! Pulse duration tau for a Gaussian-type pulse with
                                                           ! I~exp(-(t/tau)^2) [s]
  REAL                               :: WaistRadius        ! Beam waist radius (in focal spot) w_b for Gaussian-type pulse with
                                                           ! I~exp(-(r/w_b)^2) [m]
  REAL                               :: IntensityAmplitude ! Beam intensity maximum I0 Gaussian-type pulse with
                                                           ! I=I0*exp(-(t/tau)^2)exp(-(r/w_b)^2) [W/m^2]
  REAL                               :: WaveLength         ! Beam wavelength [m]
  REAL                               :: YieldSEE           ! Secondary photoelectron yield [-]
  REAL                               :: RepetitionRate     ! Pulse repetition rate [Hz]
  REAL                               :: Power              ! Average pulse power (energy of a single pulse times repetition rate) [J]
  REAL                               :: Energy             ! Single pulse energy (used when RepetitionRate and Power are not supplied [J]
  REAL                               :: Period             ! Time between the maximum intensity of two pulses [s]
  REAL                               :: tActive            ! Pulse will end at tActive (pulse time) [s]
  REAL                               :: tShift             ! Time shift for pulse corresponding to half of the Pulse width (pulse time) [s]
  INTEGER                            :: NbrOfPulses        ! Number of pulses [-]
  REAL                               :: NINT_Correction    ! nearest integer correction factor due to cut-off when converting
                                                           ! the number of particles calculated as real to integer for the
                                                           ! actual emission
  REAL                               :: WorkFunctionSEE    ! Photoelectron work function [eV]
  !REAL                               :: AngularBetaSEE
  REAL                               :: EffectiveIntensityFactor ! Scaling factor that increases I0 [-]
  INTEGER                            :: sumOfMatchedParticles    ! Sum of matched particles on all procs
  INTEGER                            :: sumOfRequestedParticles  ! Sum of requested particles on all procs
  INTEGER                            :: mySumOfMatchedParticles  ! Sum of matched particles on current proc
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
  LOGICAL                                :: CircularInflow                   ! Circular region, which can be used to define small
                                                                             ! geometry features on large boundaries
  INTEGER                                :: dir(3)                           ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                                   :: origin(2)                        ! origin in orth. coordinates of polar system
  REAL                                   :: rmax                             ! max radius of to-be inserted particles
  REAL                                   :: rmin                             ! min radius of to-be inserted particles
  INTEGER, ALLOCATABLE                   :: SurfFluxSideRejectType(:)        ! Type if parts in side can be rejected (1:SideNumber)
  LOGICAL                                :: Adaptive                         ! Is the surface flux an adaptive boundary?
  INTEGER                                :: AdaptiveType                     ! Chose the adaptive type, description in DefineParams
  REAL                                   :: AdaptiveMassflow                 ! Mass flow [kg/s], which is held constant
  REAL                                   :: AdaptivePressure                 ! Static pressure [Pa], which is held constant
  REAL, ALLOCATABLE                      :: ConstMassflowWeight(:,:,:)       ! Adaptive, Type 4: Weighting factor for SF-sides to
                                                                             ! insert the right amount of particles
  REAL, ALLOCATABLE                      :: CircleAreaPerTriaSide(:,:,:)     ! Adaptive, Type 4: Area within a triangle, determined
                                                                             ! through Monte Carlo integration (initially)
  REAL                                   :: SampledMassflow                  ! Actual mass flow rate through a surface flux boundary
  REAL, ALLOCATABLE                      :: nVFRSub(:,:)                     ! normal volume flow rate through subsubside
END TYPE

TYPE tSpecies                                                                ! Particle Data for each Species
  !General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  REAL                                   :: ChargeIC                         ! Particle Charge (without MPF)
  REAL                                   :: MassIC                           ! Particle Mass (without MPF)
  REAL                                   :: MacroParticleFactor              ! Number of Microparticle per Macroparticle
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  TYPE(typeSurfaceflux),ALLOCATABLE      :: Surfaceflux(:)                   ! Particle Data for each SurfaceFlux emission
  INTEGER                                :: nSurfacefluxBCs                  ! Number of SF emissions
#if IMPA
  LOGICAL                                :: IsImplicit
#endif
END TYPE

LOGICAL                                 :: UseCircularInflow              ! Flag is set if the circular inflow feature is used:
                                                                          ! Particle insertion only in the defined circular area
                                                                          ! on the surface of a surface flux
INTEGER, ALLOCATABLE                    :: CountCircInflowType(:,:,:)     ! Counter whether cells are inside/partially inside or
                                                                          ! outside of circular region (only with CODE_ANALYZE)
INTEGER                                  :: nSpecies                         ! number of species
TYPE(tSpecies), ALLOCATABLE              :: Species(:)  !           => NULL() ! Species Data Vector

LOGICAL                                  :: PartMeshHasPeriodicBCs
#if defined(IMPA) || defined(ROS)
LOGICAL                                  :: PartMeshHasReflectiveBCs
#endif
TYPE tParticleElementMapping
  INTEGER                , ALLOCATABLE   :: GlobalElemID(:)     ! =>NULL() ! Current global element number assigned to each Particle
  INTEGER                , ALLOCATABLE   :: LastGlobalElemID(:) ! =>NULL() ! Global element number of the old particle position

  PROCEDURE(ElemID_INTERFACE),POINTER,NOPASS :: LocalElemID !< pointer defining the mapping : global element ID -> local element ID
                                                            !< the function simply returns  : PEM%GlobalElemID(iPart) - offsetElem

  PROCEDURE(ElemID_INTERFACE),POINTER,NOPASS :: CNElemID    !< pointer defining the mapping : global element ID -> compute-node element ID
                                                            !< the function simply returns  : GlobalElem2CNTotalElem(PEM%GlobalElemID(iPart))
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
  INTEGER                , ALLOCATABLE    :: pEnd(:)           !     =>NULL()  ! End of Linked List for Particles in Element
                                                               !               ! pEnd(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pNext(:)          !     =>NULL()  ! Next Particle in same Element (Linked List)
                                                                               ! pStart(1:PIC%maxParticleNumber)
END TYPE

TYPE(tParticleElementMapping)            :: PEM

ABSTRACT INTERFACE
  PPURE INTEGER FUNCTION ElemID_INTERFACE(iPart)
    INTEGER,INTENT(IN) :: iPart
  END FUNCTION
END INTERFACE

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

INTEGER                                  :: vMPFNewPartNum
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
LOGICAL                                  :: DoSurfaceFlux                     ! Flag for emitting by SurfaceFluxBCs
LOGICAL                                  :: DoPoissonRounding                 ! Perform Poisson sampling instead of random rounding
LOGICAL                                  :: DoTimeDepInflow                   ! Insertion and SurfaceFlux w simple random rounding
LOGICAL                                  :: DoZigguratSampling                ! Sample normal randoms with Ziggurat method

TYPE tVariableTimeStep
  LOGICAL                              :: UseVariableTimeStep
  LOGICAL                              :: UseLinearScaling
  LOGICAL                              :: UseDistribution
  LOGICAL                              :: OnlyDecreaseDt
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
  REAL                                 :: TargetMaxRelaxFactor
END TYPE
TYPE(tVariableTimeStep)                :: VarTimeStep

! 2D Landmark
REAL, ALLOCATABLE :: PartPosLandmark(:,:)        ! Store particle positions during emission for placing
!                                                ! Electrons and ions at the exact same position
INTEGER           :: NbrOfParticleLandmarkMax    ! Array maximum size for storing positions
INTEGER           :: FractNbrOld,chunkSizeOld    ! Auxiliary integers for storing positions

LOGICAL           :: UseNeutralization           ! Flag for counting the charged particles impinging on a surface
CHARACTER(255)    :: NeutralizationSource        ! Name of the boundary for calculating the particle balance
INTEGER           :: NeutralizationBalance       ! Counter for charged particles (processor local): Add +1 for electrons and -1 for ions
INTEGER           :: NeutralizationBalanceGlobal ! Counter for charged particles (global): Add +1 for electrons and -1 for ions

!===================================================================================================================================
END MODULE MOD_Particle_Vars
