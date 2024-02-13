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

MODULE MOD_Particle_Vars
!===================================================================================================================================
! Contains the Particles' variables (general for all modules: PIC, DSMC, FP)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Emission_Vars
USE MOD_Particle_SurfaceFlux_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
#ifdef INTKIND8
INTEGER, PARAMETER :: IK = SELECTED_INT_KIND(18)
#else
INTEGER, PARAMETER :: IK = SELECTED_INT_KIND(8)
#endif

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
REAL    , ALLOCATABLE :: PartVeloRotRef(:,:)                                 ! (1:3,1:NParts) Velocity in the rotational reference frame
REAL    , ALLOCATABLE :: Pt(:,:)                                             ! Derivative of PartState (vx,xy,vz) only
                                                                             ! since temporal derivative of position
                                                                             ! is the velocity. Thus we can take
                                                                             ! PartState(4:6,:) as Pt(1:3)
                                                                             ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
INTEGER               :: PartDataSize                                        ! Number of entries in each line of PartData
CHARACTER(LEN=255),ALLOCATABLE :: PartDataVarNames(:)                        ! Corrensponding variable names of PartData for output/read-in
INTEGER,PARAMETER     :: PartIntSize=2                                       ! Number of entries in each line of PartInt
REAL,ALLOCATABLE      :: PartData(:,:)                                       ! PartState ordered along SFC, particle number per
                                                                             ! element given in PartInt
INTEGER,ALLOCATABLE   :: VibQuantData(:,:)                                   ! Vibrational quantization
INTEGER               :: MaxQuantNum
REAL,ALLOCATABLE      :: ElecDistriData(:,:)                                 ! Electronic excitation distribution
INTEGER               :: MaxElecQuant
REAL,ALLOCATABLE      :: AD_Data(:,:)                                        ! Ambipolar diffusion

INTEGER(KIND=IK),ALLOCATABLE :: PartInt(:,:)                                 ! Particle number per element
INTEGER(KIND=IK)             :: locnPart,offsetnPart                         ! Number and offset of particles on processors
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
INTEGER , ALLOCATABLE :: InterPlanePartIndx(:)                               ! Index list of Particles that are created
                                                                             ! at rot periodic inter plane
INTEGER               :: InterPlanePartNumber                                ! Number of Particles that are created
                                                                             ! at rot periodic inter plane
REAL    , ALLOCATABLE :: PartTimeStep(:)                                     ! (1:NParts) Variable time step
INTEGER               :: PartLorentzType
CHARACTER(LEN=256)    :: ParticlePushMethod                                  ! Type of PP-Method
INTEGER               :: nrSeeds                                             ! Number of Seeds for Random Number Generator
INTEGER , ALLOCATABLE :: seeds(:)                        !        =>NULL()   ! Seeds for Random Number Generator

TYPE tSpecies                                                                ! Particle Data for each Species
  !General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  REAL                                   :: ChargeIC                         ! Particle Charge (without MPF)
  REAL                                   :: MassIC                           ! Particle Mass (without MPF)
  REAL                                   :: MacroParticleFactor              ! Number of Microparticle per Macroparticle
  REAL                                   :: TimeStepFactor                   ! Species-specific time step factor
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  TYPE(typeSurfaceflux),ALLOCATABLE      :: Surfaceflux(:)                   ! Particle Data for each SurfaceFlux emission
  INTEGER                                :: nSurfacefluxBCs                  ! Number of SF emissions
#if IMPA
  LOGICAL                                :: IsImplicit
#endif
END TYPE

INTEGER                                  :: nSpecies                         ! number of species
CHARACTER(LEN=256)                       :: SpeciesDatabase                  ! Name of the species database
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
  INTEGER                                :: maxAllowedParticleNumber          ! Maximum allowed number of PDM%maxParticleNumber
  LOGICAL                                :: RearrangePartIDs                  ! Rearrange PartIDs during shrinking maxPartNum
#if USE_MPI
  LOGICAL                                :: UNFPafterMPIPartSend              ! UpdateNextFreePosition after MPI Part Send
#endif
  INTEGER                                :: ParticleVecLength                 ! Vector Length for Particle Push Calculation
  INTEGER                                :: ParticleVecLengthOld              ! Vector Length for Particle Push Calculation
  REAL                                   :: MaxPartNumIncrease                ! How much shall the PDM%MaxParticleNumber be increased if it is full
  INTEGER ,ALLOCATABLE                   :: nextFreePosition(:)  !  =>NULL()  ! next_free_Position(1:maxParticleNumber)
                                                                              ! List of free Positon
  LOGICAL ,ALLOCATABLE                   :: ParticleInside(:)                 ! Particle_inside (1:maxParticleNumber)
  LOGICAL ,ALLOCATABLE                   :: InRotRefFrame(:)                  ! Check for RotRefFrame (1:maxParticleNumber)
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
! Sampling of electronic excitation rates
LOGICAL                                  :: SampleElecExcitation              ! Enable sampling the electronic excitation rate per level
INTEGER                                  :: ExcitationLevelCounter            ! Counter of electronic levels to be sampled (for all species)
REAL, ALLOCATABLE                        :: ExcitationSampleData(:,:)         ! Sampled rates [1:ExcitationLevelCounter,1:nElems]
INTEGER, ALLOCATABLE                     :: ExcitationLevelMapping(:,:)       ! Mapping of collision case and level to the total electronic
                                                                              ! number of levels [1:CollInf%NumCase,1:MAXVAL(SpecXSec(:)%NumElecLevel)]
! Variable particle weighting (vMPF)
LOGICAL                                  :: usevMPF                           ! use the vMPF per particle
INTEGER, ALLOCATABLE                     :: vMPFMergeThreshold(:)             ! Max particle number per cell and (iSpec)
INTEGER, ALLOCATABLE                     :: vMPFSplitThreshold(:)             ! Min particle number per cell and (iSpec)
REAL                                     :: vMPFSplitLimit                    ! Do not split particles below this MPF threshold
LOGICAL                                  :: UseSplitAndMerge                  ! Flag for particle merge
REAL, ALLOCATABLE                        :: CellEelec_vMPF(:,:)
REAL, ALLOCATABLE                        :: CellEvib_vMPF(:,:)

! Surface flux flags
LOGICAL                                  :: DoSurfaceFlux                     ! Flag for emitting by SurfaceFluxBCs
LOGICAL                                  :: DoPoissonRounding                 ! Perform Poisson sampling instead of random rounding
LOGICAL                                  :: DoTimeDepInflow                   ! Insertion and SurfaceFlux w simple random rounding
LOGICAL                                  :: DoZigguratSampling                ! Sample normal randoms with Ziggurat method

! Variable time step
LOGICAL                                :: UseVarTimeStep
TYPE tVariableTimeStep
  LOGICAL                              :: UseLinearScaling
  LOGICAL                              :: UseDistribution
  LOGICAL                              :: UseSpeciesSpecific
  LOGICAL                              :: OnlyDecreaseDt
  LOGICAL                              :: DisableForMCC
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

! Virtual cell merge
LOGICAL                                :: DoVirtualCellMerge
INTEGER                                :: MinPartNumCellMerge
INTEGER                                :: VirtualCellMergeSpread
INTEGER                                :: MaxNumOfMergedCells
TYPE tVirtualCellMerge
  INTEGER, ALLOCATABLE                 :: MergedCellID(:)
  INTEGER                              :: NumOfMergedCells
  INTEGER                              :: MasterCell
  LOGICAL                              :: isMerged
  REAL                                 :: MergedVolume
END TYPE
TYPE (tVirtualCellMerge),ALLOCATABLE   :: VirtMergedCells(:)

! Rotational frame of reference
LOGICAL               :: UseRotRefFrame           ! flag for rotational frame of reference
INTEGER               :: RotRefFrameAxis          ! axis of rotational frame of reference (x=1, y=2, z=3)
REAL                  :: RotRefFrameFreq          ! frequency of rotational frame of reference
REAL                  :: RotRefFrameOmega(3)      ! angular velocity of rotational frame of reference
INTEGER               :: nRefFrameRegions         ! number of rotational frame of reference regions
REAL, ALLOCATABLE     :: RotRefFramRegion(:,:)    ! MIN/MAX defintion for multiple rotational frame of reference region
                                                  ! (i,RegionNumber), MIN:i=1, MAX:i=2
!===================================================================================================================================
END MODULE MOD_Particle_Vars
