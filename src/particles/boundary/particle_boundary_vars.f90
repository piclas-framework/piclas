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
MODULE MOD_Particle_Boundary_Vars
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Globals
USE mpi
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

LOGICAL                                 :: SurfTotalSideOnNode
INTEGER                                 :: SurfSampSize                  !> Energy + Force + nSpecies
INTEGER                                 :: SurfOutputSize                !> Energy + Force + nSpecies
INTEGER                                 :: SurfSpecOutputSize            !> Energy + Force + nSpecies
REAL,ALLOCPOINT,DIMENSION(:,:,:)        :: SurfSideArea                  !> Area of supersampled surface side
REAL,ALLOCPOINT,DIMENSION(:,:,:)        :: BoundaryWallTemp              !> Wall Temperature for Adaptive Case
! ====================================================================
! Mesh info
INTEGER                                 :: nGlobalSurfSides
INTEGER                                 :: nGlobalOutputSides

INTEGER                                 :: nComputeNodeSurfSides         !> Number of surface sampling sides on compute node
INTEGER                                 :: nComputeNodeSurfOutputSides   !> Number of output surface sampling sides on compute node (inner BCs only counted once and rotationally periodic BCs excluded)
INTEGER                                 :: nComputeNodeSurfTotalSides    !> Number of surface sampling sides on compute node (including halo region)
INTEGER                                 :: offsetComputeNodeSurfSide     !> elem offset of compute-node root
INTEGER                                 :: offsetComputeNodeSurfOutputSide     !> elem offset of compute-node root

! ====================================================================
! Impact statistics
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: SampWallState
REAL,ALLOCATABLE,DIMENSION(:)           :: SampWallPumpCapacity
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)   :: SampWallImpactEnergy
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)   :: SampWallImpactVector
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: SampWallImpactAngle
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: SampWallImpactNumber

! SampWallState indices for optional variables (defined in InitParticleBoundarySampling)
INTEGER                                 :: SWIVarTimeStep
INTEGER                                 :: SWIStickingCoefficient

! Output container
REAL,ALLOCATABLE                  :: MacroSurfaceVal(:,:,:,:)           !> variables,p,q,sides
REAL,ALLOCATABLE                  :: MacroSurfaceSpecVal(:,:,:,:,:)     !> Macrovalues for Species specific surface output
                                                                        !> (4,p,q,nSurfSides,nSpecies)
                                                                        !> 1: Surface Collision Counter
                                                                        !> 2: Accommodation
                                                                        !> 3: Coverage
                                                                        !> 4 (or 2): Impact energy trans
                                                                        !> 5 (or 3): Impact energy rot
                                                                        !> 6 (or 4): Impact energy vib

! ====================================================================
! MPI3 shared variables
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: GlobalSide2SurfSide           ! Mapping Global Side ID to Surf Side ID
                                                                         !> 1 - Surf SideID
                                                                         !> 2 - Surf Side proc global rank
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: SurfSide2GlobalSide           ! Inverse mapping
                                                                         !> 1 - Surf SideID
                                                                         !> 2 - Surf Side proc global rank
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: GlobalSide2SurfSide_Shared
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: SurfSide2GlobalSide_Shared

#if USE_MPI
REAL,POINTER,DIMENSION(:,:,:)           :: BoundaryWallTemp_Shared           !> Wall Temperature for Adaptive Case
INTEGER                                 :: BoundaryWallTemp_Shared_Win

REAL,POINTER,DIMENSION(:,:,:)           :: SurfSideArea_Shared           !> Area of supersampled surface side
INTEGER                                 :: SurfSideArea_Shared_Win

INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: GlobalSide2SurfHaloSide       ! Mapping Global Side ID to Surf Halo Side ID (exists only on leader procs)
                                                                         !> 1st dim: leader rank
                                                                         !> 2nd dim: Surf SideID
INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: SurfHaloSide2GlobalSide       ! Inverse mapping  (exists only on leader procs)
                                                                         !> 1st dim: leader rank
                                                                         !> 2nd dim: Surf SideID

INTEGER                                 :: GlobalSide2SurfSide_Shared_Win
INTEGER                                 :: SurfSide2GlobalSide_Shared_Win

TYPE tSurfaceMapping
  INTEGER,ALLOCATABLE                   :: RecvSurfGlobalID(:)
  INTEGER,ALLOCATABLE                   :: SendSurfGlobalID(:)
  INTEGER                               :: nSendSurfSides
  INTEGER                               :: nRecvSurfSides
  INTEGER,ALLOCATABLE                   :: RecvPorousGlobalID(:)
  INTEGER,ALLOCATABLE                   :: SendPorousGlobalID(:)
  INTEGER                               :: nSendPorousSides
  INTEGER                               :: nRecvPorousSides
END TYPE
TYPE (tSurfaceMapping),ALLOCATABLE      :: SurfMapping(:)

! ====================================================================
! Impact statistics
REAL,POINTER,DIMENSION(:,:,:,:)         :: SampWallState_Shared
REAL,POINTER,DIMENSION(:)               :: SampWallPumpCapacity_Shared
REAL,POINTER,DIMENSION(:,:,:,:,:)       :: SampWallImpactEnergy_Shared
REAL,POINTER,DIMENSION(:,:,:,:,:)       :: SampWallImpactVector_Shared
REAL,POINTER,DIMENSION(:,:,:,:)         :: SampWallImpactAngle_Shared
REAL,POINTER,DIMENSION(:,:,:,:)         :: SampWallImpactNumber_Shared

INTEGER                                 :: SampWallState_Shared_Win
INTEGER                                 :: SampWallPumpCapacity_Shared_Win
INTEGER                                 :: SampWallImpactEnergy_Shared_Win
INTEGER                                 :: SampWallImpactVector_Shared_Win
INTEGER                                 :: SampWallImpactAngle_Shared_Win
INTEGER                                 :: SampWallImpactNumber_Shared_Win
#endif /* USE_MPI */

! ====================================================================
! Rotational periodic sides
INTEGER                           :: nRotPeriodicSides         ! Number of rotational periodic sides on a compute node
INTEGER                           :: MaxNumRotPeriodicNeigh    ! Maximum number of rotationally periodic neighbours
INTEGER,ALLOCPOINT,DIMENSION(:)   :: NumRotPeriodicNeigh       ! Number of adjacent Neigbours sites in rotational periodic BC
INTEGER,ALLOCPOINT,DIMENSION(:,:) :: RotPeriodicSideMapping    ! Mapping between rotational periodic sides.
INTEGER,ALLOCPOINT,DIMENSION(:)   :: SurfSide2RotPeriodicSide  ! Mapping between surf side and periodic sides.
! ====================================================================
! Intermediate plane for multi rotational periodic sides
INTEGER,ALLOCPOINT,DIMENSION(:,:) :: InterPlaneSideMapping    ! Mapping between inter plane BC_ID and SideID.
! ====================================================================
#if USE_MPI
INTEGER,POINTER,DIMENSION(:)    :: SurfSide2RotPeriodicSide_Shared
INTEGER                         :: SurfSide2RotPeriodicSide_Shared_Win
INTEGER,POINTER,DIMENSION(:)    :: NumRotPeriodicNeigh_Shared
INTEGER                         :: NumRotPeriodicNeigh_Shared_Win
INTEGER,POINTER,DIMENSION(:)    :: Rot2Glob_temp_Shared
INTEGER                         :: Rot2Glob_temp_Shared_Win
INTEGER,POINTER,DIMENSION(:,:)  :: RotPeriodicSideMapping_temp_Shared
INTEGER                         :: RotPeriodicSideMapping_temp_Shared_Win
INTEGER,POINTER,DIMENSION(:,:)  :: RotPeriodicSideMapping_Shared
INTEGER                         :: RotPeriodicSideMapping_Shared_Win
REAL,POINTER,DIMENSION(:,:)     :: BoundingBox_Shared
INTEGER                         :: BoundingBox_Shared_Win
#endif /*USE_MPI*/

!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: NSurfSample                   ! polynomial degree of particle BC sampling
LOGICAL                                 :: CalcSurfaceImpact             ! Sample average impact energy of particles for each species
!                                                                        ! (trans, rot, vib), impact vector and angle (default=FALSE)
REAL,ALLOCATABLE                        :: XiEQ_SurfSample(:)            ! position of XiEQ_SurfSample
REAL                                    :: dXiEQ_SurfSample              ! deltaXi in [-1,1]
INTEGER                                 :: OffSetSurfSide                ! offset of local surf side
INTEGER                                 :: nSurfBC                       ! number of surface side BCs
CHARACTER(LEN=255),ALLOCATABLE          :: SurfBCName(:)                 ! names of belonging surface BC
#if USE_MPI
INTEGER                                 :: nComputeNodeInnerBCs(2)       ! Number of inner BCs with a larger global side ID on node
#endif /*USE_MPI*/

#if USE_MPI
TYPE tMPIGROUP
  INTEGER                     :: UNICATOR=MPI_COMM_NULL !< MPI communicator for surface sides (including sides inside the halo region)
  INTEGER                     :: nProcs                 !< number of MPI processes for particles
  INTEGER                     :: MyRank                 !< MyRank within communicator
END TYPE
TYPE (tMPIGROUP)              :: SurfCOMM
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! Porous BC
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: nPorousSides                       ! Number of porous sides per compute node
INTEGER,ALLOCPOINT                      :: MapSurfSideToPorousSide_Shared(:)  ! Mapping of surface side to porous side
INTEGER,ALLOCPOINT                      :: PorousBCInfo_Shared(:,:)           ! Info and mappings for porous BCs [1:3,1:nPorousSides]
                                                                              ! 1: Porous BC ID
                                                                              ! 2: SurfSide ID
                                                                              ! 3: SideType (0: inside, 1: partially)
REAL,ALLOCPOINT                         :: PorousBCProperties_Shared(:,:)     ! Properties of porous sides [1:2,1:nPorousSides]
                                                                              ! 1: Removal probability
                                                                              ! 2: Pumping speed per side
REAL,ALLOCPOINT                         :: PorousBCSampWall_Shared(:,:)       ! Sampling of impinging and deleted particles on each
                                                                              ! porous sides [1:2,1:nPorousSides]
                                                                              ! REAL variable since the particle weight is used
                                                                              ! 1: Impinging particles
                                                                              ! 2: Deleted particles
#if USE_MPI
INTEGER                                 :: MapSurfSideToPorousSide_Shared_Win
INTEGER                                 :: PorousBCInfo_Shared_Win
INTEGER                                 :: PorousBCProperties_Shared_Win
INTEGER                                 :: PorousBCSampWall_Shared_Win
#endif

REAL,ALLOCATABLE                        :: PorousBCSampWall(:,:)  ! Processor-local sampling of impinging and deleted particles
                                                                  ! REAL variable since the particle weight is used
                                                                  ! 1: Impinging particles
                                                                  ! 2: Deleted particles
!-----------------------------------------------------------------------------------------------------------------------------------
! Particle Boundary
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: RotPeriodicBC           = 6      ! = 6 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: RotPeriodicInterPlaneBC = 7      ! = 7 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryAxis            = 11     ! = 10 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)           ! Link part 1 for mapping PICLas BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)           ! Link part 2 for mapping PICLas BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)               ! Map from PICLas BCindex to Particle BC (NOT TO TYPE!)
  INTEGER              , ALLOCATABLE     :: MapToFieldBC(:)              ! Map from Particle BC (NOT TO TYPE!) to PICLas BCindex
  ! Constant wall temperature and accommodation coefficients
  REAL    , ALLOCATABLE                  :: WallTemp(:)
  REAL    , ALLOCATABLE                  :: MomentumACC(:)
  LOGICAL , ALLOCATABLE                  :: OnlySpecular(:)
  LOGICAL , ALLOCATABLE                  :: OnlyDiffuse(:)
  REAL    , ALLOCATABLE                  :: TransACC(:)
  REAL    , ALLOCATABLE                  :: VibACC(:)
  REAL    , ALLOCATABLE                  :: RotACC(:)
  REAL    , ALLOCATABLE                  :: ElecACC(:)
  ! Temperature gradient across reflective BC
  REAL    , ALLOCATABLE                  :: WallTemp2(:), WallTempDelta(:)
  REAL    , ALLOCATABLE                  :: TempGradStart(:,:), TempGradEnd(:,:), TempGradVec(:,:)
  INTEGER , ALLOCATABLE                  :: TempGradDir(:)
  ! Linear and rotational wall velocity
  REAL    , ALLOCATABLE                  :: WallVelo(:,:)
  REAL    , ALLOCATABLE                  :: PhotonEnACC(:)
  REAL    , ALLOCATABLE                  :: PhotonSEEYield(:)
  REAL    , ALLOCATABLE                  :: PhotonSEEWorkFunction(:)
  REAL    , ALLOCATABLE                  :: PhotonSEEMacroParticleFactor(:)
  INTEGER , ALLOCATABLE                  :: PhotonSEEElectronSpecies(:)
  LOGICAL , ALLOCATABLE                  :: PhotonSpecularReflection(:)
  LOGICAL , ALLOCATABLE                  :: RotVelo(:)                    ! Flag for rotating walls
  REAL    , ALLOCATABLE                  :: RotOmega(:,:)                 ! Angular velocity
  ! Species swap BCs
  INTEGER , ALLOCATABLE                  :: NbrOfSpeciesSwaps(:)          ! Number of Species to be changed at wall
  REAL    , ALLOCATABLE                  :: ProbOfSpeciesSwaps(:)         ! Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE                  :: SpeciesSwaps(:,:,:)           ! Species to be changed at wall (in, out), out=0: delete
  ! Surface models
  INTEGER , ALLOCATABLE                  :: SurfaceModel(:)               ! Model used for surface interaction (e.g. SEE models)
  LOGICAL , ALLOCATABLE                  :: Reactive(:)                   ! flag defining if surface is treated reactively
  LOGICAL , ALLOCATABLE                  :: Resample(:)                   ! Resample Equilibrium Distribution with reflection
  ! Radiative-equilibrium BC
  LOGICAL                                :: AdaptWallTemp
  LOGICAL , ALLOCATABLE                  :: UseAdaptedWallTemp(:)
  LOGICAL                                :: OutputWallTemp                ! Flag to include the wall temperature in the SurfState
                                                                          ! output, set during InitializeVariablesPartBoundary
                                                                          ! Required for the initialization of the array for the
                                                                          ! adaptive wall temperature as well
  REAL    , ALLOCATABLE                  :: RadiativeEmissivity(:)
  ! Dielectric BC
  LOGICAL , ALLOCATABLE                  :: Dielectric(:)                 ! Define if particle boundary [$] is a dielectric
                                                                          ! interface, i.e. an interface between a dielectric and
                                                                          ! a non-dielectric or a between to different dielectrics
                                                                          ! [.TRUE.] or not [.FALSE.] (requires reflective BC)
                                                                          ! (Default=FALSE.)
  ! Multi rotational periodic and interplane BCs
  LOGICAL                                :: UseRotPeriodicBC            ! Flag for rotational periodicity
  LOGICAL                                :: OutputBCDataForTesting      ! Flag to output boundary parameter which were determined
                                                                        ! automatically
  INTEGER                                :: RotPeriodicAxis             ! Axis of rotational periodicity
  REAL                                   :: RotPeriodicTol              ! Tolerance for rotationally periodic BC, angle is multiplied
                                                                        ! by 1 - RotPeriodicTol
  REAL    , ALLOCATABLE                  :: RotPeriodicAngle(:)         ! Angle and direction of rotation [1:nPartBound]
  REAL    , ALLOCATABLE                  :: RotPeriodicMin(:)           ! Min rot axi value [1:nPartBound]
  REAL    , ALLOCATABLE                  :: RotPeriodicMax(:)           ! Max rot axi value [1:nPartBound]
  LOGICAL                                :: UseInterPlaneBC             ! Flag for inter planes exist
  INTEGER , ALLOCATABLE                  :: AssociatedPlane(:)          ! Link between both coressponding intermediate planes
  INTEGER , ALLOCATABLE                  :: nSidesOnInterPlane(:)       ! Number of Sides on intermediate plane
  REAL    , ALLOCATABLE                  :: NormalizedRadiusDir(:,:)    ! Normalized vector in radius direction that is used to
                                                                        ! calculate a random position on same radius within the
                                                                        ! rot periodic segment
  REAL    , ALLOCATABLE                  :: RotAxisPosition(:)          ! Position of inter plane at rotation axis
  REAL    , ALLOCATABLE                  :: AngleRatioOfInterPlanes(:)  ! Ratio of rotation angles for the intermediate planes
  ! Boundary particle output
  LOGICAL , ALLOCATABLE                  :: BoundaryParticleOutputHDF5(:) ! Save particle position, velocity and species to
                                                                          ! PartDataBoundary container for writing to .h5 later
END TYPE

INTEGER                                  :: nPartBound                    ! number of particle boundaries
TYPE(tPartBoundary)                      :: PartBound                     ! Boundary Data for Particles

!-----------------------------------------------------------------------------------------------------------------------------------
! Boundary particle output
LOGICAL              :: DoBoundaryParticleOutputHDF5 ! Flag set automatically if particles crossing specific  boundaries are to be saved to .h5 (position of intersection, velocity, species, internal energies)
LOGICAL              :: DoBoundaryParticleOutputRay ! User-defined flag to output surface SEE or volume ionization emission particles to .h5 based on the ray tracing model
REAL, ALLOCATABLE    :: PartStateBoundary(:,:)     ! (1:11,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,Ekin,MPF,time,impact angle, BCindex
!                                                  !                 2nd index: 1 to number of boundary-crossed particles
INTEGER, PARAMETER   :: nVarPartStateBoundary=11
INTEGER              :: PartStateBoundaryVecLength ! Number of boundary-crossed particles
!===================================================================================================================================

END MODULE MOD_Particle_Boundary_Vars
