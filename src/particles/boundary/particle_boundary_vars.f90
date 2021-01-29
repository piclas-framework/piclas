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
MODULE MOD_Particle_Boundary_Vars
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

LOGICAL                                 :: SurfOnNode
INTEGER                                 :: SurfSampSize                  !> Energy + Force + nSpecies
REAL,ALLOCPOINT,DIMENSION(:,:,:)        :: SurfSideArea                  !> Area of supersampled surface side
! ====================================================================
! Mesh info
INTEGER                                 :: nSurfTotalSides
INTEGER                                 :: nOutputSides

INTEGER                                 :: nComputeNodeSurfSides         !> Number of surface sampling sides on compute node
INTEGER                                 :: nComputeNodeSurfOutputSides   !> Number of output surface sampling sides on compute node (inner BCs only counted once)
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
INTEGER,ALLOCATABLE                     :: RotPeriodicSide2GlobalSide(:) ! Mapping BC-side with PartBoundCond=6 to Global Side ID
INTEGER,ALLOCATABLE                     :: NumRotPeriodicNeigh(:)        ! Number of adjacent Neigbours sites in rotational periodic BC
INTEGER,ALLOCATABLE                     :: RotPeriodicSideMapping(:,:)   ! Mapping between rotational periodic sides.
INTEGER,ALLOCATABLE                     :: SurfSide2RotPeriodicSide(:)   ! Mapping between surf side and periodic sides.

!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: NSurfSample                   ! polynomial degree of particle BC sampling
LOGICAL                                 :: CalcSurfaceImpact             ! Sample average impact energy of particles for each species
!                                                                        ! (trans, rot, vib), impact vector and angle (default=FALSE)
REAL,ALLOCATABLE                        :: XiEQ_SurfSample(:)            ! position of XiEQ_SurfSample
REAL                                    :: dXiEQ_SurfSample              ! deltaXi in [-1,1]
INTEGER                                 :: OffSetSurfSide                ! offset of local surf side
INTEGER                                 :: OffSetInnerSurfSide           ! offset of local inner surf side
INTEGER                                 :: nSurfBC                       ! number of surface side BCs
CHARACTER(LEN=255),ALLOCATABLE          :: SurfBCName(:)                 ! names of belonging surface BC
#if USE_MPI
INTEGER,ALLOCATABLE                     :: OffSetSurfSideMPI(:)          ! integer offset for particle boundary sampling
INTEGER,ALLOCATABLE                     :: OffSetInnerSurfSideMPI(:)     ! integer offset for particle boundary sampling (innerBC)
INTEGER                                 :: nComputeNodeInnerBCs          ! Number of inner BCs with a larger global side ID on node
#endif /*USE_MPI*/

#if USE_MPI
TYPE tSurfaceSendList
  INTEGER                               :: NativeProcID
  INTEGER,ALLOCATABLE                   :: SendList(:)                   ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: RecvList(:)                   ! list containing surfsideid of sides to recv from proc

  INTEGER,ALLOCATABLE                   :: SurfDistSendList(:)           ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: SurfDistRecvList(:)           ! list containing surfsideid of sides to recv from proc
  INTEGER,ALLOCATABLE                   :: H2OSendList(:)                ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: H2ORecvList(:)                ! list containing surfsideid of sides to recv from proc
  INTEGER,ALLOCATABLE                   :: O2HSendList(:)                ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: O2HRecvList(:)                ! list containing surfsideid of sides to recv from proc

END TYPE
#endif /*USE_MPI*/

TYPE tSurfaceCOMM
  LOGICAL                               :: MPIRoot                       ! if root of mpi communicator
  INTEGER                               :: MyRank                        ! local rank in new group
  INTEGER                               :: nProcs                        ! number of processes
  LOGICAL                               :: MPIOutputRoot                 ! if root of mpi communicator
  INTEGER                               :: MyOutputRank                  ! local rank in new group
  INTEGER                               :: nOutputProcs                  ! number of output processes
#if USE_MPI
  LOGICAL                               :: InnerBCs                      ! are there InnerSides with reflective properties
  INTEGER                               :: COMM=MPI_COMM_NULL            ! communicator
  INTEGER                               :: nMPINeighbors                 ! number of processes to communicate with
  TYPE(tSurfaceSendList),ALLOCATABLE    :: MPINeighbor(:)                ! list containing all mpi neighbors
  INTEGER                               :: OutputCOMM=MPI_COMM_NULL      ! communicator for output
#endif /*USE_MPI*/
END TYPE
TYPE (tSurfaceCOMM)                     :: SurfCOMM

TYPE tSurfaceMesh
  INTEGER                               :: SampSize                      ! integer of sampsize
  INTEGER                               :: ReactiveSampSize              ! additional sample size on the surface due to use of
                                                                         ! reactive surface modelling (reactions, liquid, etc.)
  LOGICAL                               :: SurfOnProc                    ! flag if reflective boundary condition is on proc
  INTEGER                               :: nSides                        ! Number of Sides on Surface (reflective)
  INTEGER                               :: nBCSides                      ! Number of OuterSides with Surface (reflective) properties
  INTEGER                               :: nInnerSides                   ! Number of InnerSides with Surface (reflective) properties
  INTEGER                               :: nOutputSides                  ! Number of surfaces that are assigned to an MPI rank for
                                                                         ! surface sampling (MacroSurfaceVal and MacroSurfaceSpecVal)
                                                                         ! and output to .h5 (SurfData) purposes:
                                                                         ! nOutputSides = bcsides + maser_innersides
  !INTEGER                               :: nTotalSides                   ! Number of Sides on Surface incl. HALO sides
  INTEGER                               :: nGlobalSides                  ! Global number of Sides on Surfaces (reflective)
  INTEGER,ALLOCATABLE                   :: SideIDToSurfID(:)             ! Mapping of side ID to surface side ID (reflective)
  REAL, ALLOCATABLE                     :: SurfaceArea(:,:,:)            ! Area of Surface
  INTEGER,ALLOCATABLE                   :: SurfIDToSideID(:)             ! Mapping of surface side ID (reflective) to side ID
  INTEGER,ALLOCATABLE                   :: innerBCSideToHaloMap(:)       ! map of inner BC ID on slave side to corresp. HaloSide
END TYPE

TYPE (tSurfaceMesh)                     :: SurfMesh

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

TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleAnodeBC           = 4      ! = 4 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleCathodeBC         = 5      ! = 5 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: RotPeriodicBC           = 6      ! = 6 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryAxis            = 11     ! = 10 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)          ! Link part 1 for mapping PICLas BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)          ! Link part 2 for mapping PICLas BCs to Particle BC
!  INTEGER              , ALLOCATABLE     :: Map(:)                      ! Map from PICLas BCindex to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)              ! Map from PICLas BCindex to Particle BC (NOT TO TYPE!)
  !!INTEGER              , ALLOCATABLE     :: SideBCType(:)            ! list with boundary condition for each side
  REAL    , ALLOCATABLE                  :: MomentumACC(:)
  REAL    , ALLOCATABLE                  :: WallTemp(:), WallTemp2(:), WallTempDelta(:)
  REAL    , ALLOCATABLE                  :: TempGradStart(:,:), TempGradEnd(:,:), TempGradVec(:,:)
  REAL    , ALLOCATABLE                  :: TransACC(:)
  REAL    , ALLOCATABLE                  :: VibACC(:)
  REAL    , ALLOCATABLE                  :: RotACC(:)
  REAL    , ALLOCATABLE                  :: ElecACC(:)
  REAL    , ALLOCATABLE                  :: WallVelo(:,:)
  REAL    , ALLOCATABLE                  :: Voltage(:)
  LOGICAL , ALLOCATABLE                  :: RotVelo(:)                    ! Flag for rotating walls
  REAL    , ALLOCATABLE                  :: RotFreq(:)                    ! Rotation frequency of the wall
  REAL    , ALLOCATABLE                  :: RotAxi(:,:)                   ! Direction of rotation axis
  REAL    , ALLOCATABLE                  :: RotOrg(:,:)                   ! Origin of rotation axis
  INTEGER , ALLOCATABLE                  :: RotPeriodicDir(:)             ! Direction of rotation
  INTEGER , ALLOCATABLE                  :: NbrOfSpeciesSwaps(:)          ! Number of Species to be changed at wall
  REAL    , ALLOCATABLE                  :: ProbOfSpeciesSwaps(:)         ! Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE                  :: SpeciesSwaps(:,:,:)           ! Species to be changed at wall (in, out), out=0: delete
  INTEGER , ALLOCATABLE                  :: SurfaceModel(:)               ! Model used for surface interaction
                                                                            ! 0 perfect/diffusive reflection
                                                                            ! 5 SEE (secondary e- emission) by Levko2015
                                                                            ! 6 SEE (secondary e- emission) by Pagonakis2016
                                                                            !   (originally from Harrower1956)
  LOGICAL , ALLOCATABLE                  :: Reactive(:)                   ! flag defining if surface is treated reactively
  LOGICAL , ALLOCATABLE                  :: Resample(:)                   ! Resample Equilibrium Distribution with reflection
  LOGICAL , ALLOCATABLE                  :: Dielectric(:)                 ! Define if particle boundary [$] is a dielectric
!                                                                         ! interface, i.e. an interface between a dielectric and
!                                                                         ! a non-dielectric or a between to different dielectrics
!                                                                         ! [.TRUE.] or not [.FALSE.] (requires reflective BC)
!                                                                         ! (Default=FALSE.)
  LOGICAL , ALLOCATABLE                  :: BoundaryParticleOutputHDF5(:) ! Save particle position, velocity and species to
!                                                                         ! PartDataBoundary container for writing to .h5 later
END TYPE

INTEGER                                  :: nPartBound                       ! number of particle boundaries
TYPE(tPartBoundary)                      :: PartBound                         ! Boundary Data for Particles

INTEGER                                  :: nAuxBCs                     ! number of aux. BCs that are checked during tracing
LOGICAL                                  :: UseAuxBCs                     ! number of aux. BCs that are checked during tracing
CHARACTER(LEN=200), ALLOCATABLE          :: AuxBCType(:)                ! type of BC (plane, ...)
INTEGER           , ALLOCATABLE          :: AuxBCMap(:)                 ! index of AuxBC in respective Type

TYPE tAuxBC_plane
  REAL                                   :: r_vec(3)
  REAL                                   :: n_vec(3)
  REAL                                   :: radius
END TYPE tAuxBC_plane
TYPE(tAuxBC_plane), ALLOCATABLE          :: AuxBC_plane(:)

TYPE tAuxBC_cylinder
  REAL                                   :: r_vec(3)
  REAL                                   :: axis(3)
  REAL                                   :: radius
  REAL                                   :: lmin
  REAL                                   :: lmax
  LOGICAL                                :: inwards
END TYPE tAuxBC_cylinder
TYPE(tAuxBC_cylinder), ALLOCATABLE       :: AuxBC_cylinder(:)

TYPE tAuxBC_cone
  REAL                                   :: r_vec(3)
  REAL                                   :: axis(3)
  REAL                                   :: halfangle
  REAL                                   :: lmin
  REAL                                   :: lmax
  REAL                                   :: geomatrix(3,3)
  !REAL                                   :: geomatrix2(3,3)
  REAL                                   :: rotmatrix(3,3)
  LOGICAL                                :: inwards
END TYPE tAuxBC_cone
TYPE(tAuxBC_cone), ALLOCATABLE       :: AuxBC_cone(:)

TYPE tAuxBC_parabol
  REAL                                   :: r_vec(3)
  REAL                                   :: axis(3)
  REAL                                   :: zfac
  REAL                                   :: lmin
  REAL                                   :: lmax
  REAL                                   :: geomatrix4(4,4)
  REAL                                   :: rotmatrix(3,3)
  LOGICAL                                :: inwards
END TYPE tAuxBC_parabol
TYPE(tAuxBC_parabol), ALLOCATABLE       :: AuxBC_parabol(:)

TYPE tPartAuxBC
  INTEGER               :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER               :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER , ALLOCATABLE :: TargetBoundCond(:)
  REAL    , ALLOCATABLE :: MomentumACC(:)
  REAL    , ALLOCATABLE :: WallTemp(:)
  REAL    , ALLOCATABLE :: TransACC(:)
  REAL    , ALLOCATABLE :: VibACC(:)
  REAL    , ALLOCATABLE :: RotACC(:)
  REAL    , ALLOCATABLE :: ElecACC(:)
  REAL    , ALLOCATABLE :: WallVelo(:,:)
  INTEGER , ALLOCATABLE :: NbrOfSpeciesSwaps(:)  !Number of Species to be changed at wall
  REAL    , ALLOCATABLE :: ProbOfSpeciesSwaps(:) !Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE :: SpeciesSwaps(:,:,:)   !Species to be changed at wall (in, out), out=0: delete
  LOGICAL , ALLOCATABLE :: Resample(:)           !Resample Equilibirum Distribution with reflection
END TYPE
TYPE(tPartAuxBC)        :: PartAuxBC             ! auxBC Data for Particles

! Boundary particle output
LOGICAL              :: DoBoundaryParticleOutputHDF5   ! Flag set automatically if particles crossing specific
!                                                  ! boundaries are to be saved to .h5 (position of intersection,
!                                                  ! velocity, species, internal energies)
REAL, ALLOCATABLE    :: PartStateBoundary(:,:)     ! (1:10,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,Ekin,MPF,time,impact angle
!                                                  !                 2nd index: 1 to number of boundary-crossed particles
INTEGER              :: PartStateBoundaryVecLength ! Number of boundary-crossed particles
!===================================================================================================================================

END MODULE MOD_Particle_Boundary_Vars
