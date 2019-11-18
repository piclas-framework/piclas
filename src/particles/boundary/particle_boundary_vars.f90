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
  INTEGER                               :: nTotalSides                   ! Number of Sides on Surface incl. HALO sides
  INTEGER                               :: nGlobalSides                  ! Global number of Sides on Surfaces (reflective)
  INTEGER,ALLOCATABLE                   :: SideIDToSurfID(:)             ! Mapping of side ID to surface side ID (reflective)
  REAL, ALLOCATABLE                     :: SurfaceArea(:,:,:)            ! Area of Surface
  INTEGER,ALLOCATABLE                   :: SurfIDToSideID(:)             ! Mapping of surface side ID (reflective) to side ID
  INTEGER,ALLOCATABLE                   :: innerBCSideToHaloMap(:)       ! map of inner BC ID on slave side to corresp. HaloSide
END TYPE

TYPE (tSurfaceMesh)                     :: SurfMesh

TYPE tSampWall             ! DSMC sample for Wall
  ! easier to communicate
  REAL,ALLOCATABLE                      :: State(:,:,:)                ! 1-3     E_tra (pre, wall, re),
                                                                       ! 4-6     E_rot (pre, wall, re),
                                                                       ! 7-9     E_vib (pre, wall, re)
                                                                       ! 10-12   Forces in x, y, z direction
                                                                       ! 13-12+nSpecies   Wall-Collision counter
                                                                       ! 12+nSpecies+1    ParticleTimeStep
  REAL,ALLOCATABLE                      :: SurfModelState(:,:,:)       ! Sampling of reaction enthalpies and coverage
                                                                       ! first index represents
                                                                       ! 1: Heatflux from recombination reactions of two or
                                                                       !    species on the surface.
                                                                       ! 2: Heatflux from dissociation reactions of two or
                                                                       !    species on the surface.
                                                                       ! 3: Heatflux from recombination reactions of one gas
                                                                       !    species reacting at collision with another species
                                                                       !    on the surface.
                                                                       ! 4: Heatflux from dissociation reactions of one gas
                                                                       !    species reacting at collision to another species
                                                                       !    on the surface.
                                                                       ! 5: additional heatflux e.g. surface coverage
                                                                       !    reconstruction or none of the above
                                                                       ! 5+iSpecies: Coverage of iSpecies
                                                                       !    adsorption%coverage added in updatesurfacevars
  REAL,ALLOCATABLE                      :: SurfModelReactCount(:,:,:,:)! 1-2*nReact,1-nSpecies: E-R + LHrecombination coefficient
                                                                       ! (2*nReact,nSpecies,p,q)
                                                                       ! doubled entries due to adsorb and desorb direction counter
  REAL,ALLOCATABLE                      :: Accomodation(:,:,:)         ! 1-nSpecies: Accomodation
                                                                       ! (nSpecies,p,q)
  !REAL, ALLOCATABLE                    :: Energy(:,:,:)               ! 1-3 E_tra (pre, wall, re),
  !                                                                    ! 4-6 E_rot (pre, wall, re),
  !                                                                    ! 7-9 E_vib (pre, wall, re)
  !REAL, ALLOCATABLE                    :: Force(:,:,:)                ! x, y, z direction
  !REAL, ALLOCATABLE                    :: Counter(:,:,:)              ! Wall-Collision counter
  REAL                                  :: PumpCapacity                !

  REAL,ALLOCATABLE                      :: ImpactEnergy(:,:,:,:)       ! 1-nSpecies: Particle impact energy (trans, rot, vib)
  REAL,ALLOCATABLE                      :: ImpactVector(:,:,:,:)       ! 1-nSpecies: Particle impact vector (x,y,z)
  REAL,ALLOCATABLE                      :: ImpactAngle(:,:,:)          ! 1-nSpecies: Particle impact angle (angle between particle
  REAL,ALLOCATABLE                      :: ImpactNumber(:,:,:)         ! 1-nSpecies: Number of particle impacts on surface
END TYPE
TYPE(tSampWall), ALLOCATABLE            :: SampWall(:)             ! Wall sample array (number of BC-Sides)

INTEGER                                 :: nPorousBC              ! Number of porous BCs
INTEGER                                 :: nPorousBCVars = 2      ! 1: Number of particles impinging the porous BC
                                                                  ! 2: Number of particles deleted at the porous BC
INTEGER                                 :: PorousBCSampIter       !
REAL, ALLOCATABLE                       :: PorousBCMacroVal(:,:,:)!

TYPE tPorousBC
  INTEGER                               :: BC                     ! Number of the reflective BC to be used as a porous BC
  REAL                                  :: Pressure               ! Pressure at the BC [Pa], user-given
  REAL                                  :: Temperature            ! Temperature at the BC [K], user-given
  REAL                                  :: NumberDensity          ! Calculated number density [1/m3]
  REAL                                  :: PumpingSpeed           ! Given/calculated pumping speed [m3/s]
  REAL                                  :: DeltaPumpingSpeedKp    ! Proportional factor for the pumping speed controller
  REAL                                  :: DeltaPumpingSpeedKi    ! Integral factor for the pumping speed controller
  CHARACTER(LEN=50)                     :: Region                 ! Form of the porous BC: 'circular'
  LOGICAL                               :: UsingRegion            ! Use only a smaller region on the BC as a porous BC (e.g. pump)
  INTEGER                               :: dir(3)                 ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                                  :: origin(2)              ! origin in orth. coordinates of polar system
  REAL                                  :: rmax                   ! max radius of to-be inserted particles
  REAL                                  :: rmin                   ! min radius of to-be inserted particles
  INTEGER                               :: SideNumber             ! Number of BC sides for the BC
  INTEGER, ALLOCATABLE                  :: SideList(:)            ! Mapping from porous BC side list to the BC side list
  REAL, ALLOCATABLE                     :: Sample(:,:)            ! Allocated with SideNumber and nPorousBCVars
  INTEGER, ALLOCATABLE                  :: RegionSideType(:)      ! 0: side is completely inside porous region
                                                                  ! 1: side is completely outside porous region
                                                                  ! 2: side is partially inside porous region
  REAL, ALLOCATABLE                     :: RemovalProbability(:)  ! Removal probability at the porous BC
  REAL, ALLOCATABLE                     :: PressureDifference(:)  ! Removal probability at the porous BC
  REAL, ALLOCATABLE                     :: PumpingSpeedSide(:)    ! Removal probability at the porous BC
  REAL                                  :: Output(1:5)            ! 1: Counter of impinged particles on the BC
                                                                  ! 2: Measured pumping speed [m3/s] through # of deleted particles
                                                                  ! 3: Pumping speed [m3/s] used to calculate the removal prob.
                                                                  ! 4: Removal probability [0-1]
                                                                  ! 5: Pressure at the BC normalized with the user-given pressure
END TYPE
TYPE(tPorousBC), ALLOCATABLE            :: PorousBC(:)            ! Container for the porous BC, allocated with nPorousBC

INTEGER, ALLOCATABLE                    :: MapBCtoPorousBC(:)     ! Mapping the porous BC to the BC (input: BC, output: porous BC)
INTEGER, ALLOCATABLE                    :: MapSurfSideToPorousSide(:) !

TYPE tSurfColl
  INTEGER                               :: NbrOfSpecies           ! Nbr. of Species to be counted for wall collisions (def. 0: all)
  LOGICAL,ALLOCATABLE                   :: SpeciesFlags(:)        ! Species counted for wall collisions (def.: all species=T)
  LOGICAL                               :: OnlySwaps              ! count only wall collisions being SpeciesSwaps (def. F)
  LOGICAL                               :: Only0Swaps             ! count only wall collisions being delete-SpeciesSwaps (def. F)
  LOGICAL                               :: Output                 ! Print sums of all counted wall collisions (def. F)
  LOGICAL                               :: AnalyzeSurfCollis      ! Output of collided/swaped particles
                                                                  ! during Sampling period? (def. F)
END TYPE
TYPE (tSurfColl)                        :: CalcSurfCollis

TYPE tAnalyzeSurfCollis
  INTEGER                               :: maxPartNumber          ! max. number of collided/swaped particles during Sampling
  REAL, ALLOCATABLE                     :: Data(:,:)              ! Output of collided/swaped particles during Sampling period
                                                                  ! (Species,Particles,Data(x,y,z,u,v,w)
  INTEGER, ALLOCATABLE                  :: Spec(:)                ! Species of Particle in Data-array
  INTEGER, ALLOCATABLE                  :: BCid(:)                ! ID of PartBC from crossing of Particle in Data-array
  INTEGER, ALLOCATABLE                  :: Number(:)              ! collided/swaped particles per Species during Sampling period
  !REAL, ALLOCATABLE                     :: Rate(:)                ! collided/swaped particles/s per Species during Sampling period
  INTEGER                               :: NumberOfBCs            ! Nbr of BC to be analyzed (def.: 1)
  INTEGER, ALLOCATABLE                  :: BCs(:)                 ! BCs to be analyzed (def.: 0 = all)

END TYPE tAnalyzeSurfCollis
TYPE(tAnalyzeSurfCollis)                :: AnalyzeSurfCollis

TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleAnodeBC           = 4      ! = 4 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleCathodeBC         = 5      ! = 5 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryAxis            = 11     ! = 10 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: AnalyzeBC               = 100    ! = 100 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)          ! Link part 1 for mapping PICLas BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)          ! Link part 2 for mapping PICLas BCs to Particle BC
!  INTEGER              , ALLOCATABLE     :: Map(:)                      ! Map from PICLas BCindex to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)              ! Map from PICLas BCindex to Particle BC (NOT TO TYPE!)
  !!INTEGER              , ALLOCATABLE     :: SideBCType(:)            ! list with boundary condition for each side
  REAL    , ALLOCATABLE                  :: MomentumACC(:)
  REAL    , ALLOCATABLE                  :: WallTemp(:)
  REAL    , ALLOCATABLE                  :: TransACC(:)
  REAL    , ALLOCATABLE                  :: VibACC(:)
  REAL    , ALLOCATABLE                  :: RotACC(:)
  REAL    , ALLOCATABLE                  :: ElecACC(:)
  REAL    , ALLOCATABLE                  :: WallVelo(:,:)
  REAL    , ALLOCATABLE                  :: Voltage(:), Voltage_CollectCharges(:)
  INTEGER , ALLOCATABLE                  :: NbrOfSpeciesSwaps(:)          !Number of Species to be changed at wall
  REAL    , ALLOCATABLE                  :: ProbOfSpeciesSwaps(:)         !Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE                  :: SpeciesSwaps(:,:,:)           !Species to be changed at wall (in, out), out=0: delete
  INTEGER , ALLOCATABLE                  :: SurfaceModel(:)               ! Model used for surface interaction
                                                                             ! 0 perfect/diffusive reflection
                                                                             ! 1 adsorption (Kisluik) / desorption (Polanyi Wigner)
                                                                             ! 2 Recombination coefficient (Laux model)
                                                                             ! 3 adsorption/desorption + chemical interaction
                                                                             !   (SMCR with UBI-QEP, TST)
                                                                             ! 4 TODO
                                                                             ! 5 SEE (secondary e- emission) by Levko2015
                                                                             ! 6 SEE (secondary e- emission) by Pagonakis2016
                                                                             !   (orignally from Harrower1956)
                                                                             ! 101 liquid condensation coeff = 1 + evaporation
                                                                             ! 102 liquid tsuruta model
  LOGICAL , ALLOCATABLE                  :: Reactive(:)                   ! flag defining if surface is treated reactively
  LOGICAL , ALLOCATABLE                  :: SolidState(:)                 ! flag defining if reflective BC is solid or liquid
  REAL    , ALLOCATABLE                  :: SolidPartDens(:)
  REAL    , ALLOCATABLE                  :: SolidMassIC(:)
  REAL    , ALLOCATABLE                  :: SolidAreaIncrease(:)
  INTEGER , ALLOCATABLE                  :: SolidStructure(:)             ! crystal structure of solid boundary (1:fcc100 2:fcc111)
  INTEGER , ALLOCATABLE                  :: SolidCrystalIndx(:)
  LOGICAL , ALLOCATABLE                  :: Adaptive(:)
  INTEGER , ALLOCATABLE                  :: AdaptiveType(:)
  INTEGER , ALLOCATABLE                  :: AdaptiveMacroRestartFileID(:)
  REAL    , ALLOCATABLE                  :: AdaptivePressure(:)
  REAL    , ALLOCATABLE                  :: AdaptiveTemp(:)
  LOGICAL , ALLOCATABLE                  :: UseForQCrit(:)                ! Use Boundary for Q-Criterion ?
  LOGICAL , ALLOCATABLE                  :: Resample(:)                   ! Resample Equilibirum Distribution with reflection
  LOGICAL , ALLOCATABLE                  :: Dielectric(:)                 ! Define if particle boundary [$] is a dielectric
!                                                                         ! interface, i.e. an interface between a dielectric and 
!                                                                         ! a non-dielectric or a between to different dielectrics 
!                                                                         ! [.TRUE.] or not [.FALSE.] (requires reflective BC)
!                                                                         ! (Default=FALSE.)
  LOGICAL , ALLOCATABLE                  :: BoundaryParticleOutput(:)     ! Save particle position, velocity and species to 
!                                                                         ! PartDataBoundary container for writing to .h5 later
END TYPE

INTEGER                                  :: nPartBound                       ! number of particle boundaries
INTEGER                                  :: nAdaptiveBC
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
LOGICAL              :: DoBoundaryParticleOutput   ! Flag set automatically if particles crossing specific 
!                                                  ! boundaries are to be saved to .h5 (position of intersection, 
!                                                  ! velocity, species, internal energies)
REAL, ALLOCATABLE    :: PartStateBoundary(:,:)     ! (1:9,1:NParts) 1st index: x,y,z,vx,vy,vz,MPF,time,impact angle
!                                                  !                2nd index: 1 to number of boundary-crossed particles
INTEGER, ALLOCATABLE :: PartStateBoundarySpec(:)   ! Species ID of boundary-crossed particles
INTEGER              :: PartStateBoundaryVecLength ! Number of boundary-crossed particles
!===================================================================================================================================

END MODULE MOD_Particle_Boundary_Vars
