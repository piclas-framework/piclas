#include "boltzplatz.h"
MODULE MOD_Particle_Boundary_Vars
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: NSurfSample                   ! polynomial degree of particle BC sampling   
REAL,ALLOCATABLE                        :: XiEQ_SurfSample(:)            ! position of XiEQ_SurfSample
REAL                                    :: dXiEQ_SurfSample              ! deltaXi in [-1,1]
INTEGER                                 :: OffSetSurfSide                ! offset of local surf side
#ifdef MPI
INTEGER,ALLOCATABLE                     :: OffSetSurfSideMPI(:)          ! integer offset for particle boundary sampling            
#endif /*MPI*/

#ifdef MPI
TYPE tSurfaceSendList
  INTEGER                               :: NativeProcID
  INTEGER,ALLOCATABLE                   :: SendList(:)                   ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: RecvList(:)                   ! list containing surfsideid of sides to recv from proc
END TYPE
#endif /*MPI*/

TYPE tSurfaceCOMM
  LOGICAL                               :: MPIRoot                       ! if root of mpi communicator
  INTEGER                               :: MyRank                        ! local rank in new group
  INTEGER                               :: nProcs                        ! number of processes
#ifdef MPI
  INTEGER                               :: COMM                          ! communicator
  INTEGER                               :: nMPINeighbors                 ! number of processes to communicate with
  TYPE(tSurfaceSendList),ALLOCATABLE    :: MPINeighbor(:)                ! list containing all mpi neighbors
#endif /*MPI*/
END TYPE
TYPE (tSurfaceCOMM)                     :: SurfCOMM

TYPE tSurfaceMesh
  INTEGER                               :: SampSize                      ! integer of sampsize
  LOGICAL                               :: SurfOnProc                    ! flag if reflective boundary condition is on proc
  INTEGER                               :: nSides                        ! Number of Sides on Surface (reflective)
  INTEGER                               :: nTotalSides                   ! Number of Sides on Surface incl. HALO sides
  INTEGER                               :: nGlobalSides                  ! Global number of Sides on Surfaces (reflective)
  INTEGER,ALLOCATABLE                   :: SideIDToSurfID(:)             ! Mapping form the SideID to shorter side list
  REAL, ALLOCATABLE                     :: SurfaceArea(:,:,:)            ! Area of Surface 
END TYPE

TYPE (tSurfaceMesh)                     :: SurfMesh

TYPE tSampWall             ! DSMC sample for Wall                                             
  ! easier to communicate
  REAL,ALLOCATABLE                      :: State(:,:,:)                ! 1-3   E_tra (pre, wall, re),
                                                                       ! 4-6   E_rot (pre, wall, re),
                                                                       ! 7-9   E_vib (pre, wall, re)
                                                                       ! 10-12 x, y, z direction
                                                                       ! 13-12+nSpecies Wall-Collision counter
  !REAL, ALLOCATABLE                    :: Energy(:,:,:)               ! 1-3 E_tra (pre, wall, re),
  !                                                                    ! 4-6 E_rot (pre, wall, re),
  !                                                                    ! 7-9 E_vib (pre, wall, re)
  !REAL, ALLOCATABLE                    :: Force(:,:,:)                ! x, y, z direction
  !REAL, ALLOCATABLE                    :: Counter(:,:,:)              ! Wall-Collision counter
END TYPE
TYPE(tSampWall), ALLOCATABLE           :: SampWall(:)             ! Wall sample array (number of BC-Sides)


TYPE tMacroSurfaceVal                                       ! DSMC sample for Wall    
  REAL                           :: Heatflux                ! 
  REAL                           :: Force(3)                ! x, y, z direction
  REAL, ALLOCATABLE              :: Counter(:)              ! Wall-Collision counter of all Species
  REAL                           :: CounterOut              ! Wall-Collision counter for Output
END TYPE

TYPE(tMacroSurfaceVal), ALLOCATABLE     :: MacroSurfaceVal(:) ! Wall sample array (number of BC-Sides)



TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleAnodeBC           = 4      ! = 4 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleCathodeBC         = 5      ! = 5 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)          ! Link part 1 for mapping Boltzplatz BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)          ! Link part 2 for mapping Boltzplatz BCs to Particle BC
!  INTEGER              , ALLOCATABLE     :: Map(:)                      ! Map from Boltzplatz BCindex to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)              ! Map from Boltzplatz BCindex to Particle BC (NOT TO TYPE!)
  !!INTEGER              , ALLOCATABLE     :: SideBCType(:)            ! list with boundary condition for each side
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




!===================================================================================================================================

END MODULE MOD_Particle_Boundary_Vars
