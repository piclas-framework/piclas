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
INTEGER,ALLOCATABLE                     :: OffSetParticleBCSampling(:)   ! integer offset for particle boundary sampling            

TYPE tSurfaceMesh
  LOGICAL                               :: SurfSampOnProc                ! flag if reflective boundary condition is on proc
  INTEGER                               :: nSurfaceBCSides               ! Number of Sides on Surface (reflective)
  INTEGER,ALLOCATABLE                   :: SideIDToSurfaceID(:)          ! Mapping form the SideID to shorter side list
  REAL, ALLOCATABLE                     :: SurfaceArea(:,:,:)            ! Area of Surface 
END TYPE

TYPE (tSurfaceMesh)                     :: SurfMesh

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
