#include "boltzplatz.h"

MODULE MOD_Particle_Surfaces_Vars
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
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: BiLinearCoeff              ! contains the bi-linear coefficients for each side
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: SuperSampledNodes          !  
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)   :: SuperSampledBiLinearCoeff  !
INTEGER                                 :: NPartCurved                !
INTEGER                                 :: nTriangles
LOGICAL                                 :: DoPartCurved=.FALSE.       !
REAL,ALLOCATABLE,DIMENSION(:,:)         :: Vdm_CLNGeo_EquiNPartCurved
LOGICAL,ALLOCATABLE,DIMENSION(:)        :: SideIsPlanar               ! logical error if side is planar, instead of bi-linear
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SideNormVec                ! normal Vector of planar sides
REAL,ALLOCATABLE,DIMENSION(:)           :: SideDistance               ! distance of planar base from origin 
INTEGER,ALLOCATABLE,DIMENSION(:)        :: gElemBCSides               ! number of BC-Sides of element
INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: neighborElemID,neighborlocSideID
REAL                                    :: epsilonbilinear            ! bi-linear tolerance for the bi-linear - planar decision
REAL                                    :: epsilonOne                 ! epsilone for setting the boundary tolerance
REAL                                    :: epsilontol                 ! epsilone for setting the tolerance
LOGICAL                                 :: ParticleSurfaceInitIsDone=.FALSE.
!===================================================================================================================================

END MODULE MOD_Particle_Surfaces_Vars
