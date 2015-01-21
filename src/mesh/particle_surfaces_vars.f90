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
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: BiLinearCoeff                ! contains the bi-linear coefficients for each side
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: SuperSampledNodes            !  
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: BezierControlPoints3D        ! Bezier basis control points of degree equal to NGeo
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: SlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SlabIntervalls               ! intervalls beta1, beta2, beta3
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)   :: SuperSampledBiLinearCoeff    !
INTEGER                                 :: NPartCurved                  !
INTEGER                                 :: nTriangles,nQuads
LOGICAL                                 :: DoPartCurved=.FALSE.         !
REAL,ALLOCATABLE,DIMENSION(:,:)         :: Vdm_CLNGeo_EquiNPartCurved
REAL,ALLOCATABLE,DIMENSION(:,:)         :: Vdm_Bezier,sVdm_Bezier       ! 
REAL,ALLOCATABLE,DIMENSION(:,:)         :: arrayNchooseK                ! array for binomial coefficients
REAL,ALLOCATABLE,DIMENSION(:,:)         :: FacNchooseK                  ! array for binomial coefficients times prefactor
INTEGER,ALLOCATABLE,DIMENSION(:)        :: SideType                     ! integer array with side type - planar - bilinear - curved
LOGICAL,ALLOCATABLE,DIMENSION(:)        :: BoundingBoxIsEmpty
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SideNormVec                  ! normal Vector of planar sides
REAL,ALLOCATABLE,DIMENSION(:)           :: SideDistance                 ! distance of planar base from origin 
INTEGER,ALLOCATABLE,DIMENSION(:)        :: gElemBCSides                 ! number of BC-Sides of element
INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: neighborElemID,neighborlocSideID
REAL                                    :: epsilonbilinear              ! bi-linear tolerance for the bi-linear - planar decision
REAL                                    :: epsilonOne                   ! epsilone for setting the boundary tolerance
REAL                                    :: OneMepsilon
REAL                                    :: epsilontol                   ! epsilone for setting the tolerance
REAL                                    :: Mepsilontol               
LOGICAL                                 :: ParticleSurfaceInitIsDone=.FALSE.
! settings for Bezier-Clipping and definition of maximal number of intersections
REAL                                    :: ClipTolerance
REAL                                    :: SplitLimit                    ! clip if remaining area after clip is > clipforce %
INTEGER                                 :: ClipMaxInter,ClipMaxIter
REAL,ALLOCATABLE,DIMENSION(:)           :: locAlpha,locXi,locEta
REAL,ALLOCATABLE,DIMENSION(:,:)         :: XiArray,EtaArray
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: XiEtaZetaBasis
REAL,ALLOCATABLE,DIMENSION(:,:)         :: slenXiEtaZetaBasis
REAL,ALLOCATABLE,DIMENSION(:,:)         :: ElemBaryNGeo
INTEGER                                 :: MappingGuess
!===================================================================================================================================

END MODULE MOD_Particle_Surfaces_Vars
