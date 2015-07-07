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
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: BezierControlPoints3D        ! Bezier basis control points of degree equal to NGeo
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: SlabNormals                  ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SlabIntervalls               ! intervalls beta1, beta2, beta3
REAL,ALLOCATABLE,DIMENSION(:,:)         :: Vdm_Bezier,sVdm_Bezier       ! Vdm from/to Bezier Polynomial from BC representation
REAL,ALLOCATABLE,DIMENSION(:,:)         :: arrayNchooseK                ! array for binomial coefficients
REAL,ALLOCATABLE,DIMENSION(:,:)         :: FacNchooseK                  ! array for binomial coefficients times prefactor
INTEGER,ALLOCATABLE,DIMENSION(:)        :: SideType                     ! integer array with side type - planar - bilinear - curved
LOGICAL,ALLOCATABLE,DIMENSION(:)        :: BoundingBoxIsEmpty           ! logical if Side bounding box is empty
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SideNormVec                  ! normal Vector of planar sides
REAL,ALLOCATABLE,DIMENSION(:)           :: SideDistance                 ! distance of planar base from origin 
INTEGER,ALLOCATABLE,DIMENSION(:)        :: gElemBCSides                 ! number of BC-Sides of element
REAL                                    :: epsilonbilinear              ! bi-linear tolerance for the bi-linear - planar decision
REAL                                    :: epsilonOne                   ! epsilone for setting the boundary tolerance
REAL                                    :: hitEpsBi
REAL                                    :: OneMepsilon
REAL                                    :: epsilontol                   ! epsilone for setting the tolerance
REAL                                    :: Mepsilontol               
REAL                                    :: ClipHit                      ! value for clip hit
LOGICAL                                 :: ParticleSurfaceInitIsDone=.FALSE.
! settings for Bezier-Clipping and definition of maximal number of intersections
REAL                                    :: ClipTolerance                 ! tolerance for root of bezier clipping
REAL                                    :: SplitLimit                    ! clip if remaining area after clip is > clipforce %
INTEGER                                 :: ClipMaxInter                  ! maximal possible intersections for Bezier clipping
INTEGER                                 :: ClipMaxIter                   ! maximal iterations per intersections
REAL,ALLOCATABLE,DIMENSION(:)           :: locAlpha,locXi,locEta         ! position of trajectory-patch
REAL,ALLOCATABLE,DIMENSION(:,:)         :: XiArray,EtaArray              ! xi and eta history for computation of intersection
!LOGICAL                                 :: MultipleBCs                   ! allow for multiple BC during one tracking step
                                                                         ! only for do-ref-mapping required
!===================================================================================================================================

END MODULE MOD_Particle_Surfaces_Vars
