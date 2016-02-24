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
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: BezierControlPoints3DElevated! Bezier basis control points of degree equal to NGeo
REAL,ALLOCATABLE,DIMENSION(:,:)         :: ElevationMatrix              ! array for binomial coefficients used for Bezier Elevation
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: SideSlabNormals              ! normal vectors of bounding slab box (Sides)
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SideSlabIntervals            ! intervalls beta1, beta2, beta3 (Sides)
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: ElemSlabNormals              ! normal vectors of bounding slab box (Elements)
REAL,ALLOCATABLE,DIMENSION(:,:)         :: ElemSlabIntervals            ! intervalls beta1, beta2, beta3 (Elements)
REAL,ALLOCATABLE,DIMENSION(:,:)         :: Vdm_Bezier,sVdm_Bezier       ! Vdm from/to Bezier Polynomial from BC representation
REAL,ALLOCATABLE,DIMENSION(:,:)         :: D_Bezier                     ! D-Matrix of Bezier Polynomial from BC representation
REAL,ALLOCATABLE,DIMENSION(:,:)         :: arrayNchooseK                ! array for binomial coefficients
REAL,ALLOCATABLE,DIMENSION(:,:)         :: FacNchooseK                  ! array for binomial coefficients times prefactor
INTEGER,ALLOCATABLE,DIMENSION(:)        :: SideType                     ! integer array with side type - planar - bilinear - curved
LOGICAL,ALLOCATABLE,DIMENSION(:)        :: BoundingBoxIsEmpty           ! logical if Side bounding box is empty
REAL,ALLOCATABLE,DIMENSION(:,:)         :: SideNormVec                  ! normal Vector of planar sides
REAL,ALLOCATABLE,DIMENSION(:)           :: SideDistance                 ! distance of planar base from origin 
INTEGER,ALLOCATABLE,DIMENSION(:)        :: gElemBCSides                 ! number of BC-Sides of element
REAL                                    :: BezierEpsilonBilinear        ! bi-linear tolerance for the bi-linear - planar decision
REAL                                    :: BezierHitEpsBi               ! epsilon tolerance for bi-linear faces
REAL                                    :: epsilontol                   ! epsilon for setting the tolerance
REAL                                    :: OneMinusEps                  ! 1 - eps: epsilontol
REAL                                    :: OnePlusEps                   ! 1 + eps: epsilontol for setting the boundary tolerance
REAL                                    :: MinusEps                     ! - eps: epsilontol
LOGICAL                                 :: ParticleSurfaceInitIsDone=.FALSE.
! settings for Bezier-Clipping and definition of maximal number of intersections
REAL                                    :: BezierNewtonAngle            ! switch for intersection with bezier newton algorithm
                                                                        ! smallest angle of impact of particle trajectory on face
REAL                                    :: BezierClipHit                ! value for clip hit
REAL                                    :: BezierClipTolerance          ! tolerance for root of bezier clipping
REAL                                    :: BezierSplitLimit             ! clip if remaining area after clip is > clipforce %
INTEGER                                 :: BezierClipMaxIntersec        ! maximal possible intersections for Bezier clipping
INTEGER                                 :: BezierClipMaxIter            ! maximal iterations per intersections
INTEGER                                 :: BezierElevation              ! elevate polynomial degree to NGeo+BezierElevation
REAL,ALLOCATABLE,DIMENSION(:)           :: locAlpha,locXi,locEta        ! position of trajectory-patch
REAL,ALLOCATABLE,DIMENSION(:,:)         :: XiArray,EtaArray             ! xi and eta history for computation of intersection
!LOGICAL                                 :: MultipleBCs                  ! allow for multiple BC during one tracking step
                                                                        ! only for do-ref-mapping required
INTEGER                                 :: BezierSampleN                ! equidistant sampling of bezier surface for emission
REAL,ALLOCATABLE,DIMENSION(:)           :: BezierSampleXi               ! ref coordinate for equidistant bezier surface sampling
LOGICAL                                 :: BezierSampleProjection       ! do a projection in the direction of an asigned vector
REAL,DIMENSION(3)                       :: BezierSampleProjectionVec    ! Projection vector
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: SurfMeshSubSideAreas         ! areas of of sub-sides of surface mesh
                                                                        ! (1:BezierSampleN,1:BezierSampleN,1:nBCSides)
REAL,ALLOCATABLE,DIMENSION(:)           :: SurfMeshSideAreas            ! areas of of sides of surface mesh (1:nBCSides)
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: SurfMeshProjSubSideAreas     ! projected areas of of sub-sides of surface mesh
REAL,ALLOCATABLE,DIMENSION(:)           :: SurfMeshProjSideAreas        ! projected areas of of sides of surface mesh
LOGICAL                                 :: BezierSampledAreasInitIsDone
#ifdef CODE_ANALYZE
REAL                                    :: rBoundingBoxChecks           ! number of bounding box checks
REAL(KIND=16)                           :: rTotalBBChecks               ! total number of bounding box checks
REAL                                    :: rPerformBezierClip           ! number of performed bezier clips
REAL                                    :: rPerformBezierNewton         ! number of performed bezier newton intersections
REAL(KIND=16)                           :: rTotalBezierClips            ! total number of performed bezier clips
REAL(KIND=16)                           :: rTotalBezierNewton           ! total number of performed bezier newton intersections
REAL,ALLOCATABLE,DIMENSION(:)           :: SideBoundingBoxVolume        ! Bounding Box volume
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

END MODULE MOD_Particle_Surfaces_Vars
